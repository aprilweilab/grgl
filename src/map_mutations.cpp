/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2024 April Wei
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "grgl/map_mutations.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/node_data.h"
#include "grgl/visitor.h"
#include "util.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <sys/types.h>
#include <unordered_map>
#include <utility>
#include <vector>
// When enabled: garbage collects unneeded sample sets
#define CLEANUP_SAMPLE_SETS_MAPPING 1

// Every 10% we'll emit some stats about current size of the GRG.
#define EMIT_STATS_AT_PERCENT (10)
// Every 5% we compact the GRG edge lists.
#define COMPACT_EDGES_AT_PERCENT (5)
#define ONE_HUNDRED_PERCENT      (100)

// Histogram size for statistics
#define STATS_HIST_SIZE (200)

namespace grgl {

using NodeSamples = std::tuple<NodeID, size_t>;
using CandidateList = std::vector<NodeSamples>;

static bool cmpNodeSamples(const NodeSamples& ns1, const NodeSamples& ns2) {
    const size_t samplesCount1 = std::get<1>(ns1);
    const size_t& samplesCount2 = std::get<1>(ns2);
    return samplesCount1 < samplesCount2;
}

// What percentage of total samples must be covered before we switch to dense bitvecs
#define USE_DENSE_COVERAGE_PERCENT (0.001f)

class BitVecCoverageSet;
class SparseCoverageSet;

class BatchMembership {
public:
    BatchMembership() = default;

    void setOnes(size_t mutBatchSize) {
        m_numBits = mutBatchSize;
        const size_t wordCount = (mutBatchSize + 63) / 64;
        m_words.assign(wordCount, ~0ULL);

        const size_t remainderBits = mutBatchSize % 64;
        if (remainderBits != 0 && !m_words.empty()) {
            m_words.back() &= ((1ULL << remainderBits) - 1ULL);
        }
    }

    void clear() {
        for (uint64_t& word : m_words) {
            word = 0ULL;
        }
    }

    bool any() const {
        return std::any_of(m_words.begin(), m_words.end(), [](uint64_t word) { return word != 0ULL; });
    }

    void andAssign(const BatchMembership& other) {
        for (size_t wordIndex = 0; wordIndex < m_words.size(); ++wordIndex) {
            m_words[wordIndex] &= other.m_words[wordIndex];
        }
    }

    void setBit(size_t bitIndex) {
        const size_t wordIndex = bitIndex / 64;
        const size_t bitOffset = bitIndex % 64;
        m_words[wordIndex] |= (1ULL << bitOffset);
    }

    class Iterator {
    public:
        Iterator(const BatchMembership* left, const BatchMembership* right, bool atEnd)
            : m_left(left),
              m_right(right),
              m_wordIndex(atEnd ? left->m_words.size() : 0) {
            if (!atEnd) {
                ++(*this);
            }
        }

        size_t operator*() const {
            const int bitOffset = __builtin_ctzll(m_currentWord);
            return (m_wordIndex * 64) + static_cast<size_t>(bitOffset);
        }

        Iterator& operator++() {
            // Consume current bit if present
            if (m_currentWord != 0ULL) {
                m_currentWord &= (m_currentWord - 1ULL);
            }

            // Find next nonzero word
            while (m_currentWord == 0ULL && m_wordIndex < m_left->m_words.size()) {
                uint64_t word = m_left->m_words[m_wordIndex];
                if (m_right != nullptr) {
                    word ^= m_right->m_words[m_wordIndex];
                }

                m_currentWord = word;
                if (m_currentWord != 0ULL) {
                    return *this;
                }
                ++m_wordIndex;
            }

            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return m_wordIndex != other.m_wordIndex || m_currentWord != other.m_currentWord;
        }

    private:
        const BatchMembership* m_left;
        const BatchMembership* m_right; // nullptr => bits(), non-null => XOR
        size_t m_wordIndex;
        uint64_t m_currentWord{};
    };

    Iterator begin() const { return {this, nullptr, false}; }
    Iterator end() const { return {this, nullptr, true}; }

    Iterator beginDiff(const BatchMembership& other) const { return {this, &other, false}; }
    Iterator endDiff(const BatchMembership& other) const { return {this, &other, true}; }

private:
    std::vector<uint64_t> m_words;
    size_t m_numBits = 0;
};

class SampleCoverageSet {
public:
    /// the number of total samples
    virtual size_t totalSamples() const = 0;

    /// number of elements in the set, INCLUDING DUPLICATES from below nodes
    virtual size_t numSamplesCovered() const = 0;

    virtual bool empty() const = 0;
    virtual void addElem(NodeID sample) = 0;
    void mergeSampleCoverage(const SampleCoverageSet& other) {
        m_batchMembership.andAssign(other.m_batchMembership);
        mergeSampleCoverageImpl(other);
    }

    void mergeSampleCoverage(const SampleCoverageSet& other, BitVecCoverageSet& seenSet, size_t& coalescences) {
        m_batchMembership.andAssign(other.m_batchMembership);
        mergeSampleCoverageImpl(other, seenSet, coalescences);
    }

    void mergeSampleCoverage(const SampleCoverageSet& other,
                             std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                             NodeID childNode,
                             size_t& coalescences) {
        m_batchMembership.andAssign(other.m_batchMembership);
        mergeSampleCoverageImpl(other, individualToChild, childNode, coalescences);
    }

    const BatchMembership& getBatchMembership() const { return m_batchMembership; }
    BatchMembership& getBatchMembershipMut() { return m_batchMembership; }

    SampleCoverageSet() = default;
    SampleCoverageSet(SampleCoverageSet&&) = delete;
    SampleCoverageSet(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(SampleCoverageSet&&) = default;
    virtual ~SampleCoverageSet() = default;

protected:
    virtual void mergeSampleCoverageImpl(const SampleCoverageSet& other) = 0;
    virtual void
    mergeSampleCoverageImpl(const SampleCoverageSet& other, BitVecCoverageSet& seenSet, size_t& coalescences) {
        mergeSampleCoverageImpl(other);
    }
    virtual void mergeSampleCoverageImpl(const SampleCoverageSet& other,
                                         std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                                         NodeID childNode,
                                         size_t& coalescences) {
        mergeSampleCoverageImpl(other);
    }

private:
    BatchMembership m_batchMembership;
};

class SparseCoverageSet : public SampleCoverageSet {
public:
    explicit SparseCoverageSet(size_t m_totalSampleCount, size_t reserve_amount = 0)
        : m_totalSampleCount(m_totalSampleCount) {
        if (reserve_amount > 0) {
            m_sampleIdxs.reserve(reserve_amount);
        }
    }

    size_t totalSamples() const override { return m_totalSampleCount; }
    size_t numSamplesCovered() const override { return m_sampleIdxs.size(); }
    bool empty() const override { return m_sampleIdxs.empty(); }
    void addElem(NodeID sample) override { m_sampleIdxs.push_back(sample); }

protected:
    void mergeSampleCoverageImpl(const SampleCoverageSet& other) override {
        const auto* otherSparse = dynamic_cast<const SparseCoverageSet*>(&other);
        if (otherSparse == nullptr || otherSparse->m_sampleIdxs.empty()) {
            return;
        }

        appendSamples(*otherSparse);
    }

    void mergeSampleCoverageImpl(const SampleCoverageSet& other,
                                 std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                                 NodeID childNode,
                                 size_t& coalescences) override {
        const auto* otherSparse = dynamic_cast<const SparseCoverageSet*>(&other);
        if (otherSparse == nullptr || otherSparse->m_sampleIdxs.empty()) {
            mergeSampleCoverageImpl(other);
            return;
        }

        for (const NodeID sampleId : otherSparse->m_sampleIdxs) {
            m_sampleIdxs.emplace_back(sampleId);
            const auto insertPair = individualToChild.emplace(sampleId / 2, childNode);
            if (!insertPair.second && insertPair.first->second != childNode) {
                coalescences++;
            }
        }
    }

    friend class BitVecCoverageSet;

private:
    void appendSamples(const SparseCoverageSet& other) {
        m_sampleIdxs.insert(m_sampleIdxs.end(), other.m_sampleIdxs.begin(), other.m_sampleIdxs.end());
    }

    NodeIDList m_sampleIdxs;
    size_t m_totalSampleCount;
};

/// A BitVec-backed set for coverage
class BitVecCoverageSet : public SampleCoverageSet {
public:
    friend class CoalescenceVisitor;
    explicit BitVecCoverageSet(size_t totalSamples)
        : m_totalSampleCount(totalSamples),
          m_elems((totalSamples + m_blockSize - 1) / m_blockSize) {}

    explicit BitVecCoverageSet(const SparseCoverageSet& sparseCoverageSet)
        : BitVecCoverageSet(sparseCoverageSet.totalSamples()) {
        mergeSampleCoverage(sparseCoverageSet);
    }

    size_t totalSamples() const override { return m_totalSampleCount; }

    size_t numSamplesCovered() const override { return m_numSampleCoverage; }

    bool empty() const override { return m_numSamplesNonOverlapping == 0; }

    void addElem(NodeID sample) override {
        size_t vecIdx = sample / m_blockSize;
        size_t offset = sample % m_blockSize;
        m_numSampleCoverage += 1;
        // if this bit was zero before, increment by one
        m_numSamplesNonOverlapping += ~(m_elems[vecIdx] >> offset) & 1ULL;
        m_elems[vecIdx] |= (1ULL << offset);
    }

    bool contains(const NodeID sample) const {
        size_t vecIdx = sample / m_blockSize;
        size_t offset = sample % m_blockSize;
        return ((m_elems[vecIdx] >> offset) & 1ULL) != 0;
    }

    // original api does not dedup sets when computing sizes.
    size_t numSamplesNonOverlapping() const { return m_numSamplesNonOverlapping; }

    bool overlapsWith(const SampleCoverageSet& other) const {
        release_assert(totalSamples() == other.totalSamples());
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            if (numSamplesNonOverlapping() + other_v->numSamplesNonOverlapping() > totalSamples()) {
                return true;
            }
            const uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();
            for (int i = 0; i < len; i++) {
                if ((data[i] & otherData[i]) != 0) {
                    return true;
                }
            }
        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (const uint32_t sample : other_v->m_sampleIdxs) {
                size_t vecIdx = sample / m_blockSize;
                size_t offset = sample % m_blockSize;
                if ((m_elems[vecIdx] & (1ULL << offset)) != 0) {
                    return true;
                }
            }
        }
        return false;
    }

protected:
    void mergeSampleCoverageImpl(const SampleCoverageSet& other) override {
        release_assert(totalSamples() == other.totalSamples());
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            m_numSampleCoverage += other.numSamplesCovered();
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();
            for (size_t i = 0; i < len; ++i) {
                m_numSamplesNonOverlapping += __builtin_popcountll(otherData[i] & ~data[i]);
                data[i] |= otherData[i];
            }

        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (const uint32_t sample : other_v->m_sampleIdxs) {
                addElem(sample);
            }
        }
    }

    void
    mergeSampleCoverageImpl(const SampleCoverageSet& other, BitVecCoverageSet& seenSet, size_t& coalescences) override {
        release_assert(totalSamples() == other.totalSamples());
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            m_numSampleCoverage += other.numSamplesCovered();
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            uint64_t* __restrict seen_d = seenSet.m_elems.data();
            const size_t len = m_elems.size();
            constexpr uint64_t EVEN = 0x5555555555555555ULL; // all even bits
            for (size_t i = 0; i < len; ++i) {
                const uint64_t word = otherData[i];
                m_numSamplesNonOverlapping += __builtin_popcountll(word & ~data[i]);
                data[i] |= word;
                const uint64_t coalesced =
                    (word | (word >> 1ULL)) & EVEN; // pack two consecutive haplotype bits into one individual slot
                const uint64_t new_coals = coalesced & seen_d[i];
                coalescences += __builtin_popcountll(new_coals);
                seen_d[i] |= coalesced;
            }

        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            auto& seenWords = seenSet.m_elems;
            for (const uint32_t sample : other_v->m_sampleIdxs) {
                addElem(sample);
                // Collapse the two haplotypes for an individual down to their even slot in the bitset.
                const size_t base = (static_cast<size_t>(sample) / 2) * 2; // even index per individual
                const size_t vecIdx = base / m_blockSize;
                const size_t offset = base % m_blockSize;
                const uint64_t mask = (uint64_t{1} << offset);

                if ((seenWords[vecIdx] & mask) != 0) {
                    coalescences++;
                }
                seenWords[vecIdx] |= mask;
            }
        } else {
            mergeSampleCoverageImpl(other);
        }
    }

    void mergeSampleCoverageImpl(const SampleCoverageSet& other,
                                 std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                                 NodeID childNode,
                                 size_t& coalescences) override {
        mergeSampleCoverageImpl(other);
    }

private:
    constexpr static size_t m_blockSize = 64;
    std::vector<uint64_t> m_elems;
    size_t m_totalSampleCount;
    size_t m_numSampleCoverage{};
    size_t m_numSamplesNonOverlapping{};
};

class MutationBatch {
    friend class TopoCandidateCollectorVisitor;

public:
    explicit MutationBatch(NodeIDSizeT totalSampleCount, size_t batchBitCount)
        : m_batchBitCount(batchBitCount),
          m_sampleMembership(totalSampleCount) {
        for (auto& membership : m_sampleMembership) {
            membership.setOnes(m_batchBitCount);
            membership.clear();
        }
    }

    void addMutation(const Mutation& mutation, NodeIDList mutSamples) {
        size_t mutIdx = m_mutations.size();
        m_mutations.emplace_back(mutation, m_seedList.size(), mutSamples.size());
        m_seedList.insert(m_seedList.end(), mutSamples.begin(), mutSamples.end());
        m_sampleSets.push_back(mutSamples);
        for (NodeID sample : mutSamples) {
            m_sampleMembership[sample].setBit(mutIdx);
        }
    }

    const NodeIDList& seedList() const { return m_seedList; }
    const Mutation& getMutation(size_t idx) const { return m_mutations[idx].m_mutation; }
    const NodeIDList& sampleSet(size_t idx) const { return m_sampleSets[idx]; }
    size_t numMutations() const { return m_mutations.size(); }

    void clear() {
        m_sampleSets.clear();
        m_mutations.clear();
        m_seedList.clear();
        for (auto& membership : m_sampleMembership) {
            membership.clear();
        }
    }

private:
    struct MutSpan {
        MutSpan(Mutation mutation, size_t start, size_t length)
            : m_mutation(std::move(mutation)),
              m_start(start),
              m_length(length) {}
        Mutation m_mutation;
        size_t m_start;
        size_t m_length;
    };
    NodeIDList m_seedList;
    std::vector<MutSpan> m_mutations;
    std::vector<NodeIDList> m_sampleSets;
    size_t m_batchBitCount;
    std::vector<BatchMembership> m_sampleMembership;
};

class TopoCandidateCollectorVisitor : public grgl::ParallelGRGVisitor {
public:
    explicit TopoCandidateCollectorVisitor(const MutationBatch& mutBatch)
        : m_batchBitCount(mutBatch.m_batchBitCount),
          m_sampleBatchMembership(mutBatch.m_sampleMembership),
          m_collectedNodes(mutBatch.m_batchBitCount) {}

    void parallelVisit(const grgl::GRGPtr& grg,
                       const grgl::NodeIDList& nodes,
                       std::vector<bool>& results,
                       grgl::TraversalDirection direction,
                       size_t numThreads) override {
        release_assert(direction == TraversalDirection::DIRECTION_UP);
#if CLEANUP_SAMPLE_SETS_MAPPING
        if (m_refCounts.empty()) {
            m_refCounts.resize(grg->numNodes());
        }
#endif
        if (m_sampleCoverage.empty()) {
            m_sampleCoverage.resize(grg->numNodes());
        }
#pragma omp parallel for schedule(guided) num_threads(numThreads) default(none) shared(nodes, grg, results)
        for (int i = 0; i < nodes.size(); ++i) {
            results[i] = safeVisit(grg, nodes[i]);
        }

#if CLEANUP_SAMPLE_SETS_MAPPING
        for (NodeID nodeId : nodes) {
            for (const auto& childId : grg->getDownEdges(nodeId)) {
                // Skip children that aren't part of our search.
                if (m_refCounts[childId] == 0) {
                    continue;
                }
                if (--m_refCounts[childId] == 0) {
                    m_sampleCoverage[childId].reset();
                }
            }
        }
#endif
    }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass = grgl::DfsPass::DFS_PASS_NONE) override {
        release_assert(direction == TraversalDirection::DIRECTION_UP);
        release_assert(dfsPass == DfsPass::DFS_PASS_NONE);
        release_assert(grg->hasUpEdges());

#if CLEANUP_SAMPLE_SETS_MAPPING
        if (m_refCounts.empty()) {
            m_refCounts.resize(grg->numNodes());
        }
#endif

        if (m_sampleCoverage.empty()) {
            m_sampleCoverage.resize(grg->numNodes());
        }

        bool keepGoing = safeVisit(grg, nodeId);

#if CLEANUP_SAMPLE_SETS_MAPPING
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            // Skip children that aren't part of our search.
            if (m_refCounts[childId] == 0) {
                continue;
            }
            if (--m_refCounts[childId] == 0) {
                m_sampleCoverage[childId].reset();
            }
        }
#endif
        return keepGoing;
    }

    const std::unique_ptr<SampleCoverageSet>& getSamplesForCandidate(NodeID candidateId) const {
        const std::unique_ptr<SampleCoverageSet>& value = m_sampleCoverage[candidateId];
        release_assert(value);
        return value;
    }

    std::vector<std::unique_ptr<SampleCoverageSet>> m_sampleCoverage;
    std::vector<CandidateList> m_collectedNodes;

private:
    size_t m_batchBitCount;
    const std::vector<BatchMembership>& m_sampleBatchMembership;
    // shared 'thread-safe' logic
    bool safeVisit(const grgl::GRGPtr& grg, grgl::NodeID nodeId);
    std::mutex m_collectedNodesMutex;
    // These are the _total_ samples beneath each node (not restricted to current samples being searched)
#if CLEANUP_SAMPLE_SETS_MAPPING
    std::vector<NodeIDSizeT> m_refCounts;
#endif
};

bool TopoCandidateCollectorVisitor::safeVisit(const grgl::GRGPtr& grg, const grgl::NodeID nodeId) {
    release_assert(grg->hasUpEdges());
    // Note: Any modification to current node or children is safe, as we never process parent and child at the same
    // time
    const bool isRoot = grg->numUpEdges(nodeId) == 0;
    const bool isSample = grg->isSample(nodeId);
    const size_t ploidy = grg->getPloidy();
    const auto numCoals = grg->getNumIndividualCoals(nodeId);
    const bool computeCoals = !isSample && (ploidy == 2) && (COAL_COUNT_NOT_SET == numCoals);

    size_t individualCoalCount = 0;
    NodeIDList candidateNodes;

    size_t sampleCount = 0;
    for (const auto& childId : grg->getDownEdges(nodeId)) {
        const auto& childSample = m_sampleCoverage[childId];
        if (childSample) {
            sampleCount += childSample->numSamplesCovered();
        }
    }

    std::unique_ptr<SampleCoverageSet> samplesBeneathSet;
    std::unique_ptr<BitVecCoverageSet> coalTrackerSet;
    std::unordered_map<NodeIDSizeT, NodeIDSizeT> coalTrackerMap;
    double coveragePercent = (double)sampleCount / (double)grg->numSamples();
    bool dense = coveragePercent > USE_DENSE_COVERAGE_PERCENT;
    if (dense) {
        samplesBeneathSet = std::unique_ptr<SampleCoverageSet>(new BitVecCoverageSet(grg->numSamples()));
        coalTrackerSet = std::unique_ptr<BitVecCoverageSet>(new BitVecCoverageSet(grg->numSamples()));
    } else {
        size_t sparse_reserve_amount = sampleCount;
        samplesBeneathSet =
            std::unique_ptr<SampleCoverageSet>(new SparseCoverageSet(grg->numSamples(), sparse_reserve_amount));
    }

    samplesBeneathSet->getBatchMembershipMut().setOnes(m_batchBitCount);

    if (isSample) {
        samplesBeneathSet->addElem(nodeId);
        samplesBeneathSet->getBatchMembershipMut() = m_sampleBatchMembership[nodeId];
    }

#if CLEANUP_SAMPLE_SETS_MAPPING
    m_refCounts[nodeId] = grg->numUpEdges(nodeId);
#endif
    for (const auto& childId : grg->getDownEdges(nodeId)) {
        const auto& childSample = m_sampleCoverage[childId];
        if (childSample) {
            if (childSample->numSamplesCovered() > 1) {
                candidateNodes.emplace_back(childId);
            }
            if (computeCoals) {
                if (dense) {
                    samplesBeneathSet->mergeSampleCoverage(*childSample, *coalTrackerSet, individualCoalCount);
                } else {
                    samplesBeneathSet->mergeSampleCoverage(*childSample, coalTrackerMap, childId, individualCoalCount);
                }
            } else {
                samplesBeneathSet->mergeSampleCoverage(*childSample);
            }
        } else {
            samplesBeneathSet->getBatchMembershipMut().clear();
        }
    }

    // We can only record coalescence counts if there are no samples missing.
    if (samplesBeneathSet->getBatchMembership().any() && computeCoals) {
        grg->setNumIndividualCoals(nodeId, individualCoalCount);
    }
    // If we've reached the root of the graph or have missing samples beneath us, we need to stop the search
    // and emit candidate nodes to map the mutation to.
    const BatchMembership& membership = samplesBeneathSet->getBatchMembership();
    const bool keepGoing = !isRoot && membership.any();

    if (isRoot) {
        // Root can only be a candidate when it fully covers the batch for at least one bit.
        if (membership.any()) {
            std::lock_guard<std::mutex> lock(m_collectedNodesMutex);
            for (auto it = membership.begin(); it != membership.end(); ++it) {
                const size_t mutationIndex = *it;
                m_collectedNodes[mutationIndex].emplace_back(nodeId, samplesBeneathSet->numSamplesCovered());
            }
#if CLEANUP_SAMPLE_SETS_MAPPING
            // Prevent candidates from having their samplesets garbage collected.
            m_refCounts[nodeId] = MAX_GRG_NODES + 1;
#endif
            // hack to allow concurrent writing
            m_sampleCoverage.at(nodeId).swap(samplesBeneathSet);
        }

        // Emit children for bits where they have coverage that the root does not.
        for (const auto& candidate : candidateNodes) {
            const BatchMembership& childMembership = m_sampleCoverage.at(candidate)->getBatchMembership();
            auto bitIterator = childMembership.beginDiff(membership);
            auto bitEnd = childMembership.endDiff(membership);
            if (bitIterator != bitEnd) {
                std::lock_guard<std::mutex> lock(m_collectedNodesMutex);
                for (; bitIterator != bitEnd; ++bitIterator) {
                    const size_t mutationIndex = *bitIterator;
                    m_collectedNodes[mutationIndex].emplace_back(candidate,
                                                                 m_sampleCoverage.at(candidate)->numSamplesCovered());
                }
#if CLEANUP_SAMPLE_SETS_MAPPING
                // Prevent candidates from having their samplesets garbage collected.
                m_refCounts[candidate] = MAX_GRG_NODES + 1;
#endif
            }
        }

        // No parent above root.
        return false;
    }

    for (const auto& candidate : candidateNodes) {
        const BatchMembership& childMembership = m_sampleCoverage.at(candidate)->getBatchMembership();
        auto bitIterator = childMembership.beginDiff(membership);
        auto bitEnd = childMembership.endDiff(membership);
        if (bitIterator != bitEnd) {
            std::lock_guard<std::mutex> lock(m_collectedNodesMutex);
            for (; bitIterator != bitEnd; ++bitIterator) {
                const size_t mutationIndex = *bitIterator;
                m_collectedNodes[mutationIndex].emplace_back(candidate,
                                                             m_sampleCoverage.at(candidate)->numSamplesCovered());
            }
#if CLEANUP_SAMPLE_SETS_MAPPING
            // Prevent candidates from having their samplesets garbage collected.
            m_refCounts[candidate] = MAX_GRG_NODES + 1;
#endif
        }
    }
    if (keepGoing) {
        m_sampleCoverage.at(nodeId).swap(samplesBeneathSet);
    }

    return keepGoing;
}

// Tracking individual coalescence is a bit spread out, but I think it is the most efficient way to do it.
// 1. Above, when searching for candidate nodes in the existing hierarchy, any node that does not have its
//    coalescence info computed will be computed and stored.
// 2. Below, when we connect two candidate nodes together we will look for any coalescences between them.
// 3. Below, when we have "left-over" samples that did not have any candidate nodes, we will check for any
//    coalescence within the left-over samples, or between the left-over samples and nodes that have already
//    been covered.
static NodeIDList processCandidateSet(const MutableGRGPtr& grg,
                                      const Mutation& newMutation,
                                      const NodeIDList& mutSamples,
                                      const TopoCandidateCollectorVisitor& collector,
                                      CandidateList& candidates,
                                      MutationMappingStats& stats,
                                      const NodeID shapeNodeIdMax,
                                      std::pair<BpPosition, NodeID>& currentMissing) {
    std::sort(candidates.begin(), candidates.end());
    auto endOfUnique = std::unique(candidates.begin(), candidates.end());
    candidates.erase(endOfUnique, candidates.end());
    std::sort(candidates.begin(), candidates.end(), cmpNodeSamples);

    // The missingness node associated with this mutation. This relies on the fact that
    // missing data is always emitted BEFORE other data for the same site, which is a
    // property of the MutationIterator.
    const NodeID missingnessNode =
        (currentMissing.first == newMutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;

    const int ploidy = grg->getPloidy();
    if (candidates.empty()) {
        stats.mutationsWithNoCandidates++;
    } else {
        // Exact match scenario. Return early.
        const auto& candidate = candidates.back();
        const size_t candidateSetSize = std::get<1>(candidate);
        if (candidateSetSize == mutSamples.size()) {
            const auto candidateId = std::get<0>(candidate);
            stats.reusedExactly++;
            if (!newMutation.isMissing()) {
                grg->addMutation(newMutation, candidateId, missingnessNode);
            } else {
                currentMissing = {newMutation.getPosition(), candidateId};
            }
            if (candidateId >= shapeNodeIdMax) {
                stats.reusedMutNodes++;
            }
            return {};
        }
    }

    size_t individualCoalCount = 0;
    BitVecCoverageSet coalTrackingSet{grg->numSamples()};
    BitVecCoverageSet coverageSet{grg->numSamples()};

    const NodeID mutNodeId = grg->makeNode(1, true);
    if (!newMutation.isMissing()) {
        grg->addMutation(newMutation, mutNodeId, missingnessNode);
    } else {
        currentMissing = {newMutation.getPosition(), mutNodeId};
    }
    NodeIDList addedNodes;
    const size_t numMutSamples = mutSamples.size();
    while (!candidates.empty() && coverageSet.numSamplesNonOverlapping() < numMutSamples) {
        const auto& candidate = candidates.back();
        const auto candidateId = std::get<0>(candidate);
        const std::unique_ptr<SampleCoverageSet>& candidateSet = collector.getSamplesForCandidate(candidateId);
        release_assert(!candidateSet->empty());
        // Different candidates may cover different subsets of the sample set that
        // we are currently trying to cover. Those sample sets MUST be non-overlapping
        // or we will introduce a diamond into the graph:
        //  m-->n1-->s0
        //  m-->n2-->s0
        // However, there is no guarantee that there does not exist nodes (n1, n2)
        // that both point to a sample (or samples) that we care about, so we have to
        // track that here. We do that by only considering candidates that have no overlap
        // with our already-covered set.
        if (!coverageSet.overlapsWith(*candidateSet)) {
            // Mark all the sample nodes as covered.
            coverageSet.mergeSampleCoverage(*candidateSet, coalTrackingSet, individualCoalCount);
            if (candidateId >= shapeNodeIdMax) {
                stats.reusedMutNodes++;
            }

            // Use this candidate (or the nodes below it) to cover the sample subset.
            stats.reusedNodes++;
            grg->connect(mutNodeId, candidateId);
            stats.reusedNodeCoverage += candidateSet->numSamplesCovered();
            if (candidateSet->numSamplesCovered() >= stats.reuseSizeHist.size()) {
                stats.reuseSizeBiggerThanHistMax++;
            } else {
                stats.reuseSizeHist[candidateSet->numSamplesCovered()]++;
            }
        }
        candidates.pop_back();
    }

    // Any leftovers, we just connect directly from the new mutation node to the
    // samples.
    NodeIDSet uncovered;
    for (const NodeID sampleNodeId : mutSamples) {
        if (!coverageSet.contains(sampleNodeId)) {
            uncovered.emplace(sampleNodeId);
            // The individual had already been seen and >=1 of the samples was previously uncovered,
            // then the new node we create is going to be the coalescence location for that individual.
            if (ploidy == 2) {
                NodeID individualId = (sampleNodeId / 2) * 2;
                if (coalTrackingSet.contains(individualId)) {
                    individualCoalCount++;
                } else {
                    coalTrackingSet.addElem(individualId);
                }
            }
        }
    }
    if (ploidy == 2) {
        grg->setNumIndividualCoals(mutNodeId, individualCoalCount);
    }

    if (!uncovered.empty()) {
        stats.numWithSingletons++;
    }

    stats.maxSingletons = std::max(uncovered.size(), stats.maxSingletons);

    for (auto sampleNodeId : uncovered) {
        grg->connect(mutNodeId, sampleNodeId);
        stats.singletonSampleEdges++;
    }
    // This node needs to be last, for the way we update things.
    addedNodes.push_back(mutNodeId);
    return addedNodes;
}

static NodeIDList greedyAddMutation(const MutableGRGPtr& grg,
                                    const MutationBatch& mutBatch,
                                    MutationMappingStats& stats,
                                    const NodeID shapeNodeIdMax,
                                    std::pair<BpPosition, NodeID>& currentMissing,
                                    size_t numThreads) {
    // The topological order of nodeIDs is maintained through-out this algorithm, because newly added
    // nodes are only ever _root nodes_ (at the time they are added).
    release_assert(grg->nodesAreOrdered());
    const NodeIDList& mutSamples = mutBatch.seedList();

    TopoCandidateCollectorVisitor collector(mutBatch);
    if (numThreads == 1) {
        grg->visitTopo(collector, grgl::TraversalDirection::DIRECTION_UP, mutSamples);
    } else {
        grg->visitTopoParallel(collector, grgl::TraversalDirection::DIRECTION_UP, mutSamples, numThreads);
    }

    NodeIDList totalAdded{};
    std::vector<CandidateList>& candidateList = collector.m_collectedNodes;
    for (int i = 0; i < mutBatch.numMutations(); i++) {
        NodeIDList added = processCandidateSet(grg,
                                               mutBatch.getMutation(i),
                                               mutBatch.sampleSet(i),
                                               collector,
                                               candidateList[i],
                                               stats,
                                               shapeNodeIdMax,
                                               currentMissing);
        totalAdded.insert(totalAdded.end(), added.begin(), added.end());
    }
    return totalAdded;
}

MutationMappingStats
mapMutations(const MutableGRGPtr& grg, MutationIterator& mutations, bool verbose, size_t numThreads, size_t mutationBatchSize) {
    auto operationStartTime = std::chrono::high_resolution_clock::now();
#define START_TIMING_OPERATION() operationStartTime = std::chrono::high_resolution_clock::now();
#define EMIT_TIMING_MESSAGE(msg)                                                                                       \
    do {                                                                                                               \
        std::cerr << msg                                                                                               \
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - \
                                                                           operationStartTime)                         \
                         .count()                                                                                      \
                  << " ms" << std::endl;                                                                               \
    } while (0)

    if (!grg->nodesAreOrdered()) {
        throw ApiMisuseFailure(
            "mapMutations can only be used on GRGs with ordered nodes; saving/loading a GRG will do this.");
    }

    MutationMappingStats stats;
    stats.reuseSizeHist.resize(STATS_HIST_SIZE, 0);
    stats.totalMutations = mutations.countMutations();

    // For the whole graph, count the number of samples under each node.
    if (verbose) {
        std::cout << "Mapping " << stats.totalMutations << " mutations\n";
    }
    const size_t onePercent = (stats.totalMutations / ONE_HUNDRED_PERCENT) + 1;
    size_t completed = 0;

    // The low-water mark for nodes. If a NodeID is greater than or equal to this, then it
    // is a newly added (mutation) node.
    const NodeID shapeNodeIdMax = grg->numNodes();
    MutationBatch mutBatch(grg->numSamples(), mutationBatchSize);

    // For each mutation, perform a topological bottom-up traversal from the sample
    // nodes of interest, and collect all nodes that reach a subset of those nodes.
    size_t _ignored = 0;
    std::pair<BpPosition, NodeID> currentMissing = {INVALID_POSITION, INVALID_NODE_ID};
    MutationAndSamples unmapped = {Mutation(0.0, ""), NodeIDList()};
    while (mutations.next(unmapped, _ignored)) {
        const NodeIDList& mutSamples = unmapped.samples;
        if (!mutSamples.empty()) {
            stats.samplesProcessed += mutSamples.size();
            // quickly skip mutations with only one sample, as we don't explicitly handle them in our batched version
            if (mutSamples.size() == 1) {
                stats.mutationsWithOneSample++;
                stats.mutationsWithNoCandidates++;

                const NodeID sampleNodeId = mutSamples[0];
                const NodeID missingnessNode =
                    (currentMissing.first == unmapped.mutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;

                const int ploidy = grg->getPloidy();

                const NodeID mutNodeId = grg->makeNode(1, true);
                if (!unmapped.mutation.isMissing()) {
                    grg->addMutation(unmapped.mutation, mutNodeId, missingnessNode);
                } else {
                    currentMissing = {unmapped.mutation.getPosition(), mutNodeId};
                }

                // Coalescence is always 0 for singletons
                if (ploidy == 2) {
                    grg->setNumIndividualCoals(mutNodeId, 0);
                }

                stats.numWithSingletons++;
                stats.maxSingletons = std::max<size_t>(1, stats.maxSingletons);
                stats.singletonSampleEdges++;
                grg->connect(mutNodeId, sampleNodeId);

                NodeIDList singletonAdded{mutNodeId};
            } else {
                mutBatch.addMutation(unmapped.mutation, mutSamples);

                if (mutBatch.numMutations() == mutationBatchSize) {
                    NodeIDList addedNodes =
                        greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, numThreads);

                    // Update sample counts for newly added nodes.
                    mutBatch.clear();
                }
            }
        } else {
            if (mutBatch.numMutations() > 0) {
                NodeIDList addedNodes =
                    greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, numThreads);
                mutBatch.clear();
            }
            stats.emptyMutations++;
            const NodeID missingnessNode =
                (currentMissing.first == unmapped.mutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;
            grg->addMutation(unmapped.mutation, INVALID_NODE_ID, missingnessNode);
            api_exc_check(!unmapped.mutation.isMissing(), "Missing data rows cannot have no samples");
        }
        completed++;
        const size_t percentCompleted = (completed / onePercent);
        if (verbose) {
            if ((completed % onePercent == 0)) {
                std::cout << percentCompleted << "% done" << std::endl;
            }
            if ((completed % (EMIT_STATS_AT_PERCENT * onePercent) == 0)) {
                std::cout << "Last mutation sampleset size: " << mutSamples.size() << std::endl;
                std::cout << "GRG nodes: " << grg->numNodes() << std::endl;
                std::cout << "GRG edges: " << grg->numEdges() << std::endl;
                stats.print(std::cout);
            }
            if ((completed % (COMPACT_EDGES_AT_PERCENT * onePercent) == 0)) {
                START_TIMING_OPERATION();
                grg->compact();
                EMIT_TIMING_MESSAGE("Compacting GRG edges took ");
            }
        }
    }
    if (mutBatch.numMutations() > 0) {
        NodeIDList addedNodes = greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, numThreads);
    }
    return stats;
}
}; // namespace grgl
