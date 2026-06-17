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
#include <functional>
#include <iostream>
#include <memory>
#include <sys/types.h>
#include <unordered_map>
#include <utility>
#include <vector>

#if MAP_MUTS_THREADED
#include <mutex>
#define DEFINE_MUTEX(name)  std::mutex name;
#define CREATE_GUARD(mutex) std::lock_guard<std::mutex> lock(mutex);
#else
#define DEFINE_MUTEX(name)
#define CREATE_GUARD(mutex)
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif
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

struct CandidatePlan {
    bool exactMatch{false};
    NodeID exactMatchCandidate{INVALID_NODE_ID};
    NodeIDList reusedCandidates;
    NodeIDList uncoveredSamples;
    size_t individualCoalCount{0};
    MutationMappingStats stats;
};

static bool cmpNodeSamples(const NodeSamples& ns1, const NodeSamples& ns2) {
    const size_t samplesCount1 = std::get<1>(ns1);
    const size_t samplesCount2 = std::get<1>(ns2);
    return samplesCount1 < samplesCount2;
}

static void accumulateStats(MutationMappingStats& total, const MutationMappingStats& delta) {
    total.mutationsWithNoCandidates += delta.mutationsWithNoCandidates;
    total.reusedNodes += delta.reusedNodes;
    total.reusedNodeCoverage += delta.reusedNodeCoverage;
    total.reusedExactly += delta.reusedExactly;
    total.singletonSampleEdges += delta.singletonSampleEdges;
    total.numWithSingletons += delta.numWithSingletons;
    total.maxSingletons = std::max(total.maxSingletons, delta.maxSingletons);
    total.reusedMutNodes += delta.reusedMutNodes;
    total.reuseSizeBiggerThanHistMax += delta.reuseSizeBiggerThanHistMax;
    if (total.reuseSizeHist.size() < delta.reuseSizeHist.size()) {
        total.reuseSizeHist.resize(delta.reuseSizeHist.size(), 0);
    }
    for (size_t i = 0; i < delta.reuseSizeHist.size(); ++i) {
        total.reuseSizeHist[i] += delta.reuseSizeHist[i];
    }
}

// What percentage of total samples must be covered before we switch to dense bitvecs
constexpr double USE_DENSE_COVERAGE_PROPORTION(0.001f);

// SampleCoverageSet is the interface for tracking which sample nodes lie beneath a GRG node.
class SampleCoverageSet {
public:
    /// the number of total samples
    virtual size_t totalSamples() const = 0;

    /// number of elements in the set, INCLUDING DUPLICATES from below nodes
    virtual size_t numSamplesCovered() const = 0;

    virtual bool empty() const = 0;
    virtual void addElem(NodeID sample) = 0;

    void mergeSampleCoverage(const SampleCoverageSet& other) { mergeSampleCoverageImpl(other); }

    void mergeSampleCoverage(const SampleCoverageSet& other, BitVecCoverageSet& seenSet, size_t& coalescences) {
        mergeSampleCoverageImpl(other, seenSet, coalescences);
    }

    void mergeSampleCoverage(const SampleCoverageSet& other,
                             std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                             NodeID childNode,
                             size_t& coalescences) {
        mergeSampleCoverageImpl(other, individualToChild, childNode, coalescences);
    }

    virtual BatchMembership getBatchMembership() const = 0;

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
};

// membership bits stored as uint32_t words; sample list stored as NodeID (32-bit) for small coverage.
class SparseCoverageSet : public SampleCoverageSet {
public:
    explicit SparseCoverageSet(size_t batchSize, size_t totalSampleCount, size_t reserve_amount = 0)
        : m_totalSampleCount(totalSampleCount),
          m_batchSize(batchSize) {
        m_membership32.assign((m_batchSize + 31) / 32, 0u);
        if (reserve_amount > 0) {
            m_samples.reserve(reserve_amount);
        }
    }

    size_t totalSamples() const override { return m_totalSampleCount; }
    size_t numSamplesCovered() const override { return m_samples.size(); }

    void addElem(NodeID sample) override { m_samples.push_back(sample); }

    BatchMembership getBatchMembership() const override { return BatchMembership{m_membership32.data(), m_batchSize}; }

    bool empty() const override { return m_samples.empty(); }
    size_t batchSize() const { return m_batchSize; }

    const NodeID* samplesBegin() const { return m_samples.data(); }
    const NodeID* samplesEnd() const { return m_samples.data() + m_samples.size(); }

protected:
    void mergeSampleCoverageImpl(const SampleCoverageSet& other) override {
        const auto* otherSparse = dynamic_cast<const SparseCoverageSet*>(&other);
        if (otherSparse == nullptr || otherSparse->empty()) {
            return;
        }
        getBatchMembership().andAssign(otherSparse->getBatchMembership());
        appendSamples(*otherSparse);
    }

    void mergeSampleCoverageImpl(const SampleCoverageSet& other,
                                 std::unordered_map<NodeIDSizeT, NodeIDSizeT>& individualToChild,
                                 NodeID childNode,
                                 size_t& coalescences) override {
        const auto* otherSparse = dynamic_cast<const SparseCoverageSet*>(&other);
        if (otherSparse == nullptr || otherSparse->empty()) {
            return;
        }
        getBatchMembership().andAssign(otherSparse->getBatchMembership());
        for (auto it = otherSparse->samplesBegin(); it != otherSparse->samplesEnd(); ++it) {
            const NodeID sampleId = *it;
            m_samples.push_back(sampleId);
            const auto insertPair = individualToChild.emplace(sampleId / 2, childNode);
            if (!insertPair.second && insertPair.first->second != childNode) {
                coalescences++;
            }
        }
    }

    friend class BitVecCoverageSet;

private:
    void appendSamples(const SparseCoverageSet& other) {
        m_samples.insert(m_samples.end(), other.m_samples.begin(), other.m_samples.end());
    }

    mutable std::vector<uint32_t> m_membership32;
    mutable std::vector<NodeID> m_samples;
    size_t m_totalSampleCount;
    size_t m_batchSize;
};

// dense bitset representation of sample coverage and membership for large coverage sets
class BitVecCoverageSet : public SampleCoverageSet {
public:
    friend class CoalescenceVisitor;

    explicit BitVecCoverageSet(size_t totalSamples, size_t batchSize = 0)
        : m_totalSampleCount(totalSamples),
          m_batchSize(batchSize) {
        m_membership32.assign((m_batchSize + 31) / 32, 0u);
        const size_t sampleWords = (totalSamples + m_blockSize - 1) / m_blockSize;
        m_elems.assign(sampleWords, 0ull);
    }

    explicit BitVecCoverageSet(const SparseCoverageSet& sparseCoverageSet)
        : BitVecCoverageSet(sparseCoverageSet.totalSamples(), sparseCoverageSet.batchSize()) {
        mergeSampleCoverage(sparseCoverageSet);
    }

    size_t totalSamples() const override { return m_totalSampleCount; }
    bool empty() const override { return m_numSamples == 0; }

    void addElem(NodeID sample) override {
        size_t vecIdx = static_cast<size_t>(sample) / m_blockSize;
        size_t offset = static_cast<size_t>(sample) % m_blockSize;
        // if this bit was zero before, increment by one
        m_numSamples += ~(m_elems[vecIdx] >> offset) & 1ull;
        m_elems[vecIdx] |= (1ull << offset);
    }

    BatchMembership getBatchMembership() const override { return BatchMembership{m_membership32.data(), m_batchSize}; }

    bool contains(const NodeID sample) const {
        size_t vecIdx = static_cast<size_t>(sample) / m_blockSize;
        size_t offset = static_cast<size_t>(sample) % m_blockSize;
        return ((m_elems[vecIdx] >> offset) & 1ull) != 0;
    }

    // original api does not dedup sets when computing sizes.
    size_t numSamplesNonOverlapping() const { return m_numSamples; }
    size_t numSamplesCovered() const override { return m_numSamples; }

    bool overlapsWith(const SampleCoverageSet& other) const {
        release_assert(totalSamples() == other.totalSamples());
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            if (numSamplesCovered() + other_v->numSamplesCovered() > totalSamples()) {
                return true;
            }
            const uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();
            for (size_t i = 0; i < len; i++) {
                if ((data[i] & otherData[i]) != 0) {
                    return true;
                }
            }
        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (auto it = other_v->samplesBegin(); it != other_v->samplesEnd(); ++it) {
                const NodeID sample = *it;
                size_t vecIdx = static_cast<size_t>(sample) / m_blockSize;
                size_t offset = static_cast<size_t>(sample) % m_blockSize;
                if ((m_elems[vecIdx] & (1ull << offset)) != 0) {
                    return true;
                }
            }
        }
        return false;
    }

protected:
    void mergeSampleCoverageImpl(const SampleCoverageSet& other) override {
        release_assert(totalSamples() == other.totalSamples());
        getBatchMembership().andAssign(other.getBatchMembership());

        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();
            for (size_t i = 0; i < len; ++i) {
                m_numSamples += __builtin_popcountll(otherData[i] & ~data[i]);
                data[i] |= otherData[i];
            }
        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (auto it = other_v->samplesBegin(); it != other_v->samplesEnd(); ++it) {
                addElem(*it);
            }
        }
    }

    void
    mergeSampleCoverageImpl(const SampleCoverageSet& other, BitVecCoverageSet& seenSet, size_t& coalescences) override {
        release_assert(totalSamples() == other.totalSamples());
        getBatchMembership().andAssign(other.getBatchMembership());

        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            uint64_t* __restrict seen_d = seenSet.m_elems.data();
            const size_t len = m_elems.size();
            constexpr uint64_t EVEN = 0x5555555555555555ULL; // all even bits
            for (size_t i = 0; i < len; ++i) {
                const uint64_t word = otherData[i];
                m_numSamples += __builtin_popcountll(word & ~data[i]);
                data[i] |= word;
                const uint64_t coalesced =
                    (word | (word >> 1ULL)) & EVEN; // pack two consecutive haplotype bits into one individual slot
                const uint64_t new_coals = coalesced & seen_d[i];
                coalescences += __builtin_popcountll(new_coals);
                seen_d[i] |= coalesced;
            }
        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            auto& seenWords = seenSet.m_elems;
            for (const auto* it = other_v->samplesBegin(); it != other_v->samplesEnd(); ++it) {
                const NodeID sample = *it;
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
    mutable std::vector<uint32_t> m_membership32;
    std::vector<uint64_t> m_elems;
    size_t m_totalSampleCount{};
    size_t m_numSamples{};
    size_t m_batchSize{};
};

// MutationBatch groups mutations and their sample sets for batch processing and candidate collection.
class MutationBatch {
    friend class TopoCandidateCollectorVisitor;

public:
    explicit MutationBatch(NodeIDSizeT totalSampleCount, size_t batchBitCount)
        : m_batchBitCount(batchBitCount),
          m_sampleMembership(totalSampleCount) {
        const size_t numWords32 = (m_batchBitCount + 31) / 32;
        m_membershipStorage.assign(totalSampleCount * numWords32, 0u);

        for (size_t sampleIdx = 0; sampleIdx < totalSampleCount; ++sampleIdx) {
            uint32_t* base32 = m_membershipStorage.data() + (sampleIdx * numWords32);
            m_sampleMembership[sampleIdx] = BatchMembership(base32, m_batchBitCount);
            m_sampleMembership[sampleIdx].clear();
        }
    }

    void addMutation(const Mutation& mutation, NodeIDList mutSamples) {
        size_t mutIdx = m_mutations.size();
        release_assert(mutIdx < m_batchBitCount);
        m_mutations.emplace_back(mutation, m_seedList.size(), mutSamples.size());
        m_seedList.insert(m_seedList.end(), mutSamples.begin(), mutSamples.end());
        m_sampleSets.push_back(mutSamples);
        for (NodeID sample : mutSamples) {
            m_sampleMembership[sample].setBit(mutIdx);
        }
    }

    const NodeIDList& seedList() const { return m_seedList; }
    const Mutation& getMutation(size_t idx) const { return m_mutations.at(idx).m_mutation; }
    const NodeIDList& sampleSet(size_t idx) const { return m_sampleSets.at(idx); }
    size_t numMutations() const { return m_mutations.size(); }

    void clear() {
        m_sampleSets.clear();
        m_mutations.clear();
        m_seedList.clear();
        std::fill(m_membershipStorage.begin(), m_membershipStorage.end(), 0U);
    }

    MutationBatch(const MutationBatch&) = delete;
    MutationBatch& operator=(const MutationBatch&) = delete;

    MutationBatch(MutationBatch&&) = default;
    MutationBatch& operator=(MutationBatch&&) = default;

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
    std::vector<uint32_t> m_membershipStorage; // membership storage is 32-bit
};

// TopoCandidateCollectorVisitor walks the GRG to collect candidate nodes and sample coverage for each mutation bit.
class TopoCandidateCollectorVisitor : public grgl::GRGVisitor {
public:
    explicit TopoCandidateCollectorVisitor(const MutationBatch& mutBatch)
        : m_batchBitCount(mutBatch.m_batchBitCount),
          m_sampleBatchMembership(mutBatch.m_sampleMembership),
          m_collectedNodes(mutBatch.m_batchBitCount) {}
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

    const SampleCoverageSet& getSamplesForCandidate(NodeID candidateId) const {
        const std::unique_ptr<SampleCoverageSet>& value = m_sampleCoverage[candidateId];
        release_assert(value);
        return *value;
    }

    size_t batchBitCount() const { return m_batchBitCount; }

    std::vector<std::unique_ptr<SampleCoverageSet>> m_sampleCoverage;
    std::vector<CandidateList> m_collectedNodes;

private:
    size_t m_batchBitCount;
    const std::vector<BatchMembership>& m_sampleBatchMembership;
    // shared 'thread-safe' logic
    bool safeVisit(const grgl::GRGPtr& grg, grgl::NodeID nodeId);
    DEFINE_MUTEX(m_collectedNodesMutex)
#if CLEANUP_SAMPLE_SETS_MAPPING
    std::vector<NodeIDSizeT> m_refCounts;
#endif
};

bool TopoCandidateCollectorVisitor::safeVisit(const grgl::GRGPtr& grg, const grgl::NodeID nodeId) {
    release_assert(grg->hasUpEdges());
    // Note: Any modification to current node or child is safe, as we never process parent and child at the same time
    const bool isRoot = grg->numUpEdges(nodeId) == 0;
    const bool isSample = grg->isSample(nodeId);
    const size_t ploidy = grg->getPloidy();
    const auto numCoals = grg->getNumIndividualCoals(nodeId);
    const bool computeCoals = !isSample && (ploidy == PLOIDY_COAL_PROP) && (COAL_COUNT_NOT_SET == numCoals);

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
    bool dense = coveragePercent > USE_DENSE_COVERAGE_PROPORTION;
    if (dense) {
        samplesBeneathSet =
            std::unique_ptr<SampleCoverageSet>(new BitVecCoverageSet(grg->numSamples(), m_batchBitCount));
        coalTrackerSet = std::unique_ptr<BitVecCoverageSet>(new BitVecCoverageSet(grg->numSamples(), m_batchBitCount));
    } else {
        size_t sparse_reserve_amount = sampleCount;
        samplesBeneathSet = std::unique_ptr<SampleCoverageSet>(
            new SparseCoverageSet(m_batchBitCount, grg->numSamples(), sparse_reserve_amount));
    }

    samplesBeneathSet->getBatchMembership().setOnes(m_batchBitCount);

    if (isSample) {
        samplesBeneathSet->addElem(nodeId);
        auto membership = samplesBeneathSet->getBatchMembership();
        membership.copyFrom(m_sampleBatchMembership[nodeId]);
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
            samplesBeneathSet->getBatchMembership().clear();
        }
    }

    // We can only record coalescence counts if there are no samples missing.
    if (samplesBeneathSet->getBatchMembership().any() && computeCoals) {
        grg->setNumIndividualCoals(nodeId, individualCoalCount);
    }

    // If we've reached the root of the graph or have missing samples beneath us, we need to stop the search.
    const BatchMembership& membership = samplesBeneathSet->getBatchMembership();
    const bool keepGoing = !isRoot && membership.any();

    if (isRoot) {
        // Root can only be a candidate when it fully covers the batch for at least one bit.
        if (membership.any()) {
            CREATE_GUARD(m_collectedNodesMutex)
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
                CREATE_GUARD(m_collectedNodesMutex)
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
            CREATE_GUARD(m_collectedNodesMutex)
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
static CandidatePlan greedyAddMutationImmutable(const MutableGRGPtr& grg,
                                                const NodeIDList& mutSamples,
                                                const TopoCandidateCollectorVisitor& collector,
                                                CandidateList candidates,
                                                const NodeID shapeNodeIdMax) {
    CandidatePlan plan;
    plan.stats.reuseSizeHist.resize(STATS_HIST_SIZE, 0);

    std::sort(candidates.begin(), candidates.end());
    auto endOfUnique = std::unique(candidates.begin(), candidates.end());
    candidates.erase(endOfUnique, candidates.end());
    std::sort(candidates.begin(), candidates.end(), cmpNodeSamples);

    const int ploidy = grg->getPloidy();
    if (candidates.empty()) {
        plan.stats.mutationsWithNoCandidates++;
    } else {
        const auto& candidate = candidates.back();
        const size_t candidateSetSize = std::get<1>(candidate);
        if (candidateSetSize == mutSamples.size()) {
            plan.exactMatch = true;
            plan.exactMatchCandidate = std::get<0>(candidate);
            plan.stats.reusedExactly++;
            if (plan.exactMatchCandidate >= shapeNodeIdMax) {
                plan.stats.reusedMutNodes++;
            }
            return plan;
        }
    }

    const size_t batchBitCount = collector.batchBitCount();
    BitVecCoverageSet coalTrackingSet{grg->numSamples(), batchBitCount};
    BitVecCoverageSet coverageSet{grg->numSamples(), batchBitCount};
    const size_t numMutSamples = mutSamples.size();

    while (!candidates.empty() && coverageSet.numSamplesNonOverlapping() < numMutSamples) {
        const auto& candidate = candidates.back();
        const auto candidateId = std::get<0>(candidate);
        const SampleCoverageSet& candidateSet = collector.getSamplesForCandidate(candidateId);
        release_assert(!candidateSet.empty());

        if (!coverageSet.overlapsWith(candidateSet)) {
            if (ploidy == PLOIDY_COAL_PROP) {
                coverageSet.mergeSampleCoverage(candidateSet, coalTrackingSet, plan.individualCoalCount);
            } else {
                coverageSet.mergeSampleCoverage(candidateSet);
            }

            if (candidateId >= shapeNodeIdMax) {
                plan.stats.reusedMutNodes++;
            }

            plan.stats.reusedNodes++;
            plan.reusedCandidates.push_back(candidateId);
            plan.stats.reusedNodeCoverage += candidateSet.numSamplesCovered();
            if (candidateSet.numSamplesCovered() >= plan.stats.reuseSizeHist.size()) {
                plan.stats.reuseSizeBiggerThanHistMax++;
            } else {
                plan.stats.reuseSizeHist[candidateSet.numSamplesCovered()]++;
            }
        }
        candidates.pop_back();
    }

    for (const NodeID sampleNodeId : mutSamples) {
        if (!coverageSet.contains(sampleNodeId)) {
            plan.uncoveredSamples.emplace_back(sampleNodeId);
            if (ploidy == PLOIDY_COAL_PROP) {
                const NodeID individualId = (sampleNodeId / 2) * 2;
                if (coalTrackingSet.contains(individualId)) {
                    plan.individualCoalCount++;
                } else {
                    coalTrackingSet.addElem(individualId);
                }
            }
        }
    }

    if (!plan.uncoveredSamples.empty()) {
        plan.stats.numWithSingletons++;
    }
    plan.stats.maxSingletons = std::max(plan.uncoveredSamples.size(), plan.stats.maxSingletons);
    plan.stats.singletonSampleEdges += plan.uncoveredSamples.size();
    return plan;
}

static NodeIDList applyBatchModifications(const MutableGRGPtr& grg,
                                          const MutationBatch& mutBatch,
                                          const std::vector<CandidatePlan>& batchResults,
                                          std::pair<BpPosition, NodeID>& currentMissing) {
    NodeIDList added;
    const size_t batchSize = mutBatch.numMutations();
    for (size_t i = 0; i < batchSize; ++i) {
        const CandidatePlan& plan = batchResults[i];
        const Mutation& mut = mutBatch.getMutation(i);
        const NodeID missingnessNode =
            (currentMissing.first == mut.getPosition()) ? currentMissing.second : INVALID_NODE_ID;

        if (plan.exactMatch) {
            if (!mut.isMissing()) {
                grg->addMutation(mut, plan.exactMatchCandidate, missingnessNode);
            } else {
                currentMissing = {mut.getPosition(), plan.exactMatchCandidate};
            }
            continue;
        }

        const NodeID mutNodeId = grg->makeNode(1, true);
        if (!mut.isMissing()) {
            grg->addMutation(mut, mutNodeId, missingnessNode);
        } else {
            currentMissing = {mut.getPosition(), mutNodeId};
        }

        if (grg->getPloidy() == PLOIDY_COAL_PROP) {
            grg->setNumIndividualCoals(mutNodeId, plan.individualCoalCount);
        }

        for (const auto candidateId : plan.reusedCandidates) {
            grg->connect(mutNodeId, candidateId);
        }
        for (const auto sampleNodeId : plan.uncoveredSamples) {
            grg->connect(mutNodeId, sampleNodeId);
        }
        added.push_back(mutNodeId);
    }
    return added;
}

static NodeIDList processBatchPar(const MutableGRGPtr& grg,
                                  const MutationBatch& mutBatch,
                                  const TopoCandidateCollectorVisitor& collector,
                                  const NodeID shapeNodeIdMax,
                                  std::pair<BpPosition, NodeID>& currentMissing,
                                  MutationMappingStats& stats,
                                  const size_t threadCount) {
    const size_t batchSize = mutBatch.numMutations();
    std::vector<CandidatePlan> batchTasks(batchSize);
    std::vector<MutationMappingStats> localStats(batchSize);

    if (batchSize == 0) {
        return {};
    }

#if defined(_OPENMP)
    const size_t numThreads = threadCount == 0 ? static_cast<size_t>(omp_get_max_threads()) : threadCount;
    if (numThreads <= 1) {
        for (size_t i = 0; i < batchSize; ++i) {
            batchTasks[i] = greedyAddMutationImmutable(
                grg, mutBatch.sampleSet(i), collector, collector.m_collectedNodes[i], shapeNodeIdMax);
            localStats[i] = batchTasks[i].stats;
        }
        for (const auto& threadStat : localStats) {
            accumulateStats(stats, threadStat);
        }
        return applyBatchModifications(grg, mutBatch, batchTasks, currentMissing);
    }
    omp_set_num_threads(static_cast<int>(numThreads));
#pragma omp parallel
    {
#pragma omp single
        {
            for (size_t i = 0; i < batchSize; ++i) {
#pragma omp task firstprivate(i)
                {
                    batchTasks[i] = greedyAddMutationImmutable(
                        grg, mutBatch.sampleSet(i), collector, collector.m_collectedNodes[i], shapeNodeIdMax);
                    localStats[i] = batchTasks[i].stats;
                }
            }
#pragma omp taskwait
        }
    }
#else
    for (size_t i = 0; i < batchSize; ++i) {
        batchTasks[i] = greedyAddMutationImmutable(
            grg, mutBatch.sampleSet(i), collector, collector.m_collectedNodes[i], shapeNodeIdMax);
        localStats[i] = batchTasks[i].stats;
    }
#endif

    for (const auto& threadStat : localStats) {
        accumulateStats(stats, threadStat);
    }

    return applyBatchModifications(grg, mutBatch, batchTasks, currentMissing);
}

static NodeIDList greedyAddMutation(const MutableGRGPtr& grg,
                                    const MutationBatch& mutBatch,
                                    MutationMappingStats& stats,
                                    const NodeID shapeNodeIdMax,
                                    std::pair<BpPosition, NodeID>& currentMissing,
                                    const size_t threadCount) {
    // The topological order of nodeIDs is maintained throughout this algorithm.
    release_assert(grg->nodesAreOrdered());
    const NodeIDList& mutSamples = mutBatch.seedList();

    TopoCandidateCollectorVisitor collector(mutBatch);
    grg->visitTopo(collector, grgl::TraversalDirection::DIRECTION_UP, mutSamples);

    return processBatchPar(grg, mutBatch, collector, shapeNodeIdMax, currentMissing, stats, threadCount);
}

static MutationMappingStats
mapMutationsImpl(const MutableGRGPtr& grg,
                 size_t totalMutations,
                 std::function<bool(MutationAndSamples&, size_t&)> next,
                 bool verbose,
                 size_t mutationBatchSize,
                 size_t threadCount) {
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
    stats.totalMutations = totalMutations;

    if (verbose) {
        std::cout << "Mapping " << stats.totalMutations << " mutations\n";
    }
    const size_t onePercent = (stats.totalMutations / ONE_HUNDRED_PERCENT) + 1;
    size_t completed = 0;

    const NodeID shapeNodeIdMax = grg->numNodes();
    MutationBatch mutBatch(grg->numSamples(), mutationBatchSize);

    size_t _ignored = 0;
    std::pair<BpPosition, NodeID> currentMissing = {INVALID_POSITION, INVALID_NODE_ID};
    MutationAndSamples unmapped = {Mutation(0.0, ""), NodeIDList()};

    const int ploidy = grg->getPloidy();
    while (next(unmapped, _ignored)) {
        const NodeIDList& mutSamples = unmapped.samples;
        if (!mutSamples.empty()) {
            stats.samplesProcessed += mutSamples.size();
            // quickly skip mutations with only one sample
            if (mutSamples.size() == 1) {
                stats.mutationsWithOneSample++;
                stats.mutationsWithNoCandidates++;

                const NodeID sampleNodeId = mutSamples[0];
                const NodeID missingnessNode =
                    (currentMissing.first == unmapped.mutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;

                const NodeID mutNodeId = sampleNodeId;
                if (!unmapped.mutation.isMissing()) {
                    grg->addMutation(unmapped.mutation, mutNodeId, missingnessNode);
                } else {
                    currentMissing = {unmapped.mutation.getPosition(), mutNodeId};
                }

                // Coalescence is always 0 for singletons
                if (ploidy == PLOIDY_COAL_PROP) {
                    grg->setNumIndividualCoals(mutNodeId, 0);
                }
            } else {
                mutBatch.addMutation(unmapped.mutation, mutSamples);

                if (mutBatch.numMutations() == mutationBatchSize) {
                    greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, threadCount);
                    mutBatch.clear();
                }
            }
        } else {
            if (mutBatch.numMutations() > 0) {
                greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, threadCount);
                mutBatch.clear();
            }
            stats.emptyMutations++;
            const NodeID missingnessNode =
                (currentMissing.first == unmapped.mutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;
            grg->addMutation(unmapped.mutation, INVALID_NODE_ID, missingnessNode);
            api_exc_check(!unmapped.mutation.isMissing(), "Missing data rows cannot have no samples");
        }

        completed++;
        if (verbose) {
            if ((completed % onePercent) == 0) {
                const size_t percentCompleted = (completed / onePercent);
                std::cout << percentCompleted << "% done" << std::endl;
            }
            if ((completed % (EMIT_STATS_AT_PERCENT * onePercent)) == 0) {
                std::cout << "Last mutation sampleset size: " << mutSamples.size() << std::endl;
                std::cout << "GRG nodes: " << grg->numNodes() << std::endl;
                std::cout << "GRG edges: " << grg->numEdges() << std::endl;
                stats.print(std::cout);
            }
            if ((completed % (COMPACT_EDGES_AT_PERCENT * onePercent)) == 0) {
                START_TIMING_OPERATION();
                grg->compact();
                EMIT_TIMING_MESSAGE("Compacting GRG edges took ");
            }
        }
    }

    if (mutBatch.numMutations() > 0) {
        greedyAddMutation(grg, mutBatch, stats, shapeNodeIdMax, currentMissing, threadCount);
    }
    return stats;
}

MutationMappingStats
mapMutations(const MutableGRGPtr& grg,
             MutationIterator& mutations,
             bool verbose,
             size_t mutationBatchSize,
             size_t threadCount) {
    auto next = [&mutations](MutationAndSamples& unmapped, size_t& totalSamples) {
        return mutations.next(unmapped, totalSamples);
    };

    return mapMutationsImpl(grg, mutations.countMutations(), std::move(next), verbose, mutationBatchSize, threadCount);
}

MutationMappingStats mapMutations(const MutableGRGPtr& grg,
                                  const std::vector<Mutation>& mutations,
                                  const std::vector<NodeIDList>& samples,
                                  bool verbose,
                                  size_t mutationBatchSize,
                                  size_t threadCount) {
    api_exc_check(mutations.size() == samples.size(),
                  "mutations and samples must be the same length in mapMutations");

    size_t idx = 0;
    auto next = [&mutations, &samples, &idx](MutationAndSamples& unmapped, size_t& totalSamples) {
        if (idx >= mutations.size()) {
            return false;
        }
        unmapped.mutation = mutations[idx];
        unmapped.samples = samples[idx];
        totalSamples = unmapped.samples.size();
        idx++;
        return true;
    };

    return mapMutationsImpl(grg, mutations.size(), std::move(next), verbose, mutationBatchSize, threadCount);
}

} // namespace grgl
