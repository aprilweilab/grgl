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
#include <atomic>
#include <iostream>
#include <memory>
#include <omp.h>
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

static bool cmpNodeSamples(const NodeSamples& ns1, const NodeSamples& ns2) {
    const size_t samplesCount1 = std::get<1>(ns1);
    const size_t& samplesCount2 = std::get<1>(ns2);
    return samplesCount1 < samplesCount2;
}

// What percentage of total samples must be covered before we switch to dense bitvecs
#define USE_DENSE_COVERAGE_PERCENT (0.05f)

class BitVecCoverageSet;
class SparseCoverageSet;
struct CoalescenceTracker {
    friend class BitVecCoverageSet;
    friend class SparseCoverageSet;

    CoalescenceTracker() = default;
    CoalescenceTracker(BitVecCoverageSet* seen, size_t* individualCoalCount)
        : m_seenSet(seen),
          m_coalsAdded(individualCoalCount) {}

    CoalescenceTracker(std::unordered_map<NodeIDSizeT, NodeIDSizeT>* individualToChild,
                       NodeID childNode,
                       size_t* individualCoalCount)
        : m_individualToChild(individualToChild),
          m_childNode(childNode),
          m_coalsAdded(individualCoalCount) {}

private:
    BitVecCoverageSet* m_seenSet = nullptr;
    std::unordered_map<NodeIDSizeT, NodeIDSizeT>* m_individualToChild = nullptr;
    size_t* m_coalsAdded = nullptr;
    NodeID m_childNode{};
};

class SampleCoverageSet {
public:
    /// the number of total samples
    virtual size_t totalSamples() const = 0;

    /// number of elements in the set, INCLUDING DUPLICATES from below nodes
    virtual size_t numSamplesCovered() const = 0;

    virtual bool empty() const = 0;
    virtual void addElem(NodeID sample) = 0;
    virtual void mergeSampleCoverage(const SampleCoverageSet& other, CoalescenceTracker tracker = {}) = 0;

    SampleCoverageSet() = default;
    SampleCoverageSet(SampleCoverageSet&&) = delete;
    SampleCoverageSet(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(SampleCoverageSet&&) = default;
    virtual ~SampleCoverageSet() = default;
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
    void mergeSampleCoverage(const SampleCoverageSet& other, CoalescenceTracker tracker = {}) override {
        // we can't merge a bitvec coverage set into a non-bitvec coverage set
        const auto* otherSparse = dynamic_cast<const SparseCoverageSet*>(&other);
        if (otherSparse == nullptr || otherSparse->m_sampleIdxs.empty()) {
            return;
        }

        // quick append when coalescence isn't tracked
        if (tracker.m_individualToChild == nullptr) {
            m_sampleIdxs.insert(m_sampleIdxs.end(), otherSparse->m_sampleIdxs.begin(), otherSparse->m_sampleIdxs.end());
            return;
        }

        auto& individualToChild = tracker.m_individualToChild;
        for (const NodeID sampleId : otherSparse->m_sampleIdxs) {
            m_sampleIdxs.emplace_back(sampleId);
            const auto insertPair = individualToChild->emplace(sampleId / 2, tracker.m_childNode);
            if (!insertPair.second && insertPair.first->second != tracker.m_childNode) {
                (*tracker.m_coalsAdded)++;
            }
        }
    }

    friend class BitVecCoverageSet;

private:
    NodeIDList m_sampleIdxs;
    size_t m_totalSampleCount;
};

/// A BitVec-backed set for coverage
class BitVecCoverageSet : public SampleCoverageSet {
public:
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

    void mergeSampleCoverage(const SampleCoverageSet& other, CoalescenceTracker tracker = {}) override {
        m_numSampleCoverage += other.numSamplesCovered();
        release_assert(totalSamples() == other.totalSamples());
        // todo: add asserts to make sure CoalescenceTracker bitvecs (if present) are sized correctly.
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();
            if (tracker.m_seenSet == nullptr) {
                for (int i = 0; i < len; i++) {
                    m_numSamplesNonOverlapping += __builtin_popcountll(otherData[i] & ~data[i]);
                    data[i] |= otherData[i];
                }
            } else {
                auto* seenBV = tracker.m_seenSet;
                uint64_t* __restrict seen_d = seenBV->m_elems.data();

                for (size_t i = 0; i < len; ++i) {
                    uint64_t word = otherData[i];
                    m_numSamplesNonOverlapping += __builtin_popcountll(word & ~data[i]);
                    data[i] |= word;
                    constexpr uint64_t EVEN = 0x5555555555555555ULL; // all even bits
                    const uint64_t coalesced =
                        (word | (word >> 1ULL)) & EVEN; // if either haplotype has it (i.e. shifted by 1)
                    // if it's already present from another child, we coalesce at this node
                    uint64_t new_coals = coalesced & seen_d[i];
                    *tracker.m_coalsAdded += __builtin_popcountll(new_coals);
                    seen_d[i] |= coalesced;
                }
            }

        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (const uint32_t sample : other_v->m_sampleIdxs) {
                addElem(sample);
            }

            if (tracker.m_seenSet != nullptr) {
                auto* seenBV = tracker.m_seenSet;
                auto& seenWords = seenBV->m_elems;

                for (const uint32_t sample : other_v->m_sampleIdxs) {
                    const size_t base = (static_cast<size_t>(sample) / 2) * 2; // even index per individual
                    const size_t vecIdx = base / m_blockSize;
                    const size_t offset = base % m_blockSize;
                    const uint64_t mask = (uint64_t{1} << offset);

                    if ((seenWords[vecIdx] & mask) != 0) {
                        (*tracker.m_coalsAdded)++;
                    }
                    seenWords[vecIdx] |= mask;
                }
            }
        }
    }

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

private:
    constexpr static size_t m_blockSize = 64;
    std::vector<uint64_t> m_elems;
    size_t m_totalSampleCount;
    size_t m_numSampleCoverage{};
    size_t m_numSamplesNonOverlapping{};
};

class TopoCandidateCollectorVisitor : public grgl::ParallelGRGVisitor {
public:
    explicit TopoCandidateCollectorVisitor(const std::vector<NodeIDSizeT>& sampleCounts)
        : m_sampleCounts(sampleCounts) {}

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
        for (NodeID node : nodes) {
            m_sampleCoverage.emplace(node, nullptr);
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
                    m_sampleCoverage.erase(childId);
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
        m_sampleCoverage.emplace(nodeId, nullptr);
        bool keepGoing = safeVisit(grg, nodeId);

#if CLEANUP_SAMPLE_SETS_MAPPING
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            // Skip children that aren't part of our search.
            if (m_refCounts[childId] == 0) {
                continue;
            }
            if (--m_refCounts[childId] == 0) {
                m_sampleCoverage.erase(childId);
            }
        }
#endif
        return keepGoing;
    }

    std::unique_ptr<SampleCoverageSet> getSamplesForCandidate(NodeID candidateId) {
        std::unique_ptr<SampleCoverageSet> result;
        auto findIt = m_sampleCoverage.find(candidateId);
        release_assert(findIt != m_sampleCoverage.end());
        result = std::move(findIt->second);
        m_sampleCoverage.erase(findIt);
        return std::move(result);
    }

    std::unordered_map<NodeID, std::unique_ptr<SampleCoverageSet>> m_sampleCoverage;
    std::vector<NodeSamples> m_collectedNodes;

private:
    bool safeVisit(const grgl::GRGPtr& grg, const grgl::NodeID nodeId) {
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

        size_t sampleCount = m_sampleCounts[nodeId];

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

        if (isSample) {
            samplesBeneathSet->addElem(nodeId);
        }

#if CLEANUP_SAMPLE_SETS_MAPPING
        // safe, as NodeID is only handled by a single thread
        m_refCounts[nodeId] = grg->numUpEdges(nodeId);
#endif
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            const auto& childSampleIt = m_sampleCoverage.find(childId);
            if (childSampleIt != m_sampleCoverage.end()) {
                auto& childSamples = (childSampleIt->second);
                if (childSamples->numSamplesCovered() > 1) {
                    candidateNodes.emplace_back(childId);
                }

                CoalescenceTracker tracker{};
                if (computeCoals) {
                    if (dense) {
                        tracker = CoalescenceTracker(coalTrackerSet.get(), &individualCoalCount);
                    } else {
                        tracker = CoalescenceTracker(&coalTrackerMap, childId, &individualCoalCount);
                    }
                }

                samplesBeneathSet->mergeSampleCoverage(*childSamples, tracker);
            }
        }

        // Check if we had a mismatch in expected vs. total sample sets.
        release_assert(nodeId < m_sampleCounts.size());
        release_assert(m_sampleCounts[nodeId] <= grg->numSamples());
        NodeIDSizeT missing = (m_sampleCounts[nodeId] - samplesBeneathSet->numSamplesCovered());
        // We can only record coalescence counts if there are no samples missing.
        if (missing == 0 && computeCoals) {
            grg->setNumIndividualCoals(nodeId, individualCoalCount);
        }

        // If we've reached the root of the graph or have missing samples beneath us, we need to stop the search
        // and emit candidate nodes to map the mutation to.
        const bool keepGoing = (missing == 0 && !isRoot);
        if (missing == 0 && isRoot) {
            m_collectedNodesMutex.lock();
            m_collectedNodes.emplace_back(nodeId, samplesBeneathSet->numSamplesCovered()); // Root is a candidate node.
            m_collectedNodesMutex.unlock();
#if CLEANUP_SAMPLE_SETS_MAPPING
            // Prevent candidates from having their samplesets garbage collected.
            // This is safe, as we only modify the child
            m_refCounts[nodeId] = MAX_GRG_NODES + 1;
#endif
            // hack to allow concurrent writing
            m_sampleCoverage.at(nodeId).swap(samplesBeneathSet);
        } else if (!keepGoing) {
            for (const auto& candidate : candidateNodes) {
                m_collectedNodesMutex.lock();
                m_collectedNodes.emplace_back(candidate, m_sampleCoverage.at(candidate)->numSamplesCovered());
                m_collectedNodesMutex.unlock();
#if CLEANUP_SAMPLE_SETS_MAPPING
                // Prevent candidates from having their samplesets garbage collected.
                // This is safe, as we only modify the child
                m_refCounts[candidate] = MAX_GRG_NODES + 1;
#endif
            }
        } else {
            // Concurrent writing requires that we don't modify the hashmap itself
            // Note that this requires us to have added a nullptr to m_sampleCoverage
            m_sampleCoverage.at(nodeId).swap(samplesBeneathSet);
        }
        return keepGoing;
    }

    std::mutex m_collectedNodesMutex;
    // These are the _total_ samples beneath each node (not restricted to current samples being searched)
    const std::vector<NodeIDSizeT>& m_sampleCounts;
#if CLEANUP_SAMPLE_SETS_MAPPING
    std::vector<std::atomic<NodeIDSizeT>> m_refCounts;
#endif
};

// Tracking individual coalescence is a bit spread out, but I think it is the most efficient way to do it.
// 1. Above, when searching for candidate nodes in the existing hierarchy, any node that does not have its
//    coalescence info computed will be computed and stored.
// 2. Below, when we connect two candidate nodes together we will look for any coalescences between them.
// 3. Below, when we have "left-over" samples that did not have any candidate nodes, we will check for any
//    coalescence within the left-over samples, or between the left-over samples and nodes that have already
//    been covered.

static NodeIDList greedyAddMutation(const MutableGRGPtr& grg,
                                    const std::vector<NodeIDSizeT>& sampleCounts,
                                    const Mutation& newMutation,
                                    const NodeIDList& mutSamples,
                                    MutationMappingStats& stats,
                                    const NodeID shapeNodeIdMax,
                                    std::pair<BpPosition, NodeID>& currentMissing) {
    // The topological order of nodeIDs is maintained through-out this algorithm, because newly added
    // nodes are only ever _root nodes_ (at the time they are added).
    release_assert(grg->nodesAreOrdered());

    const size_t ploidy = grg->getPloidy();
    // The set of nodes that we have covered so far (greedily extended)
    // NodeIDSet covered;
    BitVecCoverageSet coverageSet{grg->numSamples()};

    TopoCandidateCollectorVisitor collector(sampleCounts);
    if (mutSamples.size() > 1) {
        grg->visitTopo(collector, grgl::TraversalDirection::DIRECTION_UP, mutSamples);
    }
    std::vector<NodeSamples>& candidates = collector.m_collectedNodes;
    std::sort(candidates.begin(), candidates.end());
    auto endOfUnique = std::unique(candidates.begin(), candidates.end());
    candidates.erase(endOfUnique, candidates.end());
    std::sort(candidates.begin(), candidates.end(), cmpNodeSamples);

    // The missingness node associated with this mutation. This relies on the fact that
    // missing data is always emitted BEFORE other data for the same site, which is a
    // property of the MutationIterator.
    const NodeID missingnessNode =
        (currentMissing.first == newMutation.getPosition()) ? currentMissing.second : INVALID_NODE_ID;

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
    // Map from an individual to which child contained it.
    // std::unordered_map<NodeIDSizeT, NodeIDSizeT> individualToChild;

    BitVecCoverageSet coalTrackingSet{grg->numSamples()};

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
        const std::unique_ptr<SampleCoverageSet> candidateSet = collector.getSamplesForCandidate(candidateId);
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
            CoalescenceTracker tracker = CoalescenceTracker(&coalTrackingSet, &individualCoalCount);
            coverageSet.mergeSampleCoverage(*candidateSet, tracker);
            /*
            for (const auto sampleId : candidateSet) {
                covered.emplace(sampleId);
                if (ploidy == 2) {
                    auto insertPair = individualToChild.emplace(sampleId / ploidy, candidateId);
                    // The individual already existed from a _different node_, so the two samples will coalesce
                    // at the new mutation node.
                    if (!insertPair.second && candidateId != insertPair.first->second) {
                        individualCoalCount++;
                    }
                }
                }*/

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
    if (uncovered.size() > stats.maxSingletons) {
        stats.maxSingletons = uncovered.size();
    }

    for (auto sampleNodeId : uncovered) {
        grg->connect(mutNodeId, sampleNodeId);
        stats.singletonSampleEdges++;
    }
    // This node needs to be last, for the way we update things.
    addedNodes.push_back(mutNodeId);
    return addedNodes;
}

MutationMappingStats mapMutations(const MutableGRGPtr& grg, MutationIterator& mutations, bool verbose) {
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
    DfsSampleCountVisitor countVisitor;
    fastCompleteDFS(grg, countVisitor);
    std::vector<NodeIDSizeT>& sampleCounts = countVisitor.m_sampleCounts;

    if (verbose) {
        std::cout << "Mapping " << stats.totalMutations << " mutations\n";
    }
    const size_t onePercent = (stats.totalMutations / ONE_HUNDRED_PERCENT) + 1;
    size_t completed = 0;

    // The low-water mark for nodes. If a NodeID is greater than or equal to this, then it
    // is a newly added (mutation) node.
    const NodeID shapeNodeIdMax = grg->numNodes();

    // For each mutation, perform a topological bottom-up traversal from the sample
    // nodes of interest, and collect all nodes that reach a subset of those nodes.
    size_t _ignored = 0;
    std::pair<BpPosition, NodeID> currentMissing = {INVALID_POSITION, INVALID_NODE_ID};
    MutationAndSamples unmapped = {Mutation(0.0, ""), NodeIDList()};
    while (mutations.next(unmapped, _ignored)) {
        const NodeIDList& mutSamples = unmapped.samples;
        if (!mutSamples.empty()) {
            stats.samplesProcessed += mutSamples.size();
            if (mutSamples.size() == 1) {
                stats.mutationsWithOneSample++;
            }
            // Step 1: Add mutation node and create edges to existing nodes.
            NodeIDList addedNodes = greedyAddMutation(
                grg, sampleCounts, unmapped.mutation, mutSamples, stats, shapeNodeIdMax, currentMissing);

            // Step 2: Add the new mutation node to the sample counts.
            sampleCounts.resize(sampleCounts.size() + addedNodes.size());
            for (const auto& nodeId : addedNodes) {
                NodeIDSizeT sumSamples = 0;
                for (const auto& childId : grg->getDownEdges(nodeId)) {
                    release_assert(sampleCounts[childId] > 0);
                    sumSamples += sampleCounts[childId];
                }
                sampleCounts[nodeId] = sumSamples;
            }
        } else {
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
    return stats;
}

}; // namespace grgl
