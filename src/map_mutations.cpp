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
#include "grgl/visitor.h"
#include "util.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>
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
#define USE_DENSE_COVERAGE_PERCENT (0.005f)

class SampleCoverageSet {
public:
    /// the number of total samples
    virtual size_t totalSamples() const = 0;

    /// number of elements in the set, INCLUDING DUPLICATES from below nodes
    virtual size_t numSamplesCovered() const = 0;

    virtual void mergeSampleCoverage(const SampleCoverageSet& other) = 0;

    SampleCoverageSet() = default;
    SampleCoverageSet(SampleCoverageSet&&) = delete;
    SampleCoverageSet(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(const SampleCoverageSet&) = default;
    SampleCoverageSet& operator=(SampleCoverageSet&&) = default;
    virtual ~SampleCoverageSet() = default;
};

class SparseCoverageSet : public SampleCoverageSet {
public:
    explicit SparseCoverageSet(size_t m_totalSampleCount)
        : m_totalSampleCount(m_totalSampleCount) {}

    size_t totalSamples() const override { return m_totalSampleCount; }
    size_t numSamplesCovered() const override { return m_sampleIdxs.size(); }

    /// we can't merge a bitvec coverage set into a
    void mergeSampleCoverage(const SampleCoverageSet& other) override {
        if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            m_sampleIdxs.insert(m_sampleIdxs.end(), other_v->m_sampleIdxs.begin(), other_v->m_sampleIdxs.begin());
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
          m_elems(totalSamples / m_blockSize) {}

    explicit BitVecCoverageSet(const SparseCoverageSet& sparseCoverageSet)
        : BitVecCoverageSet(sparseCoverageSet.totalSamples()) {
        mergeSampleCoverage(sparseCoverageSet);
    }

    void add(NodeID index) {
        size_t vecIdx = index / m_blockSize;
        size_t offset = index % m_blockSize;
        m_elems[vecIdx] |= (1U << offset);
    }

    size_t totalSamples() const override { return m_totalSampleCount; }
    size_t numSamplesCovered() const override { return m_numSampleCoverage; }

    void mergeSampleCoverage(const SampleCoverageSet& other) override {
        m_numSampleCoverage += numSamplesCovered();
        release_assert(totalSamples() == other.totalSamples());
        if (const BitVecCoverageSet* other_v = dynamic_cast<const BitVecCoverageSet*>(&other)) {
            uint64_t* __restrict data = m_elems.data();
            const uint64_t* __restrict otherData = other_v->m_elems.data();
            const size_t len = m_elems.size();

            for (int i = 0; i < len; i++) {
                data[i] |= otherData[i];
            }

        } else if (const SparseCoverageSet* other_v = dynamic_cast<const SparseCoverageSet*>(&other)) {
            for (const uint32_t& sample : other_v->m_sampleIdxs) {
                add(sample);
            }
        }
    }

private:
    constexpr static size_t m_blockSize = 64;
    std::vector<uint64_t> m_elems;
    size_t m_totalSampleCount;
    size_t m_numSampleCoverage{};
};

class TopoCandidateCollectorVisitor : public grgl::GRGVisitor {
public:
    explicit TopoCandidateCollectorVisitor(const std::vector<NodeIDSizeT>& sampleCounts)
        : m_sampleCounts(sampleCounts) {}

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
        const bool isRoot = grg->numUpEdges(nodeId) == 0;
        const bool isSample = grg->isSample(nodeId);
        const size_t ploidy = grg->getPloidy();
        const auto numCoals = grg->getNumIndividualCoals(nodeId);
        const bool computeCoals = !isSample && (ploidy == 2) && (COAL_COUNT_NOT_SET == numCoals);

        size_t individualCoalCount = 0;
        // Map from an individual to which child contained it.
        std::unordered_map<NodeIDSizeT, NodeIDSizeT> individualToChild;
        NodeIDList candidateNodes;

        NodeIDList samplesBeneath;
        if (isSample) {
            samplesBeneath.emplace_back(nodeId);
        }

        size_t childSampleCount{};
        std::vector<const SampleCoverageSet*> childSampleCoverage{};
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            const auto& childSampleIt = m_sampleCoverage.find(childId);
            if (childSampleIt != m_sampleCoverage.end()) {
                auto& sampleCoverage = (childSampleIt->second);
                childSampleCount += sampleCoverage->numSamplesCovered();
                childSampleCoverage.push_back(sampleCoverage.get());
            }
        }

        std::unique_ptr<SampleCoverageSet> samplesBeneathSet;
        double coveragePercent = (double)childSampleCount / (double)grg->numSamples();
        if (coveragePercent > USE_DENSE_COVERAGE_PERCENT) {
            samplesBeneathSet = std::unique_ptr<SampleCoverageSet>(new BitVecCoverageSet(grg->numSamples()));
        } else {
            samplesBeneathSet = std::unique_ptr<SampleCoverageSet>(new SparseCoverageSet(grg->numSamples()));
        }

#if CLEANUP_SAMPLE_SETS_MAPPING
        m_refCounts[nodeId] = grg->numUpEdges(nodeId);
#endif
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            const auto& childSampleIt = m_nodeToSamples.find(childId);
            if (childSampleIt != m_nodeToSamples.end()) {
                auto childSamples = childSampleIt->second;
                if (childSamples.size() > 1) {
                    candidateNodes.emplace_back(childId);
                }
                for (const auto childSampleId : childSamples) {
                    samplesBeneath.emplace_back(childSampleId);
                    if (computeCoals) {
                        auto insertPair = individualToChild.emplace(childSampleId / ploidy, childId);
                        // The individual already existed from a _different child_, so the two samples just coalesced.
                        if (!insertPair.second && childId != insertPair.first->second) {
                            individualCoalCount++;
                        }
                    }
                }
            }
        }
        // Check if we had a mismatch in expected vs. total sample sets.
        release_assert(nodeId < m_sampleCounts.size());
        release_assert(m_sampleCounts[nodeId] <= grg->numSamples());
        NodeIDSizeT missing = (m_sampleCounts[nodeId] - samplesBeneath.size());

        // We can only record coalescence counts if there are no samples missing.
        if (missing == 0 && computeCoals) {
            grg->setNumIndividualCoals(nodeId, individualCoalCount);
        }

        // If we've reached the root of the graph or have missing samples beneath us, we need to stop the search
        // and emit candidate nodes to map the mutation to.
        const bool keepGoing = (missing == 0 && !isRoot);
        if (missing == 0 && isRoot) {
            m_collectedNodes.emplace_back(nodeId, samplesBeneath.size()); // Root is a candidate node.
#if CLEANUP_SAMPLE_SETS_MAPPING
            // Prevent candidates from having their samplesets garbage collected.
            m_refCounts[nodeId] = MAX_GRG_NODES + 1;
#endif
            m_nodeToSamples.emplace(nodeId, std::move(samplesBeneath));
        } else if (!keepGoing) {
            for (const auto& candidate : candidateNodes) {
                m_collectedNodes.emplace_back(candidate, m_nodeToSamples[candidate].size());
#if CLEANUP_SAMPLE_SETS_MAPPING
                // Prevent candidates from having their samplesets garbage collected.
                m_refCounts[candidate] = MAX_GRG_NODES + 1;
#endif
            }
        } else {
            m_nodeToSamples.emplace(nodeId, std::move(samplesBeneath));
        }

#if CLEANUP_SAMPLE_SETS_MAPPING
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            // Skip children that aren't part of our search.
            if (m_refCounts[childId] == 0) {
                continue;
            }
            if (--m_refCounts[childId] == 0) {
                m_nodeToSamples.erase(childId);
            }
        }
#endif
        return keepGoing;
    }

    NodeIDList getSamplesForCandidate(NodeID candidateId) {
        NodeIDList result;
        auto findIt = m_nodeToSamples.find(candidateId);
        release_assert(findIt != m_nodeToSamples.end());
        result = std::move(findIt->second);
        m_nodeToSamples.erase(findIt);
        return std::move(result);
    }

    std::unordered_map<NodeID, std::unique_ptr<SampleCoverageSet>> m_sampleCoverage;
    std::vector<NodeSamples> m_collectedNodes;
    std::unordered_map<NodeID, NodeIDList> m_nodeToSamples;

private:
    // These are the _total_ samples beneath each node (not restricted to current samples being searched)
    const std::vector<NodeIDSizeT>& m_sampleCounts;
#if CLEANUP_SAMPLE_SETS_MAPPING
    std::vector<NodeIDSizeT> m_refCounts;
#endif
};

static bool setsOverlap(const NodeIDSet& alreadyCovered, const NodeIDList& candidateSet) {
    for (auto nodeId : candidateSet) {
        if (alreadyCovered.find(nodeId) != alreadyCovered.end()) {
            return true;
        }
    }
    return false;
}

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
    NodeIDSet covered;

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
    std::unordered_map<NodeIDSizeT, NodeIDSizeT> individualToChild;
    const NodeID mutNodeId = grg->makeNode(1, true);
    if (!newMutation.isMissing()) {
        grg->addMutation(newMutation, mutNodeId, missingnessNode);
    } else {
        currentMissing = {newMutation.getPosition(), mutNodeId};
    }
    NodeIDList addedNodes;
    const size_t numMutSamples = mutSamples.size();
    while (!candidates.empty() && covered.size() < numMutSamples) {
        const auto& candidate = candidates.back();
        const auto candidateId = std::get<0>(candidate);
        const NodeIDList candidateSet = collector.getSamplesForCandidate(candidateId);
        release_assert(!candidateSet.empty());
        // Different candidates may cover different subsets of the sample set that
        // we are currently trying to cover. Those sample sets MUST be non-overlapping
        // or we will introduce a diamond into the graph:
        //  m-->n1-->s0
        //  m-->n2-->s0
        // However, there is no guarantee that there does not exist nodes (n1, n2)
        // that both point to a sample (or samples) that we care about, so we have to
        // track that here. We do that by only considering candidates that have no overlap
        // with our already-covered set.
        if (!setsOverlap(covered, candidateSet)) {
            // Mark all the sample nodes as covered.
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
            }
            if (candidateId >= shapeNodeIdMax) {
                stats.reusedMutNodes++;
            }

            // Use this candidate (or the nodes below it) to cover the sample subset.
            stats.reusedNodes++;
            grg->connect(mutNodeId, candidateId);
            stats.reusedNodeCoverage += candidateSet.size();
            if (candidateSet.size() >= stats.reuseSizeHist.size()) {
                stats.reuseSizeBiggerThanHistMax++;
            } else {
                stats.reuseSizeHist[candidateSet.size()]++;
            }
        }
        candidates.pop_back();
    }

    // Any leftovers, we just connect directly from the new mutation node to the
    // samples.
    NodeIDSet uncovered;
    for (const NodeID sampleNodeId : mutSamples) {
        const auto coveredIt = covered.find(sampleNodeId);
        if (coveredIt == covered.end()) {
            uncovered.emplace(sampleNodeId);
            // The individual had already been seen and >=1 of the samples was previously uncovered,
            // then the new node we create is going to be the coalescence location for that individual.
            if (ploidy == 2) {
                auto insertPair = individualToChild.emplace(sampleNodeId / ploidy, mutNodeId);
                if (!insertPair.second) {
                    individualCoalCount++;
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
