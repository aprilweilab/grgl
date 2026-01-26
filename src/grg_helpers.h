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
#ifndef GRG_HELPERS_H
#define GRG_HELPERS_H

#include "grgl/common.h"
#include "tskit.h"

#include "common_visitors.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/ts2grg.h"
#include "grgl/visitor.h"
#include "util.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <vector>

#include <unistd.h>

#include "grgl/serialize.h"

namespace grgl {

extern char PROCESS_UNIQUE[8];

inline void setProcessUniqueID() {
    release_assert(snprintf(&PROCESS_UNIQUE[0], sizeof(PROCESS_UNIQUE), "%d", getpid()) > 0);
    PROCESS_UNIQUE[7] = 0;
}

inline const char* getProcessUniqueID() { return (const char*)&PROCESS_UNIQUE[0]; }

#define STREAM_PUID "[" << grgl::getProcessUniqueID() << "] "

inline void fastCompleteDFS(const GRGPtr& grg, GRGVisitor& visitor) {
    if (grg->nodesAreOrdered()) {
        for (grgl::NodeID i = 0; i < grg->numNodes(); i++) {
            visitor.visit(grg, i, grgl::TraversalDirection::DIRECTION_DOWN, grgl::DfsPass::DFS_PASS_BACK_AGAIN);
        }
    } else {
        grg->visitDfs(visitor, TraversalDirection::DIRECTION_DOWN, grg->getRootNodes());
    }
}

class DfsSampleCountVisitor : public grgl::GRGVisitor {
public:
    DfsSampleCountVisitor() = default;

    bool
    visit(const GRGPtr& grg, const NodeID nodeId, const TraversalDirection direction, const DfsPass dfsPass) override {
        release_assert(direction == TraversalDirection::DIRECTION_DOWN);
        if (m_sampleCounts.empty()) {
            m_sampleCounts.resize(grg->numNodes());
        }
        if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            NodeIDSizeT count = grg->isSample(nodeId) ? 1 : 0;
            for (const auto& childId : grg->getDownEdges(nodeId)) {
                count += m_sampleCounts[childId];
            }
            m_sampleCounts[nodeId] = count;
        }
        return true;
    }

    std::vector<NodeIDSizeT> m_sampleCounts;
};

static inline size_t getNaiveEdgeCount(const GRGPtr& grg) {
    DfsSampleCountVisitor countVisitor;
    grg->visitDfs(countVisitor, TraversalDirection::DIRECTION_DOWN, grg->getRootNodes());
    const std::vector<NodeIDSizeT>& sampleCounts = countVisitor.m_sampleCounts;
    size_t edgeCount = 0;
    for (const auto& tuple : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
        const NodeID& node = std::get<0>(tuple);
        if (node != INVALID_NODE_ID) {
            edgeCount += sampleCounts.at(node);
        }
        const NodeID& missingnessNode = std::get<2>(tuple);
        if (missingnessNode != INVALID_NODE_ID) {
            edgeCount += sampleCounts.at(missingnessNode);
        }
    }
    return edgeCount;
}

static inline void dumpStats(const GRGPtr& grg, const bool calculateNaiveEdges = false) {
    std::cout << "=== GRG Statistics ===" << std::endl;
    std::cout << "Nodes: " << grg->numNodes() << std::endl;
    std::cout << "Edges: " << grg->numEdges() << std::endl;
    std::cout << "Samples: " << grg->numSamples() << std::endl;
    std::cout << "Mutations: " << grg->getMutations().size() << std::endl;
    std::cout << "Ploidy: " << grg->getPloidy() << std::endl;
    std::cout << "Phased: " << (grg->isPhased() ? "true" : "false") << std::endl;
    std::cout << "Populations: " << grg->getPopulations().size() << std::endl;
    std::cout << "Range of mutations: " << grg->getBPRange().first << " - " << grg->getBPRange().second << std::endl;
    std::cout << "Specified range: " << grg->getSpecifiedBPRange().first << " - " << grg->getSpecifiedBPRange().second
              << std::endl;
    if (calculateNaiveEdges) {
        std::cout << "Naive Edges: " << getNaiveEdgeCount(grg) << std::endl;
    }
    std::cout << "======================" << std::endl;
}

/**
 * Load the GRG so that it can be modified, by adding/removing nodes and edges.
 *
 * @param[in] filename The file to load.
 * @param[in] loadUpEdges Whether to load the "up" edges (in addition to the down edges). Default: true.
 * @return A shared_ptr to the MutableGRG object.
 */
static inline MutableGRGPtr loadMutableGRG(const std::string& filename, const bool loadUpEdges = true) {
    MutableGRGPtr result;
    IFSPointer inStream = std::make_shared<std::ifstream>(filename, std::ios::binary);
    if (!inStream->good()) {
        std::cerr << "Could not read " << filename << std::endl;
        return result;
    }
    try {
        result = readMutableGrg(inStream, loadUpEdges);
    } catch (SerializationFailure& e) {
        std::cerr << "Failed to load GRG: " << e.what() << std::endl;
        return result;
    }
    return result;
}

/**
 * Load the GRG read-only. Edges and nodes cannot be changed, but the Mutations and other auxiliary
 * information still can be.
 *
 * @param[in] filename The file to load.
 * @param[in] loadUpEdges Whether to load the "up" edges (in addition to the down edges). Default: false.
 * @return A shared_ptr to the GRG object.
 */
static inline GRGPtr loadImmutableGRG(const std::string& filename, bool loadUpEdges = false) {
    GRGPtr result;
    IFSPointer inStream = std::make_shared<std::ifstream>(filename, std::ios::binary);
    if (!inStream->good()) {
        std::cerr << "Could not read " << filename << std::endl;
        return result;
    }
    try {
        result = readImmutableGrg(inStream, loadUpEdges);
    } catch (SerializationFailure& e) {
        std::cerr << "Failed to load GRG: " << e.what() << std::endl;
        return result;
    }
    return result;
}

static inline std::pair<NodeIDSizeT, EdgeSizeT>
saveGRG(const GRGPtr& theGRG, const std::string& filename, bool allowSimplify = true) {
    std::ofstream outStream(filename, std::ios::binary);
    return grgl::writeGrg(theGRG, outStream, allowSimplify);
}

static inline bool saveGRGSubset(const GRGPtr& theGRG,
                                 const std::string& filename,
                                 const TraversalDirection direction,
                                 const NodeIDList& seedList,
                                 std::pair<BpPosition, BpPosition> bpRange = {}) {
    api_exc_check(direction == TraversalDirection::DIRECTION_DOWN || theGRG->hasUpEdges(),
                  "Filtering a GRG by samples (individuals) requires loading the GRG with up edges");
    GRGOutputFilter filter(direction, seedList);
    filter.bpRange = bpRange;
    if (seedList.empty()) {
        return false;
    }
    std::ofstream outStream(filename, std::ios::binary);
    grgl::simplifyAndSerialize(theGRG, outStream, filter, true);
    return true;
}

/**
 * Visitor that visits all nodes, hashes their reachable sample sets, and stores a map
 * from every mutation to its sampleset hash.
 */
class MutationSampleHasherVisitor : public TopoSampleSetVisitor {
public:
    explicit MutationSampleHasherVisitor(bool skipRecurrent)
        : m_skipRecurrent(skipRecurrent) {}

    void processNode(const GRGPtr& grg, const NodeIDList& samplesBeneath, const NodeID nodeId) override {
        for (const auto mutId : grg->getMutationsForNode(nodeId)) {
            const auto& mutation = grg->getMutationById(mutId);
            if (m_skipRecurrent && mutation.getAllele() == Mutation::ALLELE_0) {
                continue;
            }
            auto mutIt = m_mutToSamples.emplace(mutation, NodeIDSetOrdered());
            for (const auto& sampleId : samplesBeneath) {
                mutIt.first->second.insert(sampleId);
            }
        }
    }

    std::map<Mutation, NodeIDSetOrdered, MutationLtPosAllele> m_mutToSamples;

private:
    bool m_skipRecurrent;
};

inline NodeIDSetOrdered setDifference(const NodeIDSetOrdered& set1, const NodeIDSetOrdered& set2) {
    NodeIDSetOrdered result;
    for (const auto& item : set1) {
        const auto findIt = set2.find(item);
        if (findIt == set2.end()) {
            result.emplace(item);
        }
    }
    return result;
}

/**
 * Do the two GRGs represent the same sample datasets?
 *
 * True iff:
 * - They have the same set of mutations
 * - They have the same set of samples
 * - The mapping from mutation to sample is the same
 */
static inline bool
equivalentGRGs(const GRGPtr& grg1, const GRGPtr& grg2, std::string& disagreeReason, bool skipRecurrent) {
    if (grg1->numSamples() != grg2->numSamples()) {
        disagreeReason = "Sample counts differ";
        return false;
    }
    auto operationStartTime = std::chrono::high_resolution_clock::now();
    MutationSampleHasherVisitor visitor1(skipRecurrent);
    fastCompleteDFS(grg1, visitor1);
    auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - operationStartTime)
            .count();
    std::cout << "Done collecting from first GRG: " << elapsed << " seconds" << std::endl;
    visitor1.clearSampleSets();

    operationStartTime = std::chrono::high_resolution_clock::now();
    MutationSampleHasherVisitor visitor2(skipRecurrent);
    fastCompleteDFS(grg2, visitor2);
    elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - operationStartTime)
            .count();
    std::cout << "Done collecting from second GRG: " << elapsed << " seconds" << std::endl;
    visitor2.clearSampleSets();

    std::stringstream errorStream;
    errorStream << std::fixed << std::setprecision(2);
    bool failed = false;
    if (visitor1.m_mutToSamples.size() != visitor2.m_mutToSamples.size()) {
        errorStream << "Reachable mutation counts differ: " << visitor1.m_mutToSamples.size() << " vs. "
                    << visitor2.m_mutToSamples.size() << std::endl;
        failed = true;
    }

    for (const auto& pair : visitor1.m_mutToSamples) {
        const auto& mutation = pair.first;
        const auto& samplesReached = pair.second;
        const auto& otherIt = visitor2.m_mutToSamples.find(mutation);
        if (otherIt == visitor2.m_mutToSamples.end()) {
            errorStream << "Missing mutation in the second GRG: " << mutation.getPosition() << ","
                        << mutation.getAllele() << std::endl;
            ;
            errorStream << "Has " << samplesReached.size() << " samples" << std::endl;
            failed = true;
        } else if (otherIt->second != samplesReached) {
            errorStream << "Mutation (" << mutation.getPosition() << "," << mutation.getAllele()
                        << ") has different reached samples" << std::endl;
            errorStream << "Samples in first but not second: ";
            for (const auto& sampleId : setDifference(samplesReached, otherIt->second)) {
                errorStream << sampleId << ", ";
            }
            errorStream << std::endl;
            errorStream << "Samples in second but not first: ";
            for (const auto& sampleId : setDifference(otherIt->second, samplesReached)) {
                errorStream << sampleId << ", ";
            }
            errorStream << std::endl;
            failed = true;
        }
    }
    disagreeReason = errorStream.str();
    return !failed;
}

inline MutableGRGPtr grgFromTrees(const std::string& filename,
                                  bool binaryMutations = false,
                                  bool useNodeTimes = false,
                                  bool maintainTopology = false,
                                  bool computeCoals = false) {
    tsk_treeseq_t treeSeq;
    if (0 != tsk_treeseq_load(&treeSeq, filename.c_str(), 0)) {
        throw TskitApiFailure("Failed to load treeseq file");
    }

    return grgl::convertTreeSeqToGRG(&treeSeq, binaryMutations, useNodeTimes, maintainTopology, computeCoals);
}

/**
 * Common data container used internally when building/modifying graphs. During actual construction (not modification)
 * we want to manage the life/death of node data ourselves.
 */
template <typename T> class UnmanagedNodeData {
public:
    UnmanagedNodeData() = default;

    // Direction is the direction we are currently moving. For DFS, the second pass is reversed from
    // the visitor direction.
    explicit UnmanagedNodeData(NodeIDSizeT numNodes, TraversalDirection direction)
        : m_data(numNodes),
          m_direction(direction) {}

    bool uninitialized() const { return m_data.empty(); }

    void add(const GRGPtr& grg, const NodeID node, std::unique_ptr<T> data) {
        if (node >= m_data.size()) {
            m_data.resize(node + 1);
        }
        m_data[node] = std::move(data);
    }

    const std::unique_ptr<T>& get(const NodeID node) const {
        static std::unique_ptr<T> nullT;
        if (node < m_data.size()) {
            return m_data[node];
        }
        return nullT;
    }

    void destroy(const NodeID node) {
        if (node < m_data.size()) {
            m_data[node].reset();
        }
    }

    void decr(const NodeID node) {}

    void decr_children(const NodeIDList& children) { (void)children; }

    void decr_parents(const NodeIDList& parents) { (void)parents; }

    bool hasNode(const NodeID nodeId) const { return nodeId < m_data.size(); }

    std::vector<std::unique_ptr<T>> m_data;

private:
    TraversalDirection m_direction{TraversalDirection::DIRECTION_DOWN};
};

/**
 * Common data container used internally when building/modifying graphs. We have some (storage-expensive) data
 * we want to store for each node during a traveral. When that data is no longer needed we want to release
 * it ASAP. This uses a hash table from NodeID to data (type T) to store the data and reference counts each
 * node's data based on the number of parents/children processed so far.
 */
template <typename T> class RefCountedNodeData {
public:
    RefCountedNodeData() = default;

    // Direction is the direction we are currently moving. For DFS, the second pass is reversed from
    // the visitor direction.
    explicit RefCountedNodeData(NodeIDSizeT numNodes, TraversalDirection direction)
        : m_data(numNodes),
          m_refCount(numNodes, 0),
          m_direction(direction) {}

    bool uninitialized() const { return m_refCount.empty(); }

    void add(const GRGPtr& grg, const NodeID node, std::unique_ptr<T> data) {
        if (node >= m_refCount.size()) {
            m_refCount.resize(node + 1);
            m_data.resize(node + 1);
        }
        if (m_direction == TraversalDirection::DIRECTION_UP) {
            m_refCount[node] = grg->numUpEdges(node);
        } else {
            m_refCount[node] = grg->numDownEdges(node);
        }
        // If no one needs this data, then don't store it!
        if (m_refCount[node] > 0) {
            m_data[node] = std::move(data);
        }
    }

    const std::unique_ptr<T>& get(const NodeID node) const {
        static std::unique_ptr<T> nullT;
        if (node < m_data.size()) {
            return m_data[node];
        }
        return nullT;
    }

    void destroy(const NodeID node) {
        m_refCount.at(node) = 0;
        m_data[node].reset();
    }

    void decr(const NodeID node) {
        release_assert(m_refCount.at(node) > 0);
        m_refCount[node]--;
        if (m_refCount[node] == 0) {
            m_data[node].reset();
        }
    }

    void decr_children(const NodeIDList& children) {
        release_assert(m_direction == TraversalDirection::DIRECTION_UP);
        for (const auto& child : children) {
            decr(child);
        }
    }

    void decr_parents(const NodeIDList& parents) {
        release_assert(m_direction == TraversalDirection::DIRECTION_DOWN);
        for (const auto& parent : parents) {
            decr(parent);
        }
    }

    bool hasNode(const NodeID nodeId) const { return nodeId < m_data.size(); }

    std::vector<std::unique_ptr<T>> m_data;

private:
    std::vector<NodeIDSizeT> m_refCount;
    TraversalDirection m_direction{TraversalDirection::DIRECTION_DOWN};
};

template <typename T> inline std::unique_ptr<T> make_unique(T* ptr) { return std::unique_ptr<T>(ptr); }

inline void mergeInPlace(NodeIDList::iterator start, NodeIDList::iterator middle, NodeIDList::iterator end) {
    std::inplace_merge(start, middle, end);
}

/**
 * Common coalescence calculations for a parent node with children.
 * Coalescence information lets us _exactly_ compute the mean/variance for each mutation
 * in the diploid genotype matrix. Anytime the graph is modify (including during construction)
 * we need to track this information.
 *
 * This only applies to diploid data. Each "individual ID" is just the sample ID divided by 2.
 *
 * @param[in] grg The graph.
 * @param[in] nodeToIndivs Either an UnmanagedNodeData<NodeIDList> or RefCountedNodeData<NodeIDList>
 *     that is tracking all of the uncoalesced individual information per node. This is not as memory
 *     expensive as most per-node data, because once an individual coalesces all parent nodes do not
 *     track that information. The lists of individuals are kept in vectors (NodeIDList) and sorted,
 *     as this is more than 2x faster than uses a hash table to track coalescences.
 * @param[in] nodeId The node being computed (parent).
 * @param[in] children The child nodes of nodeId.
 * @param[in,out] uncoalescedIndivs Sorted list of individuals that have not yet coalesced. All the
 *     uncoalesced individuals at the children nodes (passed in) are evaluated and any that coalesce
 *     are left out of this list. Any that remain uncoalesced remain. The list can be passed in to
 *     this function empty (the usual use case), or with individuals that you want to check coalescence
 *     against.
 * @param[in] implicitSamples Set to true if samples should be treated implicitly instead of storing
 *     them in nodeToIndivs, which is only valid when performing a full graph traversal. Partial
 *     traversal require this to be false, and samples should be handled explicitly. Default: false.
 */
template <typename NodeDataT = RefCountedNodeData<NodeIDList>>
inline NodeIDSizeT getCoalsForParent(const GRGPtr& grg,
                                     NodeDataT& nodeToIndivs,
                                     const NodeID nodeId,
                                     const NodeIDList& children,
                                     NodeIDList& uncoalescedIndivs,
                                     bool implicitSamples = false) {
    constexpr NodeIDSizeT ploidy = 2;
    release_assert(grg->getPloidy() == ploidy);
    NodeIDSizeT coalCount = 0;

    // Collect all "individuals below" each child and whenever we see one twice, count
    // it as a coalescence and remove it from the list of seen individuals.
    if (!implicitSamples && grg->isSample(nodeId)) {
        uncoalescedIndivs.emplace_back(nodeId / ploidy);
    } else {
        for (const auto& childId : children) {
            const ssize_t sortedEndPos = static_cast<ssize_t>(uncoalescedIndivs.size());
            if (implicitSamples && grg->isSample(childId)) {
                uncoalescedIndivs.emplace_back(childId / ploidy);
            } else {
                const auto& childIndivs = nodeToIndivs.get(childId);
                if (childIndivs) {
                    for (const NodeID childIndivId : *childIndivs) {
                        uncoalescedIndivs.emplace_back(childIndivId);
                    }
                    nodeToIndivs.decr(childId);
                }
            }
            if (sortedEndPos > 0) {
                mergeInPlace(uncoalescedIndivs.begin(),
                             std::next(uncoalescedIndivs.begin(), sortedEndPos),
                             uncoalescedIndivs.end());
            }
        }
        if (!uncoalescedIndivs.empty()) {
            size_t j = 0;
            NodeID prevIndiv = INVALID_NODE_ID;
            for (size_t i = 0; i < uncoalescedIndivs.size(); i++) {
                const NodeID currentIndiv = uncoalescedIndivs[i];
                if (currentIndiv == prevIndiv) {
                    release_assert(j > 0);
                    j--; // Erase the previous and skip the current
                    coalCount++;
                    prevIndiv = INVALID_NODE_ID;
                } else {
                    if (j < i) {
                        uncoalescedIndivs[j] = uncoalescedIndivs[i];
                    }
                    j++;
                    prevIndiv = currentIndiv;
                }
            }
            uncoalescedIndivs.resize(j);
        }
    }
    return coalCount;
}

// Flip a traversal direction.
inline TraversalDirection flipDir(TraversalDirection dir) {
    return (dir == TraversalDirection::DIRECTION_UP) ? TraversalDirection::DIRECTION_DOWN
                                                     : TraversalDirection::DIRECTION_UP;
}

// Convert NodeIDSet to NodeIDList
inline NodeIDList nodeSetToList(const NodeIDSet& set) {
    NodeIDList result;
    for (const auto& item : set) {
        result.emplace_back(item);
    }
    return std::move(result);
}

inline NodeIDList complementSampleList(const NodeIDSizeT numSamples, const NodeIDList& samples) {
    NodeIDList result;
    NodeID previous = INVALID_NODE_ID;
    size_t j = 0;
    for (NodeID nodeId = 0; nodeId < numSamples; nodeId++) {
        bool matched = false;
        if (j < samples.size()) {
            NodeID current = samples[j];
            api_exc_check(previous == INVALID_NODE_ID || previous < current,
                          "Sample list must be in sorted order with no duplicates");
            release_assert(current >= nodeId);
            if (nodeId == current) {
                j++;
                previous = current;
                matched = true;
            }
        }
        if (!matched) {
            result.emplace_back(nodeId);
        }
    }
    return std::move(result);
}

} // namespace grgl

#endif /* GRG_HELPERS_H */
