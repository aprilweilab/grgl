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
#include <set>
#include <sstream>
#include <vector>

#include "grgl/serialize.h"

namespace grgl {

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
    for (const auto& pair : grg->getNodeMutationPairs()) {
        if (pair.first != INVALID_NODE_ID) {
            edgeCount += sampleCounts.at(pair.first);
        }
    }
    return edgeCount;
}

static inline void dumpStats(const GRGPtr& grg) {
    std::cout << "=== GRG Statistics ===" << std::endl;
    std::cout << "Nodes: " << grg->numNodes() << std::endl;
    std::cout << "Edges: " << grg->numEdges() << std::endl;
    std::cout << "Samples: " << grg->numSamples() << std::endl;
    std::cout << "Mutations: " << grg->getMutations().size() << std::endl;
    std::cout << "Naive Edges: " << getNaiveEdgeCount(grg) << std::endl;
    std::cout << "Ploidy: " << grg->getPloidy() << std::endl;
    std::cout << "Populations: " << grg->getPopulations().size() << std::endl;
    std::cout << "Range of mutations: " << grg->getBPRange().first << " - " << grg->getBPRange().second << std::endl;
    std::cout << "Specified range: " << grg->getSpecifiedBPRange().first << " - " << grg->getSpecifiedBPRange().second
              << std::endl;
    std::cout << "======================" << std::endl;
}

static inline MutableGRGPtr loadMutableGRG(const std::string& filename) {
    MutableGRGPtr result;
    IFSPointer inStream = std::make_shared<std::ifstream>(filename, std::ios::binary);
    if (!inStream->good()) {
        std::cerr << "Could not read " << filename << std::endl;
        return result;
    }
    try {
        result = readMutableGrg(inStream);
    } catch (SerializationFailure& e) {
        std::cerr << "Failed to load GRG: " << e.what() << std::endl;
        return result;
    }
    return result;
}

static inline GRGPtr loadImmutableGRG(const std::string& filename, bool loadUpEdges = true) {
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

static inline std::pair<NodeIDSizeT, NodeIDSizeT>
saveGRG(const GRGPtr& theGRG, const std::string& filename, bool allowSimplify = true) {
    std::ofstream outStream(filename, std::ios::binary);
    return grgl::writeGrg(theGRG, outStream, allowSimplify);
}

static inline void saveGRGSubset(const GRGPtr& theGRG,
                                 const std::string& filename,
                                 const TraversalDirection direction,
                                 const NodeIDList& seedList,
                                 std::pair<BpPosition, BpPosition> bpRange = {}) {
    std::ofstream outStream(filename, std::ios::binary);
    GRGOutputFilter filter(direction, seedList);
    filter.bpRange = bpRange;
    grgl::simplifyAndSerialize(theGRG, outStream, filter, true);
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

} // namespace grgl

#endif /* GRG_HELPERS_H */
