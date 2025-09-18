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
#include "common_visitors.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"
#include "grgl/visitor.h"
#include "node_unique_hash.h"
#include "util.h"

#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace grgl {

using DigestToNode = std::unordered_map<HashDigest, NodeID>;

/**
 * Abstract class used by the sample-set based algorithms (slower).
 */
class TopoSampleHashVisitor : public TopoSampleSetVisitor {
public:
    TopoSampleHashVisitor() = default;

    virtual void processNode(const GRGPtr& grg, const HashDigest& digest, NodeID nodeId) = 0;

    void processNode(const GRGPtr& grg, const NodeIDList& samplesBeneath, NodeID nodeId) override {
        const HashDigest hash = !grg->isSample(nodeId) ? hashNodeSet(samplesBeneath) : "";
        processNode(grg, hash, nodeId);
    }
};

/**
 * Option 1: Visitor that visits all nodes, hashes the samples beneath them, and creates a map
 * from that (unique) hash to the nodeId.
 */
class SampleHasherVisitor : public TopoSampleHashVisitor {
public:
    SampleHasherVisitor() = default;

    void processNode(const GRGPtr& grg, const HashDigest& digest, const NodeID nodeId) override {
        if (!grg->isSample(nodeId)) {
            m_hashToNodeId.emplace(digest, nodeId);
        }
    }

    void clearTemporaryData() { this->clearSampleSets(); }

    DigestToNode m_hashToNodeId;
};

/**
 * Option 2: Visitor that visits all nodes, hashes their child nodeIDs, and stores a map
 * from that (unique) hash to the (parent) nodeId.
 */
class NodeHasherVisitor : public GRGVisitor {
public:
    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (dfsPass == DfsPass::DFS_PASS_NONE) {
            release_assert(direction == TraversalDirection::DIRECTION_UP);
        } else if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            release_assert(direction == TraversalDirection::DIRECTION_DOWN);
        } else {
            return true;
        }
        // Sample IDs must match between merged GRGs, no reason to hash them.
        if (!grg->isSample(nodeId)) {
            NodeIDList children = grg->getDownEdges(nodeId);
            if (!grg->edgesAreOrdered()) {
                std::sort(children.begin(), children.end());
            }
            HashDigest digest = hashNodeSet(children);
            const auto inserted = m_hashToNodeId.emplace(std::move(digest), nodeId);
            // Duplicates are possible, and ok. It just means two nodes have the exact same
            // children, and we just pick one arbitrarily to use.
            if (!inserted.second) {
                m_duplicates++;
            }
        }
        return true;
    }

    void clearTemporaryData() {}

    DigestToNode m_hashToNodeId;
    size_t m_duplicates{};
};

/**
 * Option 1: Maps nodes using unique hashes based on the set of samples reachable from the node.
 * This is more RAM/CPU intensive, and reduces the amount of hierarchy in the graph, but produces
 * smaller graphs.
 */
class NodeMapperVisitorSamples : public TopoSampleHashVisitor {
public:
    NodeMapperVisitorSamples(const GRGPtr& sourceGrg,
                             MutableGRG& targetGrg,
                             DigestToNode& targetHashToNodeId,
                             bool combineNodes)
        : m_targetGrg(targetGrg),
          m_targetHashToNodeId(targetHashToNodeId),
          m_combineNodes(combineNodes) {
        m_nodeIdToTargetNodeId.resize(sourceGrg->numNodes(), INVALID_NODE_ID);
    }

    void processNode(const GRGPtr& grg, const HashDigest& digest, const NodeID nodeId) override {
        if (grg->isSample(nodeId)) {
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, nodeId, mutAndMiss.second);
            }
            m_mappedExactly++;
            return;
        }
        // See if this node maps exactly to a node in the target GRG
        const auto mappedNodeIt = m_targetHashToNodeId.find(digest);
        if (m_combineNodes && mappedNodeIt != m_targetHashToNodeId.end()) {
            const auto targetNodeId = mappedNodeIt->second;
            // Copy all the mutations for this node to the target one.
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, targetNodeId, mutAndMiss.second);
            }
            release_assert(nodeId < m_nodeIdToTargetNodeId.size());
            m_nodeIdToTargetNodeId[nodeId] = targetNodeId;
            m_mappedExactly++;
        } else {
            // We have to make a node in the target GRG.
            const NodeID targetNodeId = m_targetGrg.makeNode();
            release_assert(nodeId < m_nodeIdToTargetNodeId.size());
            m_nodeIdToTargetNodeId[nodeId] = targetNodeId;
            m_targetHashToNodeId[digest] = targetNodeId;
            // Copy node data.
            m_targetGrg.setNumIndividualCoals(targetNodeId, grg->getNumIndividualCoals(nodeId));

            // And reconnect all the child nodes appropriately. We only do child because
            // we are doing a bottom-up topological graph search.
            for (const auto& childId : grg->getDownEdges(nodeId)) {
                if (grg->isSample(childId)) {
                    m_targetGrg.connect(targetNodeId, childId);
                } else {
                    // The topological order assures that we _must_ have a mapping already.
                    const auto targetChildId = m_nodeIdToTargetNodeId[childId];
                    release_assert(INVALID_NODE_ID != targetChildId);

                    m_targetGrg.connect(targetNodeId, targetChildId);
                }
            }

            // Copy all the mutations for this node to the target one.
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, targetNodeId, mutAndMiss.second);
            }
        }
    }

    size_t m_mappedExactly{};
    std::vector<NodeID> m_nodeIdToTargetNodeId;

private:
    MutableGRG& m_targetGrg;
    DigestToNode& m_targetHashToNodeId;
    bool m_combineNodes;
};

/**
 * Option 2: Maps nodes using unique hashes based on the set of immediate children beneath the node.
 * This is RAM/CPU efficient, maintains more hierachy, but produces 5-8% larger graphs.
 */
class NodeMapperVisitor : public GRGVisitor {
public:
    NodeMapperVisitor(const GRGPtr& sourceGrg,
                      MutableGRG& targetGrg,
                      DigestToNode& targetHashToNodeId,
                      bool combineNodes)
        : m_targetGrg(targetGrg),
          m_targetHashToNodeId(targetHashToNodeId),
          m_combineNodes(combineNodes) {
        m_nodeIdToTargetNodeId.resize(sourceGrg->numNodes(), INVALID_NODE_ID);
    }

    // Returns true if the result is in sorted order, false otherwise.
    bool convertToTargets(NodeIDList& nodeList) {
        bool isSorted = true;
        NodeID prevValue = 0;
        for (size_t i = 0; i < nodeList.size(); i++) {
            nodeList[i] = m_nodeIdToTargetNodeId[nodeList[i]];
            release_assert(INVALID_NODE_ID != nodeList[i]);
            if (nodeList[i] < prevValue) {
                isSorted = false;
            }
            prevValue = nodeList[i];
        }
        return isSorted;
    }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (dfsPass == DfsPass::DFS_PASS_NONE) {
            release_assert(direction == TraversalDirection::DIRECTION_UP);
        } else if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            release_assert(direction == TraversalDirection::DIRECTION_DOWN);
        } else {
            return true;
        }
        // Previously we merged based on covered samples, now we merge based on (immediate) children nodes.
        // The difference is that we previously merged more nodes, but now we maintain more hierarchy. In
        //
        // There may be post-processing steps we can do later to "cleanup" duplication that involves hierarchy
        // under one node, and no hierarchy under the other, but both cover the same sample set.
        if (grg->isSample(nodeId)) {
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, nodeId, mutAndMiss.second);
            }
            m_nodeIdToTargetNodeId[nodeId] = nodeId;
            m_mappedExactly++;
            return true;
        }
        NodeIDList children = grg->getDownEdges(nodeId);
        const bool alreadySorted = convertToTargets(children);
        if (!alreadySorted) {
            std::sort(children.begin(), children.end());
        }
        const HashDigest digest = hashNodeSet(children);

        // See if this node maps exactly to a node in the target GRG
        const auto mappedNodeIt = m_targetHashToNodeId.find(digest);
        if (m_combineNodes && mappedNodeIt != m_targetHashToNodeId.end()) {
            const auto targetNodeId = mappedNodeIt->second;
            // Copy all the mutations for this node to the target one.
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, targetNodeId, mutAndMiss.second);
            }
            release_assert(nodeId < m_nodeIdToTargetNodeId.size());
            m_nodeIdToTargetNodeId[nodeId] = targetNodeId;
            m_mappedExactly++;
        } else {
            // We have to make a node in the target GRG.
            const NodeID targetNodeId = m_targetGrg.makeNode();
            release_assert(nodeId < m_nodeIdToTargetNodeId.size());
            m_nodeIdToTargetNodeId[nodeId] = targetNodeId;
            const auto inserted = m_targetHashToNodeId.emplace(digest, targetNodeId);
            release_assert(inserted.second); // Every hash must be unique.
            // Copy node data.
            m_targetGrg.setNumIndividualCoals(targetNodeId, grg->getNumIndividualCoals(nodeId));

            // And reconnect all the child nodes appropriately. We only do child because
            // we are doing a bottom-up topological graph search.
            for (const auto& childId : children) {
                m_targetGrg.connect(targetNodeId, childId);
            }

            // Copy all the mutations for this node to the target one.
            for (const auto mutAndMiss : grg->getMutationsForNode<GRG::MutAndNode>(nodeId)) {
                const auto& mutation = grg->getMutationById(mutAndMiss.first);
                m_targetGrg.addMutation(mutation, targetNodeId, mutAndMiss.second);
            }
        }
        return true;
    }

    size_t m_mappedExactly{};
    std::vector<NodeID> m_nodeIdToTargetNodeId;

private:
    MutableGRG& m_targetGrg;
    DigestToNode& m_targetHashToNodeId;
    bool m_combineNodes;
};

template <typename Hasher, typename Mapper>
void mergeHelper(MutableGRG& grg, const std::list<std::string>& otherGrgFiles, bool combineNodes, bool verbose) {
    // Compute hashes on the GRG that we are mapping to.
    Hasher hashVisitor;
    if (combineNodes) {
        grg.visitDfs(hashVisitor, TraversalDirection::DIRECTION_DOWN, grg.getRootNodes());
        hashVisitor.clearTemporaryData();
    }
    DigestToNode& hashToNodeId = hashVisitor.m_hashToNodeId;

    size_t seenMutations = grg.numMutations();
    for (const auto& otherGrgFile : otherGrgFiles) {
        const auto otherGrg = loadImmutableGRG(otherGrgFile, /*loadUpEdges=*/false);
        if (!otherGrg) {
            std::stringstream err;
            err << "Could not load GRG from " << otherGrgFile;
            throw ApiMisuseFailure(err.str().c_str());
        }
        if (grg.numSamples() != otherGrg->numSamples()) {
            std::stringstream err;
            err << "Sample count mismatch: " << grg.numSamples() << " vs. " << otherGrg->numSamples();
            throw ApiMisuseFailure(err.str().c_str());
        }
        seenMutations += otherGrg->numMutations();

        // We need to adjust these mutations (their missingness nodes) after we map them into the target.
        const size_t adjustMutsAfter = grg.numMutations();

        // Do the actual node mapping and copy relevant nodes/mutations/edges to the target GRG.
        Mapper mapperVisitor(otherGrg, grg, hashToNodeId, combineNodes);
        fastCompleteDFS(otherGrg, mapperVisitor);
        if (verbose) {
            std::cout << "Mapped exactly: " << mapperVisitor.m_mappedExactly << std::endl;
        }
        // Copy any left-over mutations that were not associated with nodes
        for (const auto mutAndMiss : otherGrg->getUnmappedMutations<GRG::MutAndNode>()) {
            const auto& mutation = otherGrg->getMutationById(mutAndMiss.first);
            grg.addMutation(mutation, INVALID_NODE_ID, mutAndMiss.second);
        }
        // Adjust all the missingness node ID's to use the new mapping.
        grg.adjustMissingnessNodeIds(adjustMutsAfter, mapperVisitor.m_nodeIdToTargetNodeId);

        // The range has to be contiguous, even if there is a "gap" between the two merged GRGs.
        const auto& otherRange = otherGrg->getSpecifiedBPRange();
        auto myRange = grg.getSpecifiedBPRange();
        if (otherRange.first < myRange.first) {
            myRange.first = otherRange.first;
        }
        if (otherRange.second > myRange.second) {
            myRange.second = otherRange.second;
        }
        grg.setSpecifiedBPRange(myRange);
    }
    release_assert(grg.numMutations() == seenMutations);
}

void MutableGRG::merge(const std::list<std::string>& otherGrgFiles,
                       bool combineNodes,
                       bool useSampleSets,
                       bool verbose) {
    if (useSampleSets) {
        mergeHelper<SampleHasherVisitor, NodeMapperVisitorSamples>(*this, otherGrgFiles, combineNodes, verbose);
    } else {
        mergeHelper<NodeHasherVisitor, NodeMapperVisitor>(*this, otherGrgFiles, combineNodes, verbose);
    }
}

} // namespace grgl
