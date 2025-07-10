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
 * Visitor that visits all nodes, hashes their child nodeIDs, and stores a map
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

    DigestToNode m_hashToNodeId;
    size_t m_duplicates{};
};

/**
 * Visitor that visits all nodes, hashes their reachable sample sets, and tries to map
 * them to another GRG that has already had its nodes hashed.
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
        // Previousl we merged based on covered samples, now we merge based on (immediate) children nodes.
        // The difference is that we previously merged more nodes, but now we maintain more hierarchy. In
        //
        // There may be post-processing steps we can do later to "cleanup" duplication that involves hierarchy
        // under one node, and no hierarchy under the other, but both cover the same sample set.
        if (grg->isSample(nodeId)) {
            for (const auto mutId : grg->getMutationsForNode(nodeId)) {
                const auto& mutation = grg->getMutationById(mutId);
                m_targetGrg.addMutation(mutation, nodeId);
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
            for (const auto mutId : grg->getMutationsForNode(nodeId)) {
                const auto& mutation = grg->getMutationById(mutId);
                m_targetGrg.addMutation(mutation, targetNodeId);
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
            for (const auto mutId : grg->getMutationsForNode(nodeId)) {
                const auto& mutation = grg->getMutationById(mutId);
                m_targetGrg.addMutation(mutation, targetNodeId);
            }
        }
        return true;
    }

    size_t m_mappedExactly{};

private:
    MutableGRG& m_targetGrg;
    std::vector<NodeID> m_nodeIdToTargetNodeId;
    DigestToNode& m_targetHashToNodeId;
    bool m_combineNodes;
};

void MutableGRG::merge(const std::list<std::string>& otherGrgFiles, bool combineNodes) {
    // Compute hashes on the GRG that we are mapping to.
    NodeHasherVisitor hashVisitor;
    if (combineNodes) {
        this->visitDfs(hashVisitor, TraversalDirection::DIRECTION_DOWN, this->getRootNodes());
    }
    DigestToNode& hashToNodeId = hashVisitor.m_hashToNodeId;

    for (const auto& otherGrgFile : otherGrgFiles) {
        const auto otherGrg = loadImmutableGRG(otherGrgFile);
        if (!otherGrg) {
            std::stringstream err;
            err << "Could not load GRG from " << otherGrgFile;
            throw ApiMisuseFailure(err.str().c_str());
        }
        if (this->numSamples() != otherGrg->numSamples()) {
            std::stringstream err;
            err << "Sample count mismatch: " << this->numSamples() << " vs. " << otherGrg->numSamples();
            throw ApiMisuseFailure(err.str().c_str());
        }
        // Do the actual node mapping and copy relevant nodes/mutations/edges to the target GRG.
        NodeMapperVisitor mapperVisitor(otherGrg, *this, hashToNodeId, combineNodes);
        fastCompleteDFS(otherGrg, mapperVisitor);
        std::cout << "Mapped exactly: " << mapperVisitor.m_mappedExactly << std::endl;
        // Copy any left-over mutations that were not associated with nodes
        for (const auto mutId : otherGrg->getUnmappedMutations()) {
            const auto& mutation = otherGrg->getMutationById(mutId);
            this->addMutation(mutation, INVALID_NODE_ID);
        }
        // The range has to be contiguous, even if there is a "gap" between the two merged GRGs.
        const auto& otherRange = otherGrg->getSpecifiedBPRange();
        if (otherRange.first < this->m_specifiedRange.first) {
            this->m_specifiedRange.first = otherRange.first;
        }
        if (otherRange.second > this->m_specifiedRange.second) {
            this->m_specifiedRange.second = otherRange.second;
        }
    }
}

} // namespace grgl
