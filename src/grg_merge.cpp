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
 * Constructs a map from hash(reached-samples) to nodeId, for fast lookup of exactly
 * matching nodes.
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
 * Visitor that visits all nodes, hashes their reachable sample sets, and stores a map
 * from that (unique) hash to the nodeId.
 */
class NodeHasherVisitor : public TopoSampleHashVisitor {
public:
    void processNode(const GRGPtr& grg, const HashDigest& digest, const NodeID nodeId) override {
        if (!grg->isSample(nodeId)) {
            m_hashToNodeId.emplace(digest, nodeId);
        }
    }

    DigestToNode m_hashToNodeId;
};

/**
 * Visitor that visits all nodes, hashes their reachable sample sets, and tries to map
 * them to another GRG that has already had its nodes hashed.
 */
class NodeMapperVisitor : public TopoSampleHashVisitor {
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

    void processNode(const GRGPtr& grg, const HashDigest& digest, const NodeID nodeId) override {
        if (grg->isSample(nodeId)) {
            for (const auto mutId : grg->getMutationsForNode(nodeId)) {
                const auto& mutation = grg->getMutationById(mutId);
                m_targetGrg.addMutation(mutation, nodeId);
            }
            m_mappedExactly++;
            return;
        }
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
            for (const auto mutId : grg->getMutationsForNode(nodeId)) {
                const auto& mutation = grg->getMutationById(mutId);
                m_targetGrg.addMutation(mutation, targetNodeId);
            }
        }
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
        hashVisitor.clearSampleSets();
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
