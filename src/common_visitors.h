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
#ifndef GRG_COMMON_VISITORS_H
#define GRG_COMMON_VISITORS_H

#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/visitor.h"
#include "util.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace grgl {

// There are likely other memory improvements that would have more impact.
#define CLEANUP_SAMPLE_SETS

/**
 * Constructs a map from hash(reached-samples) to nodeId, for fast lookup of exactly
 * matching nodes.
 */
class TopoSampleSetVisitor : public GRGVisitor {
public:
    TopoSampleSetVisitor() = default;

    virtual void processNode(const grgl::GRGPtr& grg, const NodeIDList& samplesBeneath, NodeID nodeId) = 0;

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        // DFS "down" is equivalent to Topo "up", for our concerns. The DFS is faster if you want to
        // visit the whole graph. The topo is needed if you only want to visit a subset of the graph.
        if (dfsPass == DfsPass::DFS_PASS_NONE) {
            release_assert(direction == TraversalDirection::DIRECTION_UP);
        } else if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            release_assert(direction == TraversalDirection::DIRECTION_DOWN);
        } else {
            return true;
        }
        if (m_sampleLists.empty()) {
            m_sampleLists.resize(grg->numNodes());
#ifdef CLEANUP_SAMPLE_SETS
            release_assert(m_refCounts.empty());
            m_refCounts.resize(grg->numNodes());
#endif
        }
        release_assert(nodeId <= m_sampleLists.size());
#ifdef CLEANUP_SAMPLE_SETS
        m_refCounts[nodeId] = grg->numUpEdges(nodeId);
#endif
        auto& samplesForThisNode = m_sampleLists[nodeId];
        for (const auto& childId : grg->getDownEdges(nodeId)) {
            if (grg->isSample(childId)) {
                samplesForThisNode.push_back(childId);
            } else {
                for (const auto& sampleId : m_sampleLists[childId]) {
                    samplesForThisNode.push_back(sampleId);
                }
            }
#ifdef CLEANUP_SAMPLE_SETS
            // Memory cleanup.
            release_assert(m_refCounts[childId] > 0);
            if (--m_refCounts[childId] == 0) {
                m_sampleLists[childId].clear();
            }
#endif
        }
        std::sort(samplesForThisNode.begin(), samplesForThisNode.end());
        processNode(grg, samplesForThisNode, nodeId);
        return true;
    }

    void clearSampleSets() { m_sampleLists.clear(); }

private:
#ifdef CLEANUP_SAMPLE_SETS
    std::vector<NodeIDSizeT> m_refCounts;
#endif
    std::vector<NodeIDList> m_sampleLists;
};

/**
 * Get the frontier of a list of starting nodes (seeds). The frontier is the first node on each
 * path that is reached by all seeds. Some paths may not have any nodes that are reached by all
 * seeds, in which case no nodes from that path will be included in the result.
 *
 * This is a sparse implementation, as the frontier is most valuable when there is a lot
 * of sharing.
 */
class FrontierVisitor : public GRGVisitor {
public:
    explicit FrontierVisitor(const grgl::NodeIDList& seedList)
        : m_numSeeds(seedList.size()) {
        for (const NodeID seed : seedList) {
            m_reached.emplace(seed, 1);
        }
    }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        // Does not support DFS, just topo.
        assert(dfsPass != grgl::DfsPass::DFS_PASS_THERE);

        const NodeIDSizeT reached = m_reached.at(nodeId);
        // All seeds have reached the node, add to frontier and terminate traversal
        // along this particular path.
        if (reached == m_numSeeds) {
            m_frontier.emplace_back(nodeId);
            return false;
        }
        const auto& successors = (direction == DIRECTION_UP) ? grg->getUpEdges(nodeId) : grg->getDownEdges(nodeId);
        for (const auto& succId : successors) {
            auto insertIt = m_reached.emplace(succId, 0);
            insertIt.first->second += reached;
        }
        return true;
    }

    NodeIDList m_frontier;

private:
    std::unordered_map<NodeID, NodeIDSizeT> m_reached;
    const NodeIDSizeT m_numSeeds;
};

} // namespace grgl

#endif /* GRG_COMMON_VISITORS_H */
