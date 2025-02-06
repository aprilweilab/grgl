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

} // namespace grgl

#endif /* GRG_COMMON_VISITORS_H */
