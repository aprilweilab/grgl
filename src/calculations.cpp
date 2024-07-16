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
 * should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "calculations.h"

#include <iostream>
#include <unordered_map>
#include <vector>

#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/visitor.h"
#include "util.h"

/**
 * Visitor that computes allele frequency. Can either be used via downward depth-first search
 * (start at mutation ndoes) or upward via topological order (start at sample nodes).
 */
class AlleleFreqVisitor : public grgl::GRGVisitor {
public:
    AlleleFreqVisitor() = default;

    bool visit(const grgl::ConstGRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (m_samplesBeneath.empty()) {
            m_samplesBeneath.resize(grg->numNodes());
        }
        if (dfsPass == grgl::DfsPass::DFS_PASS_BACK_AGAIN) {
            // Depth-first search must go down
            release_assert(direction == grgl::TraversalDirection::DIRECTION_DOWN);
            grgl::NodeIDSizeT samplesBeneath = 0;
            if (grg->isSample(nodeId)) {
                samplesBeneath++;
            }
            for (const auto& child : grg->getDownEdges(nodeId)) {
                samplesBeneath += m_samplesBeneath[child];
            }
            m_samplesBeneath[nodeId] = samplesBeneath;
            release_assert(samplesBeneath <= grg->numSamples());
        } else if (dfsPass == grgl::DfsPass::DFS_PASS_NONE) {
            // Topological order must go up.
            release_assert(direction == grgl::TraversalDirection::DIRECTION_UP);
            for (const auto& parent : grg->getUpEdges(nodeId)) {
                if (grg->isSample(nodeId)) {
                    m_samplesBeneath[parent]++;
                } else {
                    m_samplesBeneath[parent] += m_samplesBeneath[nodeId];
                }
            }
        }
        return true;
    }

    std::vector<grgl::NodeIDSizeT> m_samplesBeneath;
};

void emitAlleleFrequency(grgl::GRGPtr& grg,
                         std::ostream& outStream,
                         std::pair<uint32_t, uint32_t> bpRange,
                         const grgl::NodeIDList& onlySamples) {
    static constexpr char SEP = '\t';
    AlleleFreqVisitor visitorForDfs;
    if (bpRange.first == bpRange.second && onlySamples.empty()) {
        fastCompleteDFS(grg, visitorForDfs);
    } else if (!onlySamples.empty()) {
        grg->visitTopo(visitorForDfs, grgl::TraversalDirection::DIRECTION_UP, onlySamples);
    } else {
        grgl::NodeIDList seeds;
        for (const auto& pair : grg->getNodeMutationPairs()) {
            const grgl::Mutation& mut = grg->getMutationById(pair.second);
            if (mut.getPosition() >= bpRange.first && mut.getPosition() < bpRange.second &&
                pair.first != INVALID_NODE_ID) {
                seeds.push_back(pair.first);
            }
        }
        if (seeds.empty()) {
            std::cout << "No variant in range" << std::endl;
            return;
        }
        grg->visitDfs(visitorForDfs, grgl::TraversalDirection::DIRECTION_DOWN, seeds);
    }
    std::map<grgl::Mutation, grgl::NodeIDSizeT> counts;
    for (const auto& pair : grg->getNodeMutationPairs()) {
        const grgl::Mutation& mut = grg->getMutationById(pair.second);
        if (bpRange.first != bpRange.second &&
            (bpRange.first > mut.getPosition() || bpRange.second <= mut.getPosition())) {
            continue;
        }
        counts.emplace(mut, 0);
        if (pair.first != INVALID_NODE_ID) {
            counts.at(mut) += visitorForDfs.m_samplesBeneath[pair.first];
        }
    }
    outStream << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP << "TOTAL" << std::endl;
    for (const auto& mutAndCount : counts) {
        const auto& mut = mutAndCount.first;
        outStream << mut.getPosition() << SEP << mut.getRefAllele() << SEP << mut.getAllele() << SEP
                  << mutAndCount.second << SEP << grg->numSamples() << std::endl;
    }
}
