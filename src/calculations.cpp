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
#include "calculations.h"

#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>

#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/visitor.h"
#include "util.h"

using namespace grgl;

/**
 * Visitor that computes allele frequency. Can either be used via downward depth-first search
 * (start at mutation ndoes) or upward via topological order (start at sample nodes).
 */
class AlleleFreqVisitor : public grgl::GRGVisitor {
public:
    AlleleFreqVisitor() = default;

    bool visit(const grgl::GRGPtr& grg,
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

void emitAlleleCounts(grgl::GRGPtr& grg,
                      std::ostream& outStream,
                      std::pair<uint32_t, uint32_t> bpRange,
                      const grgl::NodeIDList& onlySamples) {
    static constexpr char SEP = '\t';
    AlleleFreqVisitor visitorForDfs;
    if (bpRange.first == bpRange.second && onlySamples.empty()) {
        fastCompleteDFS(grg, visitorForDfs);
    } else if (!onlySamples.empty()) {
        if (bpRange.first != bpRange.second) {
            throw ApiMisuseFailure("--region and --sample-subset cannot be combined");
        }
        grg->visitTopo(visitorForDfs, grgl::TraversalDirection::DIRECTION_UP, onlySamples);
    } else {
        grgl::NodeIDList seeds;
        for (const auto& tuple : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
            const MutationId& mutId = std::get<1>(tuple);
            const grgl::Mutation& mut = grg->getMutationById(mutId);
            if (mut.getPosition() >= bpRange.first && mut.getPosition() < bpRange.second) {
                const NodeID& node = std::get<0>(tuple);
                if (node != INVALID_NODE_ID) {
                    seeds.push_back(node);
                }
                const NodeID& missingnessNode = std::get<2>(tuple);
                if (missingnessNode != INVALID_NODE_ID) {
                    seeds.push_back(missingnessNode);
                }
            }
        }
        if (seeds.empty()) {
            std::cout << "No variant in range" << std::endl;
            return;
        }
        grg->visitDfs(visitorForDfs, grgl::TraversalDirection::DIRECTION_DOWN, seeds);
    }
    const auto& mutsNodesMisses = grg->getMutationsToNodeOrdered<GRG::MutNodeMiss>();
    outStream << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP << "MISS COUNT" << SEP
              << "TOTAL" << std::endl;
    size_t i = 0;
    while (i < mutsNodesMisses.size()) {
        size_t samplesMissing = 0;
        size_t samplesWithMut = 0;
        const grgl::MutationId& mutId = std::get<0>(mutsNodesMisses[i]);
        const grgl::Mutation* mut = &grg->getMutationById(mutId);
        const grgl::BpPosition position = mut->getPosition();
        const std::string& ref = mut->getRefAllele();
        const std::string& alt = mut->getAllele();
        if (bpRange.first != bpRange.second && (bpRange.first > position || bpRange.second <= position)) {
            i++;
            continue;
        }
        // Accumulate values for all Mutations capturing the same (position, alleles) values.
        do {
            const grgl::NodeID& nodeId = std::get<1>(mutsNodesMisses[i]);
            if (nodeId != INVALID_NODE_ID) {
                samplesWithMut += visitorForDfs.m_samplesBeneath[nodeId];
            }
            const grgl::NodeID& missingnessNode = std::get<2>(mutsNodesMisses[i]);
            if (missingnessNode != INVALID_NODE_ID) {
                samplesMissing += visitorForDfs.m_samplesBeneath[missingnessNode];
            }
            i++;
            if (i < mutsNodesMisses.size()) {
                const grgl::MutationId& nextMutId = std::get<0>(mutsNodesMisses[i]);
                mut = &grg->getMutationById(nextMutId);
            }
        } while (i < mutsNodesMisses.size() && position == mut->getPosition() && ref == mut->getRefAllele() &&
                 alt == mut->getAllele());
        outStream << position << SEP << ref << SEP << alt << SEP << samplesWithMut << SEP << samplesMissing << SEP
                  << grg->numSamples() << std::endl;
    }
}

/**
 * Visitor that computes the number of heterozygous and homozygous individuals below nodes.
 */
class ZygosityInfoVisitor : public grgl::GRGVisitor {
public:
    ZygosityInfoVisitor() = default;

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (m_samplesBeneath.empty()) {
            m_samplesBeneath.resize(grg->numNodes());
            m_homozygousBeneath.resize(grg->numNodes());
        }
        if (dfsPass == grgl::DfsPass::DFS_PASS_BACK_AGAIN) {
            // Depth-first search must go down
            release_assert(direction == grgl::TraversalDirection::DIRECTION_DOWN);

            grgl::NodeIDSizeT homozygBeneath = 0;
            grgl::NodeIDSizeT samplesBeneath = 0;
            if (grg->isSample(nodeId)) {
                samplesBeneath++;
            } else {
                // Start with the number of individuals that coalesce at exactly this node.
                homozygBeneath = grg->getNumIndividualCoals(nodeId);
                release_assert(homozygBeneath != COAL_COUNT_NOT_SET && homozygBeneath <= grg->numIndividuals());

                for (const auto& child : grg->getDownEdges(nodeId)) {
                    samplesBeneath += m_samplesBeneath[child];
                    homozygBeneath += m_homozygousBeneath[child];
                }
            }
            m_samplesBeneath[nodeId] = samplesBeneath;
            m_homozygousBeneath[nodeId] = homozygBeneath;
            release_assert(samplesBeneath <= grg->numSamples());
        }
        release_assert(dfsPass != grgl::DfsPass::DFS_PASS_NONE);
        return true;
    }

    std::vector<grgl::NodeIDSizeT> m_homozygousBeneath;
    std::vector<grgl::NodeIDSizeT> m_samplesBeneath;
};

void emitZygosityInfo(grgl::GRGPtr& grg,
                      std::ostream& outStream,
                      std::pair<uint32_t, uint32_t> bpRange,
                      const grgl::NodeIDList& onlySamples) {
    static constexpr char SEP = '\t';
    if (bpRange.first != bpRange.second || !onlySamples.empty()) {
        std::cerr << "TODO: support range/sample subsets for zygosity info calculation." << std::endl;
        return;
    }
    if (grg->getPloidy() != 2) {
        std::cerr << "Calculating zygosity information does not work for non-diploid data." << std::endl;
        return;
    }
    ZygosityInfoVisitor visitorForDfs;
    fastCompleteDFS(grg, visitorForDfs);
    std::map<grgl::Mutation, std::pair<grgl::NodeIDSizeT, grgl::NodeIDSizeT>> counts;
    for (const auto& tuple : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
        const grgl::Mutation& mut = grg->getMutationById(std::get<1>(tuple));
        counts.insert({mut, {0, 0}});
        const NodeID& node = std::get<0>(tuple);
        if (node != INVALID_NODE_ID) {
            auto& countPair = counts.at(mut);
            countPair.first += visitorForDfs.m_samplesBeneath[node];
            countPair.second += visitorForDfs.m_homozygousBeneath[node];
        }
    }
    const size_t totalIndividuals = grg->numIndividuals();
    const auto& mutsNodesMisses = grg->getMutationsToNodeOrdered<GRG::MutNodeMiss>();
    outStream << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "AA" << SEP << "Aa" << SEP << "aa" << SEP
              << "MISS" << std::endl;
    size_t i = 0;
    while (i < mutsNodesMisses.size()) {
        size_t samplesMissing = 0;
        size_t samplesWithMut = 0;
        size_t homozygWithMut = 0;
        const grgl::MutationId mutId = std::get<0>(mutsNodesMisses[i]);
        const grgl::Mutation* mut = &grg->getMutationById(mutId);
        const grgl::BpPosition position = mut->getPosition();
        const std::string& ref = mut->getRefAllele();
        const std::string& alt = mut->getAllele();
        // Accumulate values for all Mutations capturing the same (position, alleles) values.
        do {
            const grgl::NodeID& nodeId = std::get<1>(mutsNodesMisses[i]);
            if (nodeId != INVALID_NODE_ID) {
                samplesWithMut += visitorForDfs.m_samplesBeneath[nodeId];
                homozygWithMut += visitorForDfs.m_homozygousBeneath[nodeId];
            }
            const grgl::NodeID& missingnessNode = std::get<2>(mutsNodesMisses[i]);
            if (missingnessNode != INVALID_NODE_ID) {
                samplesMissing += visitorForDfs.m_samplesBeneath[missingnessNode];
            }
            i++;
            if (i < mutsNodesMisses.size()) {
                const grgl::MutationId& nextMutId = std::get<0>(mutsNodesMisses[i]);
                mut = &grg->getMutationById(nextMutId);
            }
        } while (i < mutsNodesMisses.size() && position == mut->getPosition() && ref == mut->getRefAllele() &&
                 alt == mut->getAllele());

        const size_t heteroWithMut = samplesWithMut - (2 * homozygWithMut);
        const size_t neitherWithMut = totalIndividuals - (homozygWithMut + heteroWithMut);
        assert(neitherWithMut + heteroWithMut + homozygWithMut == totalIndividuals);
        outStream << position << SEP << ref << SEP << alt << SEP << neitherWithMut << SEP << heteroWithMut << SEP
                  << homozygWithMut << SEP << samplesMissing << std::endl;
    }
}
