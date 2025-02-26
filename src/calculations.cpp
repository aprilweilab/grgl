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

#if GSL_ENABLED
extern "C" {
#include "cdf/gsl_cdf.h"
}
#endif

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

void emitAlleleFrequency(grgl::GRGPtr& grg,
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
    const auto& mutIdAndNodes = grg->getMutationsToNodeOrdered();
    outStream << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "ALT COUNT" << SEP << "TOTAL" << std::endl;
    size_t i = 0;
    while (i < mutIdAndNodes.size()) {
        size_t samplesWithMut = 0;
        const grgl::MutationId mutId = mutIdAndNodes[i].first;
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
            const grgl::NodeID& nodeId = mutIdAndNodes[i].second;
            if (nodeId != INVALID_NODE_ID) {
                samplesWithMut += visitorForDfs.m_samplesBeneath[nodeId];
            }
            i++;
            if (i < mutIdAndNodes.size()) {
                const grgl::MutationId& nextMutId = mutIdAndNodes[i].first;
                mut = &grg->getMutationById(nextMutId);
            }
        } while (i < mutIdAndNodes.size() && position == mut->getPosition() && ref == mut->getRefAllele() &&
                 alt == mut->getAllele());
        outStream << position << SEP << ref << SEP << alt << SEP << samplesWithMut << SEP << grg->numSamples()
                  << std::endl;
    }
}

void randomPhenotypeData(const size_t seed,
                         const double mean,
                         const double stddev,
                         std::vector<double>& phenVector,
                         double& totalPheno,
                         double& YY) {
    static std::mt19937 generator;
    generator.seed(seed);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<double> dist{mean, stddev};

    for (auto& value : phenVector) {
        value = dist(generator);
        totalPheno += value;
        YY += value * value;
    }
}

void loadPhenotypeData(const std::string& phenotypeTextFile,
                       std::vector<double>& phenVector,
                       double& totalPheno,
                       double& YY) {
    std::ifstream inStream(phenotypeTextFile);
    if (!inStream.good()) {
        std::stringstream ssErr;
        ssErr << "Could not read " << phenotypeTextFile;
        throw grgl::BadInputFileFailure(ssErr.str().c_str());
    }
    std::string myText;
    size_t i = 0;
    while (getline(inStream, myText) && i < phenVector.size()) {
        std::vector<std::string> v = split(myText, ' ');
        if (v.size() != 3) {
            std::stringstream ssErr;
            ssErr << "Each line in phenotype file must have three space-separated columns." << " Line " << i
                  << " failed this check.";
            throw grgl::BadInputFileFailure(ssErr.str().c_str());
        }
        phenVector[i] = std::stof(v[2]);
        YY += phenVector[i] * phenVector[i];
        totalPheno += phenVector[i];
        i += 1;
    }
    if (i < phenVector.size()) {
        std::cerr << "WARNING: Too few individual phenotype values in " << phenotypeTextFile
                  << ", remaining individuals will have 0.0 for phenotype value." << std::endl;
    }
}

void emitBeta(const grgl::GRGPtr& grg, const std::string& phenotype, std::ostream& outStream, const bool betaOnly) {
    const NodeIDSizeT num_nodes = grg->numNodes();
    const NodeIDSizeT num_samples = grg->numSamples() / 2;

    if (grg->getPloidy() != 2) {
        std::cerr << "GWAS only works for diploid data (your GRG has ploidy=" << grg->getPloidy() << ")" << std::endl;
        return;
    }

#if !GSL_ENABLED
    if (!betaOnly) {
        std::cerr << "WARNING: There will be no p-values calculated. Enable GSL (GNU Scientific Library) "
                     "in your GRGL build if you want p-values"
                  << std::endl;
    }
#endif

    double totalPheno = 0;
    double YY = 0;
    std::vector<double> phen(num_samples);
    if (phenotype == USE_RANDOM_PHENOTYPE) {
        randomPhenotypeData(0, 0, 10, phen, totalPheno, YY);
    } else {
        loadPhenotypeData(phenotype, phen, totalPheno, YY);
    }

    std::vector<NodeIDSizeT> frequencyMap(num_nodes);
    std::vector<double> nodeXYcount(num_nodes);
    std::vector<size_t> nodeXXcount(num_nodes);

    for (grgl::NodeID node = 0; node < num_nodes; node++) {
        if (grg->isSample(node)) {
            nodeXXcount[node] = 1;
            frequencyMap[node] = 1;
            nodeXYcount[node] = phen.at(node / 2);
        } else {
            NodeIDSizeT frequency = 0;
            size_t nodeXX = 0;
            double nodeXY = 0;
            for (const auto& child : grg->getDownEdges(node)) {
                nodeXX += nodeXXcount[child];
                frequency += frequencyMap[child];
                nodeXY += nodeXYcount[child];
            }
            frequencyMap[node] = frequency;
            nodeXXcount[node] = nodeXX + 2 * grg->getNumIndividualCoals(node);
            nodeXYcount[node] = nodeXY;
        }
    }

    const auto& mutIdAndNodes = grg->getMutationsToNodeOrdered();
    const double dof = num_samples - 2;
    static constexpr char SEP = '\t';
    if (betaOnly) {
        outStream << "BP" << SEP << "COUNT" << SEP << "BETA" << std::endl;
    } else {
        outStream << "BP" << SEP << "COUNT" << SEP << "BETA" << SEP << "B0" << SEP << "SE" << SEP << "R2" << SEP << "T"
                  << SEP << "P" << std::endl;
    }

    size_t i = 0;
    while (i < mutIdAndNodes.size()) {
        double freq = 0;
        double countXY = 0;
        double countXX = 0;
        const grgl::MutationId mutId = mutIdAndNodes[i].first;
        // Accumulate values for all nodes associated with the Mutation. They are guaranteed to be consecutive by
        // the getMutationsToNodeOrder() method.
        do {
            const grgl::NodeID& nodeId = mutIdAndNodes[i].second;
            if (nodeId != INVALID_NODE_ID) {
                countXY += nodeXYcount[nodeId];
                countXX += (double)nodeXXcount[nodeId];
                freq += (double)frequencyMap[nodeId];
            }
            i++;
        } while (i < mutIdAndNodes.size() && mutIdAndNodes[i].first == mutId);
        const grgl::Mutation& mut = grg->getMutationById(mutId);

        if (freq == 0) {
            if (betaOnly) {
                outStream << mut.getPosition() << SEP << 0 << SEP << 0 << std::endl;
            } else {
                outStream << mut.getPosition() << SEP << 0 << SEP << 0 << SEP << 0 << SEP << 0 << SEP << 0 << SEP << 0
                          << SEP << 0 << std::endl;
            }
        } else {
            const double freqNormalized = freq / (double)num_samples;
            const double nodeXY = countXY - freqNormalized * (double)totalPheno;
            const double nodeXX = countXX - freq * freqNormalized;
            const double beta = nodeXY / nodeXX;
            if (betaOnly) {
                outStream << mut.getPosition() << SEP << freq << SEP << beta << std::endl;
            } else {
                const double b0 = (double)totalPheno / (double)num_samples - freqNormalized * beta;
                const double err = YY - 2 * b0 * totalPheno - 2 * beta * countXY + num_samples * b0 * b0 +
                                   2 * b0 * beta * freq + beta * beta * countXX;
                double se = std::sqrt(err / (num_samples - 2) / (nodeXX));
                const double t_val = beta / se;
                const double s_tot = YY - totalPheno * totalPheno / num_samples;
                const double r2 = 1 - (err / s_tot);
                double p_val = std::nan("1");

#if GSL_ENABLED
                double cdf = gsl_cdf_tdist_P(t_val, dof);
                if (t_val > 0) {
                    p_val = 2 * (1.0 - cdf);
                } else {
                    p_val = 2 * cdf;
                }
#endif

                outStream << mut.getPosition() << SEP << freq << SEP << beta << SEP << b0 << SEP << se << SEP
                          << std::scientific << r2 << std::fixed << std::setprecision(4) << SEP << t_val << SEP << p_val
                          << std::endl;
            }
        }
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
    for (const auto& pair : grg->getNodeMutationPairs()) {
        const grgl::Mutation& mut = grg->getMutationById(pair.second);
        counts.insert({mut, {0, 0}});
        if (pair.first != INVALID_NODE_ID) {
            auto& countPair = counts.at(mut);
            countPair.first += visitorForDfs.m_samplesBeneath[pair.first];
            countPair.second += visitorForDfs.m_homozygousBeneath[pair.first];
        }
    }
    const size_t totalIndividuals = grg->numIndividuals();
    const auto& mutIdAndNodes = grg->getMutationsToNodeOrdered();
    outStream << "POSITION" << SEP << "REF" << SEP << "ALT" << SEP << "AA" << SEP << "Aa" << SEP << "aa" << std::endl;
    size_t i = 0;
    while (i < mutIdAndNodes.size()) {
        size_t samplesWithMut = 0;
        size_t homozygWithMut = 0;
        const grgl::MutationId mutId = mutIdAndNodes[i].first;
        const grgl::Mutation* mut = &grg->getMutationById(mutId);
        const grgl::BpPosition position = mut->getPosition();
        const std::string& ref = mut->getRefAllele();
        const std::string& alt = mut->getAllele();
        // Accumulate values for all Mutations capturing the same (position, alleles) values.
        do {
            const grgl::NodeID& nodeId = mutIdAndNodes[i].second;
            if (nodeId != INVALID_NODE_ID) {
                samplesWithMut += visitorForDfs.m_samplesBeneath[nodeId];
                homozygWithMut += visitorForDfs.m_homozygousBeneath[nodeId];
            }
            i++;
            if (i < mutIdAndNodes.size()) {
                const grgl::MutationId& nextMutId = mutIdAndNodes[i].first;
                mut = &grg->getMutationById(nextMutId);
            }
        } while (i < mutIdAndNodes.size() && position == mut->getPosition() && ref == mut->getRefAllele() &&
                 alt == mut->getAllele());

        const size_t heteroWithMut = samplesWithMut - (2 * homozygWithMut);
        const size_t neitherWithMut = totalIndividuals - (homozygWithMut + heteroWithMut);
        assert(neitherWithMut + heteroWithMut + homozygWithMut == totalIndividuals);
        outStream << position << SEP << ref << SEP << alt << SEP << neitherWithMut << SEP << heteroWithMut << SEP
                  << homozygWithMut << std::endl;
    }
}
