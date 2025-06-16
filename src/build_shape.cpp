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
#include "build_shape.h"

#include <array>
#include <chrono>
#include <deque>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "grgl/mutation.h"
#include "hap_index.h"
#include "similarity/bf_hash.h"
#include "util.h"

namespace grgl {

static NodeIDList addGrgShapeFromHashing(const MutableGRGPtr& grg,
                                         NodeToHapVect& nodeHashes,
                                         const NodeIDList& initialNodes,
                                         const size_t tripletLevels) {
    // Lambda that does hamming distance between nodes. Instead of passing around the vectors contains the
    // bloom filters / hashes, we use a BK-tree that has a distance callback between elements.
    auto compareNodeIds = [&](const NodeID& node1, const NodeID& node2) {
        return bitwiseHamming(nodeHashes.at(node1), nodeHashes.at(node2));
    };
    HaplotypeIndex hashIndex(compareNodeIds);

    NodeIDSet covered;
    NodeIDList levelNodes;

    auto operationStartTime = std::chrono::high_resolution_clock::now();
    for (NodeID nodeId : initialNodes) {
        const auto& genotypeHash = nodeHashes.at(nodeId);
        hashIndex.add(nodeId);
        levelNodes.push_back(nodeId);
    }
    std::cout << "** Building index took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       operationStartTime)
                     .count()
              << " ms\n";
    operationStartTime = std::chrono::high_resolution_clock::now();

    size_t level = 0;
    bool createdNodes = true;
    while (createdNodes) {
        std::cout << "Pass " << level << " -- " << levelNodes.size() << std::endl;
        NodeIDList nextLevelNodes;

        createdNodes = false;
        while (!levelNodes.empty()) {
            auto nodeId = levelNodes.back();
            levelNodes.pop_back();

            // Skip nodes we've already covered. Note: we could relax this if we wanted to create
            // more complex graphs.
            auto coveredIt = covered.find(nodeId);
            if (coveredIt != covered.end()) {
                continue;
            }
            covered.insert(nodeId);

            // Find all similar nodes that create a "cluster" with this node, based on mutations

            // NOTE: this allows using nodes from the next level (recently created) as well as any previous
            // level, IF the node has not already been used as a child. This is CRITICAL to the performance
            // of this algorithm; specifically using nodes that were just recently created is important (think
            // of it as: unbalanced trees are important)
            NodeIDList similar = hashIndex.getMostSimilarNodes(nodeId, false);
            if (!similar.empty()) {
                // Pick the first element from similar as the representative. They should all be identical.
                NodeID first = similar.front();
                size_t isize = 0;
                HaplotypeVector vect1 = nodeHashes.at(nodeId);
                HaplotypeVector vect2 = nodeHashes.at(first);
                const size_t vect1Size = countBits(vect1);
                const size_t vect2Size = countBits(vect2);
                HaplotypeVector intersect = bitwiseIntersect(vect1, vect2, isize);
                // We create up to three nodes: (A & B), (A - B), (B - A) for any of them that are non-empty.
                if (isize > 0) {
                    NodeID newNodeId = grg->makeNode();
                    nodeHashes.emplace_back(intersect);
                    nextLevelNodes.push_back(newNodeId);
                    hashIndex.add(newNodeId);

                    grg->connect(newNodeId, nodeId);
                    for (const NodeID similarId : similar) {
                        grg->connect(newNodeId, similarId);
                    }
                }
                if ((level < tripletLevels) && (isize < vect1Size)) {
                    const size_t subSize = bitwiseSubtract(vect1, intersect);
                    if (subSize > 0) {
                        NodeID newNodeId = grg->makeNode();
                        nodeHashes.emplace_back(std::move(vect1));
                        nextLevelNodes.push_back(newNodeId);
                        hashIndex.add(newNodeId);
                        grg->connect(newNodeId, nodeId);
                    }
                }
                if ((level < tripletLevels) && (isize < vect2Size)) {
                    const size_t subSize = bitwiseSubtract(vect2, intersect);
                    if (subSize > 0) {
                        NodeID newNodeId = grg->makeNode();
                        nodeHashes.emplace_back(std::move(vect2));
                        nextLevelNodes.push_back(newNodeId);
                        hashIndex.add(newNodeId);
                        for (const NodeID similarId : similar) {
                            grg->connect(newNodeId, similarId);
                        }
                    }
                }

                for (const NodeID similarId : similar) {
                    covered.insert(similarId);
                }
                covered.insert(nodeId);

                createdNodes = true;
            } else {
                nextLevelNodes.push_back(nodeId);
                // getMostSimilarNodes "removes" nodeId, so we have to add it back
                hashIndex.add(nodeId);
            }
        }
        level++;
        levelNodes = std::move(nextLevelNodes);
        hashIndex.emitStats();
    }
    std::cout << "** Constructing tree took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       operationStartTime)
                     .count()
              << " ms\n";
    hashIndex.emitStats();
    return levelNodes;
}

static void getMutStats(MutationIterator& iterator,
                        size_t& avgMutsPerSample,
                        size_t& mutCutoff25Percentile,
                        size_t& numSamples,
                        const size_t dropBelowCount) {
    size_t totalMutRefs = 0;
    std::vector<size_t> mutSampleCounts;
    iterator.reset();
    MutationAndSamples mutAndSamples;
    while (iterator.next(mutAndSamples, numSamples)) {
        const size_t mutRefs = mutAndSamples.samples.size();
        if (mutRefs < dropBelowCount) {
            continue;
        }
        totalMutRefs += mutRefs;
        mutSampleCounts.push_back(mutRefs);
    }
    if (!mutSampleCounts.empty()) {
        avgMutsPerSample = totalMutRefs / numSamples;
        std::sort(mutSampleCounts.begin(), mutSampleCounts.end());
        const size_t idx25Percent = mutSampleCounts.size() / 4;
        release_assert(idx25Percent < mutSampleCounts.size());
        mutCutoff25Percentile = mutSampleCounts[idx25Percent];
    } else {
        avgMutsPerSample = 0;
        mutCutoff25Percentile = 0;
    }
}

uint16_t genotypeHashIndex(MutationIterator& mutIterator,
                           grgl::NodeToHapVect& hashIndex,
                           const size_t bitsPerMutation,
                           const double dropBelowThreshold) {
    size_t ploidy = 0;
    size_t numIndividuals = 0;
    bool isPhased = false;
    mutIterator.getMetadata(ploidy, numIndividuals, isPhased);
    const size_t numSamples = ploidy * numIndividuals;

    size_t dropBelowCount = 0;
    if (dropBelowThreshold < 1.0) {
        dropBelowCount = (size_t)((double)numSamples * dropBelowThreshold);
    } else {
        dropBelowCount = (size_t)dropBelowThreshold;
    }

    size_t _ignore = 0;
    size_t avgMutsPerSample = 0;
    size_t mutCutoff25Percentile = 0;
    getMutStats(mutIterator, avgMutsPerSample, mutCutoff25Percentile, _ignore, dropBelowCount);
    const size_t bitsPerElement = (8 * sizeof(HapVectorT));
    // TODO: there can be some really sparse scenarios, mostly with filtered data (like only the homozygous
    // mutations), and this helps. It would be better to compute this based on total variant count.
    constexpr size_t MIN_ELEMENTS = 10;
    const size_t vectorElemCount =
        std::max<size_t>(MIN_ELEMENTS, (bitsPerMutation * avgMutsPerSample) / bitsPerElement);
    const size_t vectorBitSize = vectorElemCount * bitsPerElement;

    std::cout << "Using genotype hashes of length " << vectorElemCount << std::endl;

    std::vector<BFHash> bloomFilters(numSamples, BFHash(vectorBitSize));

    MutationAndSamples mutAndSamples = {Mutation(0.0, ""), NodeIDList()};
    mutIterator.reset();
    size_t dropped = 0;
    size_t variantCount = 0;
    while (mutIterator.next(mutAndSamples, _ignore)) {
        if (mutAndSamples.samples.size() < dropBelowCount) {
            dropped++;
            continue;
        }
        // std::hash() of an integer is just the integer, and we get better results with a little
        // bit more random-like behavior.
        const size_t hashInput = hash_combine(std::hash<size_t>{}(variantCount), 42);
        for (auto sampleId : mutAndSamples.samples) {
            bloomFilters.at(sampleId).addHash(hashInput);
        }
        variantCount++;
    }
    std::cout << "Creating hash index with " << bloomFilters.size() << " entries for " << variantCount << " variants"
              << std::endl;
    std::cout << "Dropped " << dropped << " variants" << std::endl;
    hashIndex = grgl::NodeToHapVect(bloomFilters.size());
    for (size_t sampleId = 0; sampleId < bloomFilters.size(); sampleId++) {
        hashIndex[sampleId] = std::move(bloomFilters[sampleId].stealVector());
    }
    if (mutIterator.numFlippedAlleles() > 0) {
        std::cout << "Flipped " << mutIterator.numFlippedAlleles() << " reference alleles (to the major allele)"
                  << std::endl;
    }
    return static_cast<uint16_t>(ploidy);
}

MutableGRGPtr createEmptyGRGFromSamples(const std::string& sampleFile,
                                        FloatRange& genomeRange,
                                        size_t bitsPerMutation,
                                        GrgBuildFlags buildFlags,
                                        const double dropBelowThreshold,
                                        const std::map<std::string, std::string>& indivIdToPop,
                                        const size_t tripletLevels) {
    MutableGRGPtr result;
    NodeToHapVect hashIndex;

#define GRGBS_LOG_OUTPUT(msg)                                                                                          \
    do {                                                                                                               \
        if (static_cast<bool>(buildFlags & GBF_VERBOSE_OUTPUT)) {                                                      \
            std::cerr << msg;                                                                                          \
        }                                                                                                              \
    } while (0)

    GRGBS_LOG_OUTPUT("Building genotype hash index..." << std::endl);
    const bool useBinaryMuts = static_cast<bool>(buildFlags & GBF_USE_BINARY_MUTS);
    const bool emitMissingData = static_cast<bool>(buildFlags & GBF_EMIT_MISSING_DATA);
    const bool flipRefMajor = static_cast<bool>(buildFlags & GBF_FLIP_REF_MAJOR);
    std::shared_ptr<grgl::MutationIterator> mutationIterator =
        makeMutationIterator(sampleFile, genomeRange, useBinaryMuts, emitMissingData, flipRefMajor);
    auto operationStartTime = std::chrono::high_resolution_clock::now();
    uint16_t ploidy = genotypeHashIndex(*mutationIterator, hashIndex, bitsPerMutation, dropBelowThreshold);
    GRGBS_LOG_OUTPUT("** Hashing input took " << std::chrono::duration_cast<std::chrono::milliseconds>(
                                                     std::chrono::high_resolution_clock::now() - operationStartTime)
                                                     .count()
                                              << " ms" << std::endl);

    GRGBS_LOG_OUTPUT("Done" << std::endl);
    result = std::make_shared<MutableGRG>(hashIndex.size(), ploidy);

    std::vector<std::string> indivIds;
    if (!indivIdToPop.empty()) {
        size_t ploidy = 0;
        size_t numIndividuals = 0;
        bool isPhased = false;
        mutationIterator->getMetadata(ploidy, numIndividuals, isPhased);
        indivIds = mutationIterator->getIndividualIds();
        if (!indivIds.empty()) {
            release_assert(indivIds.size() == numIndividuals);
            std::map<std::string, size_t> popDescriptionMap;
            for (NodeID individual = 0; individual < indivIds.size(); individual++) {
                const auto& stringId = indivIds[individual];
                const auto& findIt = indivIdToPop.find(stringId);
                if (findIt == indivIdToPop.end()) {
                    std::stringstream ssErr;
                    ssErr << "Could not find population mapping for individual " << stringId;
                    throw std::runtime_error(ssErr.str());
                }
                const auto& popDescription = findIt->second;
                const size_t nextPopId = popDescriptionMap.size();
                const auto& findPopIt = popDescriptionMap.emplace(popDescription, nextPopId);
                const auto popId = findPopIt.first->second;
                if (findPopIt.second) {
                    release_assert(popId == nextPopId);
                    result->addPopulation(popDescription);
                } else {
                    release_assert(popId != nextPopId);
                }
                for (NodeID offset = 0; offset < ploidy; offset++) {
                    const NodeID sampleId = (individual * ploidy) + offset;
                    release_assert(sampleId < result->numSamples());
                    result->setPopulationId(sampleId, popId);
                }
            }
        } else {
            throw ApiMisuseFailure("No individual IDs present in input; required for population mapping");
        }
    }
    GRGBS_LOG_OUTPUT("Adding GRG shape from genotype hashes..." << std::endl);
    addGrgShapeFromHashing(result, hashIndex, result->getSampleNodes(), tripletLevels);
    const auto actualRange = mutationIterator->getBpRange();
    result->setSpecifiedBPRange({(BpPosition)actualRange.start(), (BpPosition)actualRange.end()});
    GRGBS_LOG_OUTPUT("Done" << std::endl);

    // Add individual identifiers.
    if (!static_cast<bool>(buildFlags & GBF_NO_INDIVIDUAL_IDS)) {
        if (indivIds.empty()) {
            indivIds = mutationIterator->getIndividualIds();
        }
        if (indivIds.empty()) {
            GRGBS_LOG_OUTPUT("WARNING: No individual identifiers present in the input file" << std::endl);
        }
        for (const auto& ident : indivIds) {
            result->addIndividualId(ident);
        }
    }

#undef GRGBS_LOG_OUTPUT

    return result;
}

} // namespace grgl
