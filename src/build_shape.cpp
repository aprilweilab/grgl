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
#include "build_shape.h"

#include <array>
#include <chrono>
#include <deque>
#include <limits>
#include <unordered_map>

#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "hap_index.h"
#include "similarity/bf_hash.h"
#include "similarity/mmh3.h"
#include "util.h"

namespace grgl {

static NodeID
createParentForNodes(const MutableGRGPtr& grg, NodeToHapVect& nodeHashes, NodeIDList& cluster, NodeIDSet& covered) {
    std::list<HaplotypeVector> similarHashes;
    // Now create a node to represent this new cluster, and connect its children.
    const auto newNodeId = grg->makeNode();
    for (const NodeID similarId : cluster) {
        covered.insert(similarId);
        grg->connect(newNodeId, similarId);
        const auto& similarHash = nodeHashes.at(similarId);
        similarHashes.push_back(similarHash);
    }
    auto newGenotypeHash = bitwiseIntersect(similarHashes);
    assert(nodeHashes.size() == newNodeId);
    nodeHashes.emplace_back(std::move(newGenotypeHash));
    return newNodeId;
}

static NodeIDList addGrgShapeFromHashing(const MutableGRGPtr& grg,
                                         NodeToHapVect& nodeHashes,
                                         const NodeIDList& initialNodes,
                                         const bool binaryTrees) {
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
            NodeIDList similar = hashIndex.getMostSimilarNodes(nodeId, !binaryTrees);
            if (!similar.empty()) {
                similar.push_back(nodeId);
                NodeID newNodeId = createParentForNodes(grg, nodeHashes, similar, covered);
                nextLevelNodes.push_back(newNodeId);
                hashIndex.add(newNodeId);
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
        const size_t mutRefs = mutAndSamples.second.size();
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

bool genotypeHashIndex(MutationIterator& mutIterator,
                       grgl::NodeToHapVect& hashIndex,
                       const size_t bitsPerMutation,
                       const double dropBelowThreshold) {
    constexpr size_t MAX_ALLELE_SIZE_FOR_HASH = 4;
    struct PosAndAllele {
        size_t position;
        char allele[MAX_ALLELE_SIZE_FOR_HASH];
    };

    size_t ploidy = 0;
    size_t numIndividuals = 0;
    bool isPhased = false;
    mutIterator.getMetadata(ploidy, numIndividuals, isPhased);
    const size_t numSamples = ploidy * numIndividuals;

    size_t dropBelowCount = 0;
    if (dropBelowThreshold < 1.0) {
        dropBelowCount = numSamples * dropBelowThreshold;
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
    size_t numVariants = 0;
    while (mutIterator.next(mutAndSamples, _ignore)) {
        if (mutAndSamples.second.size() < dropBelowCount) {
            dropped++;
            continue;
        }
        const size_t position = static_cast<size_t>(mutAndSamples.first.getPosition());
        PosAndAllele hashInput = {position, {}};
        strncpy(&hashInput.allele[0], mutAndSamples.first.getAllele().c_str(), MAX_ALLELE_SIZE_FOR_HASH);
        for (auto sampleId : mutAndSamples.second) {
            bloomFilters.at(sampleId).addItem(hashInput);
        }
        numVariants++;
    }
    std::cout << "Creating hash index with " << bloomFilters.size() << " entries for " << numVariants << " variants"
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
    return true;
}

MutableGRGPtr createEmptyGRGFromSamples(const std::string& sampleFile,
                                        FloatRange& genomeRange,
                                        size_t bitsPerMutation,
                                        const bool useBinaryMuts,
                                        const bool emitMissingData,
                                        const bool flipRefMajor,
                                        const double dropBelowThreshold) {
    MutableGRGPtr result;
    NodeToHapVect hashIndex;
    std::cout << "Building genotype hash index..." << std::endl;
    std::shared_ptr<grgl::MutationIterator> mutationIterator =
        makeMutationIterator(sampleFile, genomeRange, useBinaryMuts, emitMissingData, flipRefMajor);
    auto operationStartTime = std::chrono::high_resolution_clock::now();
    bool indexedOk = genotypeHashIndex(*mutationIterator, hashIndex, bitsPerMutation, dropBelowThreshold);
    if (!indexedOk) {
        std::cerr << "Failed constructing hash index\n";
        return result;
    }
    std::cout << "** Hashing input took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                       operationStartTime)
                     .count()
              << " ms\n";

    std::cout << "Done" << std::endl;
    result = std::make_shared<MutableGRG>(hashIndex.size());
    std::cout << "Adding GRG shape from genotype hashes..." << std::endl;
    const bool buildBinaryTrees = true;
    addGrgShapeFromHashing(result, hashIndex, result->getSampleNodes(), buildBinaryTrees);
    std::cout << "Done" << std::endl;

    return result;
}

} // namespace grgl
