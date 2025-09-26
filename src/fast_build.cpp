/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2025 April Wei
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
#include "fast_build.h"

#include <array>
#include <cassert>
#include <chrono>
#include <deque>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "grgl/mutation.h"
#include "grgl/windowing.h"
#include "util.h"

#include "hap_segment.h"

#ifdef DEBUG_SANITY_CHECKS
#include "grgl/visitor.h"

#define DEBUG_PRINT(msg) std::cout << msg
#else
#define DEBUG_PRINT(msg)
#endif

namespace grgl {

inline bool hasBSFlag(GrgBuildFlags flags, GrgBuildFlags flag) { return (bool)(flags & flag); }

void addExtraInfoToGRG(MutableGRGPtr& grg,
                       grgl::MutationIterator& mutationIterator,
                       GrgBuildFlags buildFlags,
                       const std::map<std::string, std::string>& indivIdToPop) {
#define GRGBS_LOG_OUTPUT(msg)                                                                                          \
    do {                                                                                                               \
        if (static_cast<bool>(buildFlags & GBF_VERBOSE_OUTPUT)) {                                                      \
            std::cerr << msg;                                                                                          \
        }                                                                                                              \
    } while (0)

    // Population mapping.
    std::vector<std::string> indivIds;
    if (!indivIdToPop.empty()) {
        size_t ploidy = 0;
        size_t numIndividuals = 0;
        bool isPhased = false;
        mutationIterator.getMetadata(ploidy, numIndividuals, isPhased);
        indivIds = mutationIterator.getIndividualIds();
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
                    grg->addPopulation(popDescription);
                } else {
                    release_assert(popId != nextPopId);
                }
                for (NodeID offset = 0; offset < ploidy; offset++) {
                    const NodeID sampleId = (individual * ploidy) + offset;
                    release_assert(sampleId < grg->numSamples());
                    grg->setPopulationId(sampleId, popId);
                }
            }
        } else {
            throw ApiMisuseFailure("No individual IDs present in input; required for population mapping");
        }
    }

    // Specified range.
    const auto actualRange = mutationIterator.getBpRange();
    grg->setSpecifiedBPRange({(BpPosition)actualRange.start(), (BpPosition)actualRange.end()});

    // Individual identifiers.
    if (!static_cast<bool>(buildFlags & GBF_NO_INDIVIDUAL_IDS)) {
        if (indivIds.empty()) {
            indivIds = mutationIterator.getIndividualIds();
        }
        if (indivIds.empty()) {
            GRGBS_LOG_OUTPUT("WARNING: No individual identifiers present in the input file" << std::endl);
        }
        for (const auto& ident : indivIds) {
            grg->addIndividualId(ident);
        }
    }

#undef GRGBS_LOG_OUTPUT
}

HaplotypeWindow createWindow(const size_t firstMutationIdx,
                             const size_t numMutations,
                             std::vector<HaplotypeVector>& window,
                             // This is the running vector of haplotype values per sample. So for example, sample 10
                             // would be sampleHapVects[10] = [hapIdx for window 0, hapIdx for window 1, ...]
                             std::vector<HapIdxList>& sampleHapVects,
                             size_t maxPrecomputed = 1024) {
    HaplotypeWindow result = {firstMutationIdx, firstMutationIdx + numMutations, {}, {}, {}, 0, {}};
    HapIdx nextHapIdx = 1;

    // The hap index 0 is always the empty vector (no mutations)
    result.hapMap.emplace(HaplotypeVector(window[0].size()), 0);

    for (NodeID sampleId = 0; sampleId < window.size(); sampleId++) {
        const auto inserted = result.hapMap.emplace(window[sampleId], nextHapIdx);
        // If we inserted it, it is new (unique so far).
        if (inserted.second) {
            nextHapIdx++;
        }
        const HapIdx samplesHapIdx = inserted.first->second;
        sampleHapVects[sampleId].emplace_back(samplesHapIdx);

        // Clear out the bitvector for this sample so we can reuse it for the next window.
        memset(window[sampleId].data(), 0, window[sampleId].size() * sizeof(HapVectorT));
    }

    const size_t numHaps = result.hapMap.size();
    release_assert(numHaps < MAX_FAST_HAPLOTYPES);
    result.distMatrix.resize(maxPrecomputed * maxPrecomputed, INVALID_HAP_DIST);

    // Move the haplotypes to a vector and order them by the haplotype index.
    for (const auto& hap : result.hapMap) {
        result.orderedHaps.emplace_back(hap.first);
    }
    release_assert(result.orderedHaps.size() == numHaps);
    auto hapsLessThan = [&](const HaplotypeVector& hap1, const HaplotypeVector& hap2) {
        return result.hapMap.at(hap1) < result.hapMap.at(hap2);
    };
    std::sort(result.orderedHaps.begin(), result.orderedHaps.end(), hapsLessThan);

    // Precompute the distances between all observed haplotypes.
    result.orderedHapCounts.resize(numHaps);
    for (size_t i = 0; i < numHaps; i++) {
        result.orderedHapCounts.at(i) = countBits(result.orderedHaps[i]);
    }
    result.precomputed = maxPrecomputed;

    DEBUG_PRINT("=== REGION ===\n");
    DEBUG_PRINT("Unique haplotypes: " << numHaps << std::endl);
    return std::move(result);
}

#ifdef DEBUG_SANITY_CHECKS
class UniqueSamplesCheck : public GRGVisitor {
public:
    UniqueSamplesCheck() = default;

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            release_assert(direction == TraversalDirection::DIRECTION_DOWN);

            NodeIDSet mySampleSet;
            NodeIDList mySamples;
            size_t failures = 0;
            for (const auto& childId : grg->getDownEdges(nodeId)) {
                if (grg->isSample(childId)) {
                    mySamples.push_back(childId);
                } else {
                    for (const auto s : m_sampleLists.at(childId)) {
                        if (mySampleSet.find(s) != mySampleSet.end()) {
                            failures++;
                            // std::cout << "FAIL! nodeId=" << nodeId << ", child=" << childId << "\n";
                        }
                        mySamples.push_back(s);
                        mySampleSet.emplace(s);
                    }
                }
            }
            m_sampleLists.emplace(nodeId, std::move(mySamples));
            if (failures > 0) {
                m_failures.emplace(nodeId, failures);
            }
        }
        return true;
    }
    std::unordered_map<grgl::NodeID, NodeIDSizeT> m_failures;

private:
    std::unordered_map<grgl::NodeID, NodeIDList> m_sampleLists;
};
#endif

std::vector<size_t> hapsToMutIndices(const HapIdxList& hapList, const std::vector<HaplotypeWindow>& windowInfo) {
    std::vector<size_t> result;
    release_assert(hapList.size() == windowInfo.size());
    for (size_t i = 0; i < hapList.size(); i++) {
        const HapIdx hapIdx = hapList[i];
        DEBUG_PRINT("CHECK HAP (" << i << "," << hapIdx << "), window_first=" << windowInfo[i].firstMutationIdx
                                  << ": ");
        const HaplotypeVector& mutBitVect = windowInfo[i].orderedHaps.at(hapIdx);
        const auto& mutList = getBitsAsList(mutBitVect);
        release_assert(hapIdx != 0 || mutList.empty());
        for (const size_t mutIdx : mutList) {
            DEBUG_PRINT((mutIdx + windowInfo[i].firstMutationIdx) << ", ");
            result.push_back(mutIdx + windowInfo[i].firstMutationIdx);
        }
        DEBUG_PRINT("\n");
    }
    return std::move(result);
}

// Just a context object to simplify passing this information around.
struct HapWindowContext {
    // All mutations within the region spanned by the hap segments.
    std::vector<Mutation> allMutations;
    // The mutations that are filtered out specifically to be directly mapped (low frequency).
    std::vector<MutationAndSamples> directMap;
    // The information for each of the hap segments (windows). There must be at least one.
    std::vector<HaplotypeWindow> windowInfo;
    // The actual haplotypes as bit vectors, this only contains _unique_ haplotypes.
    std::vector<HapIdxList> sampleHapVects;
};

/**
 * Use the MutationIterator to get all hap segments of the given length.
 */
void getHapSegments(MutationIterator& mutIterator,
                    const size_t hapLength,
                    const size_t numSamples,
                    const double dropBelowThreshold,
                    const size_t stopAfter,
                    const bool stopAfterIsMutCount,
                    HapWindowContext& hapContext) {
    release_assert(stopAfter > 0);
    size_t windex = 0;

    // If the threshold is <1.0, it is a frequency.
    size_t directMapBelowCt = 0;
    if (dropBelowThreshold < 1.0) {
        directMapBelowCt = (size_t)((double)numSamples * dropBelowThreshold);
    } else {
        directMapBelowCt = (size_t)dropBelowThreshold;
    }

    hapContext.sampleHapVects.resize(numSamples);
    std::vector<HaplotypeVector> window(numSamples, createHapVect(hapLength));
    size_t lastMutStart = 0;
    size_t windows = 0;

    // This is either mutations processed so far, or unique hap segments (sum) so far,
    // depending on stopAfterIsMutCount.
    size_t soFar = 0;

    bool done = false;
    size_t _ignore = 0;
    MutationAndSamples mutAndSamples;
    // The missing data Mutations always come first, so this tells us if the position of any _non-missing-data_
    // Mutation also has missing data.
    size_t lastPositionWithMissing = std::numeric_limits<size_t>::max();

    while (!done && mutIterator.next(mutAndSamples, _ignore)) {
        // Stopping when we reach a certain mutation threshold.
        if (stopAfterIsMutCount) {
            soFar += 1;
            if (stopAfter == soFar) {
                done = true;
            }
        }

        const size_t position = mutAndSamples.mutation.getPosition();
        if (mutAndSamples.mutation.isMissing()) {
            lastPositionWithMissing = position;
        }

        // We will only direct map mutations if there is no missing data at the site. Otherwise, tracking
        // the missingness NodeIDs across regular and direct mapping is a nightmare.
        if (mutAndSamples.samples.size() < directMapBelowCt && lastPositionWithMissing != position) {
            hapContext.directMap.emplace_back(std::move(mutAndSamples));
            continue;
        }

        hapContext.allMutations.push_back(std::move(mutAndSamples.mutation));
        for (size_t i = 0; i < mutAndSamples.samples.size(); i++) {
            setBit(window[mutAndSamples.samples[i]], windex);
        }
        windex++;

        if (windex == hapLength) {
            DEBUG_PRINT("WINDOW " << lastMutStart << "\n");
            hapContext.windowInfo.emplace_back(
                createWindow(lastMutStart, hapLength, window, hapContext.sampleHapVects));
            lastMutStart = hapContext.allMutations.size();
            windex = 0;
            windows++;

            // Stopping when we reach a certain unique haplotype sum.
            if (!stopAfterIsMutCount) {
                soFar += hapContext.windowInfo.back().orderedHaps.size();
                if (soFar >= stopAfter) {
                    done = true;
                }
            }
        }
    }
    if (windex > 0) {
        hapContext.windowInfo.push_back(createWindow(lastMutStart, windex, window, hapContext.sampleHapVects));
        windows++;
    }
}

NodeIDSizeT getCoalsForParent(MutableGRGPtr& grg,
                              std::unordered_map<NodeID, NodeIDList>& nodeToIndivs,
                              const NodeIDList& children,
                              std::unordered_set<NodeIDSizeT>& seenIndivs,
                              bool cleanup) {
    constexpr NodeIDSizeT ploidy = 2;
    NodeIDSizeT coalCount = 0;

    // Collect all "individuals below" each child and whenever we see one twice, count
    // it as a coalescence and remove it from the list of seen individuals.
    for (const NodeID child : children) {
        if (grg->isSample(child)) {
            const NodeIDSizeT indiv = child / ploidy;
            auto insertIt = seenIndivs.insert(indiv);
            if (!insertIt.second) {
                seenIndivs.erase(insertIt.first);
                coalCount++;
            }
        } else {
            auto findIt = nodeToIndivs.find(child);
            release_assert(findIt != nodeToIndivs.end());
            const NodeIDList& rightIndividuals = findIt->second;
            for (const NodeID indiv : findIt->second) {
                auto insertIt = seenIndivs.insert(indiv);
                if (!insertIt.second) {
                    seenIndivs.erase(insertIt.first);
                    coalCount++;
                }
            }
            if (cleanup) {
                nodeToIndivs.erase(findIt);
            }
        }
    }
    return coalCount;
}

// A GRG stores information about coalescences, in order to compute X^T*X efficiently for GWAS
// calculations. This information is pretty cheap to calculate, since we construct the graph
// bottom-up and only use each node once.
void propagateCoalInformation(MutableGRGPtr& grg,
                              std::unordered_map<NodeID, NodeIDList>& nodeToIndivs,
                              const NodeID intersectionParent,
                              const NodeID leftParent,
                              const NodeID rightParent,
                              const NodeID leftNode,
                              const NodeIDList& rightNodes) {
    // We only track coalescences for diploids right now. Haploids don't need it, and polyploids are
    // tricky because we have to track more than just a list below (we need a map).
    constexpr NodeID ploidy = 2;
    if (grg->getPloidy() != ploidy) {
        return;
    }
    release_assert(intersectionParent != INVALID_NODE_ID);
    release_assert(leftNode != INVALID_NODE_ID);

    // Get the RIGHT side individuals, from maybe multiple nodes. We count the coalescences
    // as we go. Everytime an individual is seen twice, it has coalesced (diploid), so we then
    // remove it since the parent nodes do not need to track the coalescence any further.
    std::unordered_set<NodeIDSizeT> nextIndivs;
    NodeIDList intersectIndividuals;
    NodeIDSizeT rightCoalCount =
        getCoalsForParent(grg, nodeToIndivs, rightNodes, nextIndivs, /*cleanup=*/(rightParent != INVALID_NODE_ID));

    // At this point:
    //  nodeToIndivs has the individuals that have not coalesced yet at rightParent
    //  rightCoalCount has the count of individuals that coalesced at rightParent
    //
    // intersectionParent also contains all of the above, so we'll reuse this information for that later.
    NodeIDList rightIndivs;
    rightIndivs.reserve(nextIndivs.size());
    for (const auto& indiv : nextIndivs) {
        rightIndivs.emplace_back(indiv);
        intersectIndividuals.emplace_back(indiv);
    }
    if (rightParent != INVALID_NODE_ID) {
        grg->setNumIndividualCoalsGrow(rightParent, rightCoalCount);
        nodeToIndivs.emplace(rightParent, std::move(rightIndivs));
    }

    // Now process the LEFT side individuals, which will coalescence with the RIGHT side individuals
    // at the intersectionParent.
    // There are never any coalescences in the left parent, because it has a single child.
    NodeIDSizeT intersectCoalCount = 0;
    if (grg->isSample(leftNode)) {
        const NodeIDSizeT indiv = leftNode / ploidy;
        auto insertIt = nextIndivs.insert(indiv);
        if (!insertIt.second) {
            nextIndivs.erase(insertIt.first);
            intersectCoalCount++;
        } else {
            intersectIndividuals.emplace_back(indiv);
        }
        if (leftParent != INVALID_NODE_ID) {
            nodeToIndivs.insert({leftParent, {indiv}});
            grg->setNumIndividualCoalsGrow(leftParent, 0);
        }
    } else {
        auto findIt = nodeToIndivs.find(leftNode);
        release_assert(findIt != nodeToIndivs.end());
        for (const NodeID indiv : findIt->second) {
            auto insertIt = nextIndivs.insert(indiv);
            if (!insertIt.second) {
                nextIndivs.erase(insertIt.first);
                intersectCoalCount++;
            } else {
                intersectIndividuals.emplace_back(indiv);
            }
        }
        if (leftParent != INVALID_NODE_ID) {
            nodeToIndivs.emplace(leftParent, std::move(findIt->second));
            grg->setNumIndividualCoalsGrow(leftParent, 0);
            nodeToIndivs.erase(findIt);
        }
    }
    grg->setNumIndividualCoalsGrow(intersectionParent, intersectCoalCount + rightCoalCount);
    nodeToIndivs.emplace(intersectionParent, std::move(intersectIndividuals));
}

MutableGRGPtr buildTree(HapWindowContext& context,
                        const size_t numSamples,
                        const size_t ploidy,
                        const bool isPhased,
                        GrgBuildFlags buildFlags,
                        double rebuildProportion = 0.25) {
    // Lambda that does hamming distance between nodes. Instead of passing around the vectors contains the
    // bloom filters / hashes, we use a BK-tree that has a distance callback between elements.
    auto compareNodeIds = [&](const NodeID& node1, const NodeID& node2) {
        const auto& hapVect1 = context.sampleHapVects.at(node1);
        const auto& hapVect2 = context.sampleHapVects.at(node2);
        size_t dist = 0;
        assert(hapVect1.size() == hapVect2.size());
        assert(hapVect1.size() == context.windowInfo.size());
        for (size_t i = 0; i < hapVect1.size(); i++) {
            dist += context.windowInfo[i].getDistance(hapVect1[i], hapVect2[i]);
        }
        return dist;
    };
    HaplotypeIndex hashIndex(compareNodeIds, rebuildProportion);

    NodeIDSet covered;
    NodeIDList levelNodes;

    MutableGRGPtr result = std::make_shared<MutableGRG>(numSamples, ploidy, isPhased);

    auto operationStartTime = std::chrono::high_resolution_clock::now();
    for (NodeID nodeId = 0; nodeId < numSamples; nodeId++) {
        hashIndex.add(nodeId);
        levelNodes.push_back(nodeId);
    }
    if (hasBSFlag(buildFlags, GBF_VERBOSE_OUTPUT)) {
        std::cout << "** Building index took "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                           operationStartTime)
                         .count()
                  << " ms\n";
    }
    operationStartTime = std::chrono::high_resolution_clock::now();

    /* The algorithm generates multitrees that look like this:
     *
     *              leftParent  intersectNode     rightParent
     *                 \          |              /
     *                  \        / \           / \
     *                   \---\  /   --------> B   |
     *                         A         \------> C
     *
     * Where:
     *   intersectNode represents the ancestral haplotype shared by A, B, C
     *   rightParent represents the ancestral haplotype shared by B, C, but not in intersectNode
     *   leftParent represents the ancestral haplotype of A not in intersectNode
     */

    std::unordered_map<NodeID, NodeIDList> nodeToIndivs;
    size_t level = 0;
    bool createdNodes = true;
    while (createdNodes) {
        if (hasBSFlag(buildFlags, GBF_VERBOSE_OUTPUT)) {
            std::cout << "Pass " << level << " -- " << levelNodes.size() << std::endl;
        }
        NodeIDList nextLevelNodes;

        createdNodes = false;
        while (levelNodes.size() > 1) {
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
                size_t h1OnlySize = 0;
                size_t h2OnlySize = 0;
                HapIdxList hap1 = context.sampleHapVects.at(nodeId);
                HapIdxList hap2 = context.sampleHapVects.at(first);
                release_assert(hap1.size() == hap2.size());
                release_assert(hap1.size() == context.windowInfo.size());
                HapIdxList h1Only;
                HapIdxList h2Only;
                HapIdxList intersect;

                // Generate three output HapIdxLists: intersect, h1Only, h2Only.
                for (size_t k = 0; k < hap1.size(); k++) {
                    const HapIdx h1 = hap1[k];
                    const HapIdx h2 = hap2[k];
                    if (h1 == h2) {
                        intersect.emplace_back(h1);
                        h1Only.emplace_back(0);
                        h2Only.emplace_back(0);
                        isize += context.windowInfo[k].orderedHapCounts.at(h1);
                        assert(isize > 0 || h1 == 0);
                    } else {
                        HapIdx intersectHapIdx = INVALID_HAP_IDX;
                        HapIdx h1OnlyIdx = INVALID_HAP_IDX;
                        HapIdx h2OnlyIdx = INVALID_HAP_IDX;
                        size_t intersectBits = 0;
                        const HaplotypeVector& h1Hap = context.windowInfo[k].orderedHaps.at(h1);
                        const HaplotypeVector& h2Hap = context.windowInfo[k].orderedHaps.at(h2);
                        HaplotypeVector newHapVec = bitwiseIntersect(h1Hap, h2Hap, intersectBits);
                        if (intersectBits > 0) {
                            intersectHapIdx =
                                context.windowInfo[k].addOrGetHaplotype(std::move(newHapVec), intersectBits);

                            HaplotypeVector h1OnlyHap = context.windowInfo[k].orderedHaps.at(h1);
                            const size_t h1OnlyBits =
                                bitwiseSubtract(h1OnlyHap, context.windowInfo[k].orderedHaps.at(h2));
                            if (0 == h1OnlyBits) {
                                h1OnlyIdx = 0;
                            } else {
                                h1OnlyIdx = context.windowInfo[k].addOrGetHaplotype(std::move(h1OnlyHap), h1OnlyBits);
                            }
                            h1OnlySize += h1OnlyBits;

                            HaplotypeVector h2OnlyHap = context.windowInfo[k].orderedHaps.at(h2);
                            const size_t h2OnlyBits =
                                bitwiseSubtract(h2OnlyHap, context.windowInfo[k].orderedHaps.at(h1));
                            if (0 == h2OnlyBits) {
                                h2OnlyIdx = 0;
                            } else {
                                h2OnlyIdx = context.windowInfo[k].addOrGetHaplotype(std::move(h2OnlyHap), h2OnlyBits);
                            }
                            h2OnlySize += h2OnlyBits;

                            isize += intersectBits;
                        } else {
                            h1OnlyIdx = h1;
                            h1OnlySize += context.windowInfo[k].orderedHapCounts[h1];
                            h2OnlyIdx = h2;
                            h2OnlySize += context.windowInfo[k].orderedHapCounts[h2];

                            intersectHapIdx = 0;
                        }
                        intersect.emplace_back(intersectHapIdx);
                        release_assert(h1OnlyIdx != h2OnlyIdx || (h1OnlyIdx == 0 && h2OnlyIdx == 0));
                        h1Only.emplace_back(h1OnlyIdx);
                        h2Only.emplace_back(h2OnlyIdx);
                    }
                }
                release_assert(hap1.size() == intersect.size());
                release_assert(hap1.size() == h1Only.size());
                release_assert(hap1.size() == h2Only.size());

                // We create up to three nodes: (A & B), (A - B), (B - A) for any of them that are non-empty.
                if (isize > 0) {
                    const NodeID iNodeId = result->makeNode();
                    release_assert(iNodeId >= numSamples);
                    context.sampleHapVects.emplace_back(intersect);
                    nextLevelNodes.push_back(iNodeId);
                    hashIndex.add(iNodeId);

                    result->connect(iNodeId, nodeId);
                    for (const NodeID similarId : similar) {
                        release_assert(covered.find(similarId) == covered.end());
                        result->connect(iNodeId, similarId);
                    }

                    // Add nodes for "only h1" and "only h2", i.e. things not in the intersection.
                    // We only look at this if intersection is non-zero, because otherwise "only h1" is just h1, etc.
                    NodeID h1SubNodeId = INVALID_NODE_ID;
                    if (h1OnlySize > 0) {
                        h1SubNodeId = result->makeNode();
                        context.sampleHapVects.emplace_back(h1Only);
                        result->connect(h1SubNodeId, nodeId);
                    }
                    NodeID h2SubNodeId = INVALID_NODE_ID;
                    if (h2OnlySize > 0) {
                        h2SubNodeId = result->makeNode();
                        context.sampleHapVects.emplace_back(h2Only);
                        for (const NodeID similarId : similar) {
                            result->connect(h2SubNodeId, similarId);
                        }
                    }
                    // Only the intersection node is added to the hashIndex (to become a child of another node),
                    // so we only need to save individual info for it.
                    propagateCoalInformation(result, nodeToIndivs, iNodeId, h1SubNodeId, h2SubNodeId, nodeId, similar);
                }

                for (const NodeID similarId : similar) {
                    covered.insert(similarId);
                }
                createdNodes = true;
            }
        }
        level++;
        levelNodes = std::move(nextLevelNodes);
        if (hasBSFlag(buildFlags, GBF_VERBOSE_OUTPUT)) {
            hashIndex.emitStats(std::cerr);
        }
    }

    size_t totalHaps = 0;
    for (size_t i = 0; i < context.windowInfo.size(); i++) {
        totalHaps += context.windowInfo[i].hapMap.size();
    }
    if (hasBSFlag(buildFlags, GBF_VERBOSE_OUTPUT)) {
        std::cout << "Unique haplotype segments: " << totalHaps << "\n";
        std::cout << "** Constructing tree took "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() -
                                                                           operationStartTime)
                         .count()
                  << " ms\n";
        hashIndex.emitStats(std::cerr);
    }

    // Unless requested not to, perform tree-based mutation mapping. This is cheaper than calling
    // MapMutations(), but does not create as good of hierarchy.
    if (!hasBSFlag(buildFlags, GBF_NO_TREE_MAP)) {

        // Collect the mapping from mutation to the root nodes needing it.
        std::vector<NodeIDList> mutIndexToNodes(context.allMutations.size());
        for (const NodeID rootId : result->getRootNodes()) {
            size_t haps = 0;
            assert(nodeToIndivs.find(rootId) != nodeToIndivs.end() || rootId < numSamples);
            std::vector<size_t> mutIndices = hapsToMutIndices(context.sampleHapVects[rootId], context.windowInfo);
            DEBUG_PRINT("Root " << rootId << " has " << mutIndices.size() << " mutations\n");
            for (const size_t mutIndex : mutIndices) {
                mutIndexToNodes[mutIndex].push_back(rootId);
            }
        }

        // Now create the mutation nodes and connect them. TODO: in the future we might want to call
        // MapMutations with non-sample nodes as the seeds and let it do the work here, which might
        // give us much better hierarchy.
        std::pair<BpPosition, NodeID> currentMissing = {INVALID_POSITION, INVALID_NODE_ID};
        BpPosition currentPosition = INVALID_POSITION;
        for (size_t mutIdx = 0; mutIdx < mutIndexToNodes.size(); mutIdx++) {
            const NodeID newNode = result->makeNode();
            const Mutation& mut = context.allMutations.at(mutIdx);
            const auto pos = mut.getPosition();
            if (mut.isMissing()) {
                // Any missing Mutation must be the _first_ Mutation at the given bp position.
                release_assert(currentPosition != pos);
                currentMissing = {pos, newNode};
                DEBUG_PRINT("Adding missingness node " << mut.getPosition() << " to ");
            } else {
                const NodeID missNode = (currentMissing.first == pos) ? currentMissing.second : INVALID_NODE_ID;
                result->addMutation(mut, newNode, missNode);
                DEBUG_PRINT("Adding mutation " << mut.getPosition() << ", " << mut.getAllele() << " to ");
            }
            const auto& nodeList = mutIndexToNodes[mutIdx];
            for (const NodeID childNode : nodeList) {
                DEBUG_PRINT(childNode << ", ");
                result->connect(newNode, childNode);
            }
            DEBUG_PRINT("\n");

            // Update the coalescence information for the newly create mutation node.
            if (ploidy == 2) {
                std::unordered_set<NodeIDSizeT> uncoalescedIndividuals;
                NodeIDSizeT mutNodeCoals =
                    getCoalsForParent(result, nodeToIndivs, nodeList, uncoalescedIndividuals, false);
                result->setNumIndividualCoalsGrow(newNode, mutNodeCoals);
            }
            currentPosition = mut.getPosition();
        }
    }
    return result;
}

// Just map the given mutations and their samples directly as edges. This is equivalent to the
// sparse matrix representation of these mutations+samples.
void directMap(MutableGRGPtr& grg, const std::vector<MutationAndSamples>& toBeMapped) {
    std::unordered_map<NodeID, NodeIDList> nodeToIndivs;
    for (const MutationAndSamples& pair : toBeMapped) {
        const NodeID mutNode = grg->makeNode();
        grg->addMutation(pair.mutation, mutNode);
        for (const NodeID sampleId : pair.samples) {
            grg->connect(mutNode, sampleId);
        }
        if (grg->getPloidy() == 2) {
            std::unordered_set<NodeIDSizeT> uncoalescedIndividuals;
            NodeIDSizeT mutNodeCoals =
                getCoalsForParent(grg, nodeToIndivs, pair.samples, uncoalescedIndividuals, false);
            grg->setNumIndividualCoalsGrow(mutNode, mutNodeCoals);
        }
    }
}

// Given a number of samples, compute the optimal value for the sum(unique hap segments) per tree.
// This comes from experiments on multiple datasets, and then performing a linear regression.
double getOptimalHapSegSum(const size_t numSamples, const double multiplier) {
    // We bound the result to keep it sane on outlier cases (like small test case datasets, e.g.)
    constexpr double MIN_HAPSUM = 1000;
    constexpr double MAX_HAPSUM = 20000;

    // The regression was between x=log(numSamples) and y=sum(unique hap segments), because we
    // expect the unique variants in the dataset to grow logarithmically w/r/t number of samples.
    constexpr double slope = 2383.46;
    constexpr double intercept = -18143.49;
    const double predictedValue = (slope * std::log(numSamples)) + intercept;

    // The multiplier is used to choose non-optimal values, e.g. for speeding up an analysis at
    // the cost of producing less compressed files.
    return std::min(MAX_HAPSUM, std::max(MIN_HAPSUM, predictedValue * multiplier));
}

std::string getTreeGRGName(const std::string& totalOutputName, size_t treeNumber) {
    std::stringstream ssOut;
    ssOut << totalOutputName << ".tree_" << treeNumber;
    return ssOut.str();
}

/**
 * Create a GRG from a tabular data file (IGD, BGEN, VCF, ...).
 *
 * @param[in] filePrefix A file prefix that is gauranteed to be unique in some meaningful way to the
 *      user. For example, if the user builds to GRGs from the same input, but gives them different
 *      output names, then you can use the output name as a unique prefix.
 * @param[in] sampleFile The input file containing tabular data.
 * @param[in] genomeRange The range of the genome, as defined by the input file, that this GRG will
 *      cover.
 * @param[in] hapLength TODO remove?
 * @param[in] buildFlags
 * @param[in] itFlags
 * @param[in] treeCount The number of trees to span the region with, or 0 if it should automatically
 *      detect the optimal number of trees.
 * @param[in] noTreeBelowThreshold Don't both trying to create tree hierarchy for mutations that have
 *      a count or frequency below this threshold. If this threshold is less than 1.0 then this is
 *      treated as a frequency. Otherwise it is treated as a count.
 * @param[in] indivIdToPop Map from individual identifier (string) to population identifier (string).
 * @param[in] rebuildProportion When the proportion of deleted nodes in the BK-Tree reaches this
 *      threshold, rebuild the BK-Tree.
 */
MutableGRGPtr fastGRGFromSamples(const std::string& filePrefix,
                                 const std::string& sampleFile,
                                 FloatRange& genomeRange,
                                 GrgBuildFlags buildFlags,
                                 MutationIteratorFlags itFlags,
                                 size_t treeCount,
                                 double noTreeBelowThreshold,
                                 const std::map<std::string, std::string>& indivIdToPop,
                                 const double rebuildProportion) {
    // TODO:
    // 1. Make the datatype for these haplotype segments a fixed size array. Now that we are using
    //    this for automatic parameter detection, we can't change it without extra work.
    // 2. Since length is 128, look into doing calculations using SSE/AVX.
    constexpr size_t hapLength = FIXED_HAP_LENGTH;
    static_assert(MAX_FAST_HAP_LENGTH == 254, "Fix below exception text if changing this size.");
    release_assert(hapLength <= MAX_FAST_HAP_LENGTH);

#define FAST_GRG_OUTPUT(msg)                                                                                           \
    do {                                                                                                               \
        if (hasBSFlag(buildFlags, GBF_VERBOSE_OUTPUT)) {                                                               \
            std::cerr << msg;                                                                                          \
        }                                                                                                              \
    } while (0)
#if 0
#define DUMP_HAP_SEG(hapContext)                                                                                       \
    do {                                                                                                               \
        std::cout << "Unique Hapsegments: \n";                                                                         \
        for (size_t i = 0; i < (hapContext).windowInfo.size(); i++) {                                                  \
            std::cout << i << ": " << (hapContext).windowInfo[i].hapMap.size() << "\n";                                \
        }                                                                                                              \
    } while (0)
#else
#define DUMP_HAP_SEG(hapContext)
#endif
#define CANNOT_EXIST(filename)                                                                                         \
    do {                                                                                                               \
        if (pathExists(filename)) {                                                                                    \
            std::stringstream ssErr;                                                                                   \
            ssErr << "Will not overwrite existing file: " << (filename);                                               \
            throw std::runtime_error(ssErr.str());                                                                     \
        }                                                                                                              \
    } while (0)

    std::shared_ptr<grgl::MutationIterator> mutIterator = makeMutationIterator(sampleFile, genomeRange, itFlags);

    // Algorithm:
    // 1. Get the hap segments.
    //    1b. Unphased data only: pseudo-phase the hap segments by swapping mutations.
    // 2. Create a BK-Tree that computes distance over and between hap segments.
    //    2b. Unphased data only: pseudo-phase the haplotypes by swpping hap segments.
    // 3. Build the graph using the hap lists instead of mut lists
    // 4. Map mutations to the relevant root nodes of the graph
    //
    // If we filter out LF mutations, we expect haplotype sharing to be even larger (at lower levels of the tree).
    // We could build a GRG for HF and LF separately and then merge them.
    // We could filter out LF mutations and just use MapMutations to add them (though we will likely lose all
    // kinds of hierarchy doing this). The LF mutations create more unique haplotypes, but can be useful.
    size_t numIndividuals = 0;
    bool isPhased = false;
    size_t ploidy = 0;
    mutIterator->getMetadata(ploidy, numIndividuals, isPhased);
    const size_t numSamples = numIndividuals * ploidy;
    const size_t totalMuts = mutIterator->countMutations();
    size_t remainingMuts = totalMuts;

    std::list<std::string> treeFiles;
    // Read in the haplotypes for the first tree, and then estimate the number of trees by looking at
    // the span (in number of variants) that the first tree covered vs. the total range we need all trees
    // to cover.
    if (treeCount == 0 && remainingMuts > 0) {
        const double hapSumMultiplier =
            hasBSFlag(buildFlags, GBF_TREES_FASTER1) ? 0.75 : (hasBSFlag(buildFlags, GBF_TREES_FASTER2) ? 0.5 : 1.0);
        const size_t targetHapSegSum = static_cast<size_t>(getOptimalHapSegSum(numSamples, hapSumMultiplier));
        HapWindowContext hapContext;
        getHapSegments(*mutIterator,
                       hapLength,
                       numSamples,
                       noTreeBelowThreshold,
                       /*stopAfter=*/targetHapSegSum,
                       /*stopAfterIsMutCount=*/false,
                       hapContext);
        FAST_GRG_OUTPUT("Windows: " << hapContext.windowInfo.size() << "\n");
        DUMP_HAP_SEG(hapContext);

        const size_t handledSoFar = hapContext.allMutations.size() + hapContext.directMap.size();
        release_assert(handledSoFar <= totalMuts);
        // We've already processed the first tree, so subtract 1
        treeCount = (roundUpToMultiple<size_t>(totalMuts, handledSoFar) / handledSoFar) - 1;
        release_assert(treeCount <= totalMuts); // Int underflow sanity check (should not be possible)
        FAST_GRG_OUTPUT("Auto-detected a total tree count of " << treeCount + 1 << "\n");

        MutableGRGPtr tree = buildTree(hapContext, numSamples, ploidy, isPhased, buildFlags, rebuildProportion);
        directMap(tree, hapContext.directMap);
        const std::string treeFilename = getTreeGRGName(filePrefix, 0);
        CANNOT_EXIST(treeFilename);
        saveGRG(tree, treeFilename);
        treeFiles.push_back(treeFilename);

        remainingMuts = totalMuts - handledSoFar;
    }
    release_assert(treeCount != 0 || remainingMuts == 0);

    const size_t mutsPerTree =
        (remainingMuts == 0) ? 0 : roundUpToMultiple<size_t>(remainingMuts, treeCount) / treeCount;
    for (size_t i = 1; i <= treeCount; i++) {
        HapWindowContext hapContext;
        getHapSegments(*mutIterator,
                       hapLength,
                       numSamples,
                       noTreeBelowThreshold,
                       /*stopAfter=*/mutsPerTree,
                       /*stopAfterIsMutCount=*/true,
                       hapContext);
        FAST_GRG_OUTPUT("Windows: " << hapContext.windowInfo.size() << "\n");
        DUMP_HAP_SEG(hapContext);

        MutableGRGPtr tree = buildTree(hapContext, numSamples, ploidy, isPhased, buildFlags, rebuildProportion);
        directMap(tree, hapContext.directMap);
        const std::string treeFilename = getTreeGRGName(filePrefix, i);
        CANNOT_EXIST(treeFilename);
        saveGRG(tree, treeFilename);
        treeFiles.push_back(treeFilename);
    }

#ifndef NDEBUG
    size_t _ignore_debug = 0;
    MutationAndSamples mutAndSamples;
    assert(!mutIterator->next(mutAndSamples, _ignore_debug));
#endif

#ifdef DEBUG_SANITY_CHECKS
    UniqueSamplesCheck tmp;
    result->visitDfs(tmp, TraversalDirection::DIRECTION_DOWN, result->getRootNodes());
    for (auto pair : tmp.m_failures) {
        std::cout << "FAIL: " << pair.first << " has " << pair.second << " failures\n";
    }
#endif

    MutableGRGPtr result;
    if (treeFiles.empty()) {
        result = std::make_shared<MutableGRG>(numSamples, ploidy, isPhased);
    } else {
        result = loadMutableGRG(treeFiles.front());
    }
    if (treeFiles.size() > 1) {
        std::list<std::string> toMerge = treeFiles;
        toMerge.pop_front();
        // We ignore range violations because the way we build trees above may overlap at the first/last
        // mutation of each tree.
        result->merge(
            toMerge, /*combineNodes=*/true, /*useSampleSets=*/false, /*verbose=*/false, /*ignoreRangeViolations=*/true);
    }
    for (const auto& filename : treeFiles) {
        deleteFile(filename);
    }
    addExtraInfoToGRG(result, *mutIterator, buildFlags, indivIdToPop);

#undef FAST_GRG_OUTPUT
    return result;
}

} // namespace grgl
