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
#include "grgl/ts2grg.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "tskit/core.h"
#include "tskit/trees.h"
#include "tskit_util.h"
#include "util.h"

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <tskit.h>

#define TSKIT_OK_OR_THROW(ok, msg)                                                                                     \
    do {                                                                                                               \
        int tskit_ok_val = (ok);                                                                                       \
        if (tskit_ok_val != 0) {                                                                                       \
            std::stringstream errMsg;                                                                                  \
            errMsg << (msg) << ": " << tsk_strerror(tskit_ok_val);                                                     \
            throw TskitApiFailure(msg);                                                                                \
        }                                                                                                              \
    } while (0)

namespace grgl {

class TsToGrgContext {
public:
    explicit TsToGrgContext(size_t numNodes)
        : m_pendingSplit(numNodes) {}

    void addNodeMapping(tsk_id_t tsNodeId, NodeID grgNodeId) {
        auto resultPair =
            this->m_argToGrg.insert(std::pair<tsk_id_t, std::list<NodeID>>(tsNodeId, std::list<NodeID>()));
        resultPair.first->second.push_back(grgNodeId);
    }

    NodeID getLatestMappedNode(tsk_id_t tsNodeId) {
        auto findIt = this->m_argToGrg.find(tsNodeId);
        return findIt->second.back();
    }

    inline bool isPendingSplit(const tsk_id_t tsNodeId) {
        if (tsNodeId < m_pendingSplit.size()) {
            return m_pendingSplit[tsNodeId];
        }
        return false;
    }

    inline bool hasExistingNode(const tsk_id_t tsNodeId) {
        auto itExists = m_argToGrg.find(tsNodeId);
        return itExists != m_argToGrg.end();
    }

    inline void removePendingSplit(const tsk_id_t tsNodeId) {
        if (tsNodeId < m_pendingSplit.size()) {
            m_pendingSplit[tsNodeId] = false;
        }
    }

    inline void addPendingSplit(const tsk_id_t tsNodeId) { m_pendingSplit.at(tsNodeId) = true; }

    /**
     * Find the most recent common ancester in a forest with times on nodes.
     * Either (or both) of the two nodes can be NULL, which implies the (non-existent)
     * node above the root in the tree.
     */
    tsk_id_t markPendingAndGetMRCA(const tsk_tree_t* tree, const tsk_id_t origNodeA, const tsk_id_t origNodeB) {
        tsk_id_t nodeA = origNodeA;
        tsk_id_t nodeB = origNodeB;
        double timeA = -1.0;
        double timeB = -1.0;
        if (nodeA != TSK_NULL) {
            addPendingSplit(nodeA);
            TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeA, &timeA), "Getting node time");
        }
        if (nodeB != TSK_NULL) {
            addPendingSplit(nodeB);
            TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeB, &timeB), "Getting node time");
        }
        while (nodeA != nodeB) {
            if (nodeA != TSK_NULL && (nodeB == TSK_NULL || timeA < timeB)) {
                TSKIT_OK_OR_THROW(tsk_tree_get_parent(tree, nodeA, &nodeA), "Getting parent");
                // Can occur if there are multiple roots to the tree.
                if (TSK_NULL != nodeA) {
                    TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeA, &timeA), "Getting node time");
                    addPendingSplit(nodeA);
                }
            } else {
                TSKIT_OK_OR_THROW(tsk_tree_get_parent(tree, nodeB, &nodeB), "Getting parent");
                // Can occur if there are multiple roots to the tree.
                if (TSK_NULL != nodeB) {
                    TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeB, &timeB), "Getting node time");
                    addPendingSplit(nodeB);
                }
            }
        }
        release_assert(nodeA == nodeB);
        return nodeA;
    }

    void markPendingUntilMRCA(const tsk_tree_t* tree, const std::list<tsk_id_t>& mrcaLeaves) {
        std::list<tsk_id_t> workList;
        for (const tsk_id_t node : mrcaLeaves) {
            if (!tsIsDisconnected(tree, node)) {
                workList.push_back(node);
            }
        }
        while (workList.size() >= 2) {
            const tsk_id_t nodeA = workList.back();
            workList.pop_back();
            const tsk_id_t nodeB = workList.back();
            workList.pop_back();
            const tsk_id_t ancestor = markPendingAndGetMRCA(tree, nodeA, nodeB);
            workList.push_back(ancestor);
        }
    }

    void markPendingUntilRoot(const tsk_tree_t* tree, const std::list<tsk_id_t>& leaves) {
        std::unordered_set<tsk_id_t> seen;
        for (const tsk_id_t node : leaves) {
            tsk_id_t activeNode = node;
            while (activeNode != TSK_NULL) {
                addPendingSplit(activeNode);
                auto insertPair = seen.emplace(activeNode);
                // We didn't actually insert b/c it was already there -- we can stop traversing.
                if (!insertPair.second) {
                    break;
                }
                TSKIT_OK_OR_THROW(tsk_tree_get_parent(tree, activeNode, &activeNode), "Getting parent");
            }
        }
    }

private:
    std::unordered_map<tsk_id_t, std::list<NodeID>> m_argToGrg;
    std::vector<bool> m_pendingSplit;
};

// TODO make this non-recursive.
static NodeID addMutationFromTree(TsToGrgContext& context,
                                  MutableGRGPtr& grg,
                                  const tsk_treeseq_t* treeSeq,
                                  const tsk_tree_t* tree,
                                  tsk_id_t tsNodeId) {
    const bool isSample = tsk_treeseq_is_sample(treeSeq, tsNodeId);
    const bool canReuse = context.hasExistingNode(tsNodeId) && (isSample || !context.isPendingSplit(tsNodeId));
    if (canReuse) {
        return context.getLatestMappedNode(tsNodeId);
    }
    context.removePendingSplit(tsNodeId);
    release_assert(!isSample);

    NodeID newNodeId = grg->makeNode();
    for (tsk_id_t childNodeId = tree->left_child[tsNodeId]; childNodeId != TSK_NULL;
         childNodeId = tree->right_sib[childNodeId]) {
        const NodeID childId = addMutationFromTree(context, grg, treeSeq, tree, childNodeId);
        grg->connect(newNodeId, childId);
    }
    context.addNodeMapping(tsNodeId, newNodeId);
    return newNodeId;
}

MutableGRGPtr
convertTreeSeqToGRG(const tsk_treeseq_t* treeSeq, bool binaryMutations, bool useNodeTimes, bool maintainTopology) {
    const size_t initialNodeCount = tsk_treeseq_get_num_nodes(treeSeq);
    const size_t numSamples = tsk_treeseq_get_num_samples(treeSeq);
    const size_t numIndividuals = tsk_treeseq_get_num_individuals(treeSeq);
    const size_t ploidy = (numIndividuals == 0) ? 2 : (numSamples / numIndividuals);
    MutableGRGPtr grg = std::make_shared<MutableGRG>(numSamples, static_cast<uint16_t>(ploidy), initialNodeCount);
    tsk_tree_t currentTree;
    TsToGrgContext constructionContext(initialNodeCount);

    // Copy over the population data from the samples.
    const tsk_size_t numPopulations = tsk_treeseq_get_num_populations(treeSeq);
    for (size_t i = 0; i < numPopulations; i++) {
        tsk_population_t population;
        TSKIT_OK_OR_THROW(tsk_treeseq_get_population(treeSeq, i, &population), "Failed getting population");
        grg->addPopulation(std::string(population.metadata, population.metadata_length));
    }
    // Setup the GRG sample nodes to have the first block of IDs.
    const tsk_id_t* sampleMap = tsk_treeseq_get_samples(treeSeq);
    for (size_t i = 0; i < numSamples; i++) {
        const tsk_id_t tsNodeId = sampleMap[i];
        tsk_node_t node;
        TSKIT_OK_OR_THROW(tsk_treeseq_get_node(treeSeq, tsNodeId, &node), "Failed getting node");
        constructionContext.addNodeMapping(tsNodeId, i);

        auto& grgNodeData = grg->getNodeData((NodeID)i);
        if (node.population != TSK_NULL) {
            grgNodeData.populationId = node.population;
        }
    }

    // Compute which nodes are involved in edge additions/deletions for each interval.
    std::unordered_map<double, std::list<tsk_id_t>> positionToModifiedNodes;
    for (size_t i = 0; i < treeSeq->tables->edges.num_rows; i++) {
        const double edgeAddedPosition = treeSeq->tables->edges.left[i];
        const double edgeDeletedPosition = treeSeq->tables->edges.right[i];
        auto addIt = positionToModifiedNodes.insert(
            std::pair<double, std::list<tsk_id_t>>(edgeAddedPosition, std::list<tsk_id_t>()));
        addIt.first->second.push_back(treeSeq->tables->edges.parent[i]);
        addIt.first->second.push_back(treeSeq->tables->edges.child[i]);
        auto deleteIt = positionToModifiedNodes.insert(
            std::pair<double, std::list<tsk_id_t>>(edgeDeletedPosition, std::list<tsk_id_t>()));
        deleteIt.first->second.push_back(treeSeq->tables->edges.parent[i]);
        deleteIt.first->second.push_back(treeSeq->tables->edges.child[i]);
    }

    size_t treeCount = 0;
    TSKIT_OK_OR_THROW(tsk_tree_init(&currentTree, treeSeq, 0), "Failed to init tree");
    int treeItResult = 0;
    for (treeItResult = tsk_tree_first(&currentTree); TSK_TREE_OK == treeItResult;
         treeItResult = tsk_tree_next(&currentTree)) {
        const double startOfTree = currentTree.interval.left;
        const double endOfTree = currentTree.interval.right;
        // Add all the new mutations to the GRG, and the relevant nodes/edges.
        tsk_size_t num_sites = 0;
        const tsk_site_t* sites = nullptr;
        int siteItResult = tsk_tree_get_sites(&currentTree, &sites, &num_sites);
        TSKIT_OK_OR_THROW(siteItResult, "Failed to get sites from tree");
        for (tsk_size_t i = 0; i < num_sites; i++) {
            const tsk_site_t* site = &sites[i];
            std::string ancestralState = std::string(site->ancestral_state, site->ancestral_state_length);
            for (tsk_size_t j = 0; j < site->mutations_length; j++) {
                const tsk_id_t tsMutNode = site->mutations[j].node;
                EXTERNAL_ID originalId = NO_ORIGINAL_ID;
                std::string derivedState =
                    std::string(site->mutations[j].derived_state, site->mutations[j].derived_state_length);
                if (binaryMutations) {
                    if (derivedState != ancestralState) {
                        derivedState = Mutation::ALLELE_1;
                    } else {
                        // In the case that the tree sequence has recurrent mutations back to the reference
                        // allele, we need to still capture that. This makes comparison with other GRGs,
                        // such as from VCF, a pain, but there isn't a great alternative.
                        derivedState = Mutation::ALLELE_0;
                    }
                } else {
                    originalId = site->mutations[j].id;
                }
                double mutTime = site->mutations[j].time;
                if (useNodeTimes) {
                    TSKIT_OK_OR_THROW(tsk_tree_get_time(&currentTree, tsMutNode, &mutTime), "Failed to get node time");
                }
                const Mutation theMutation = Mutation(
                    (uint64_t)site->position, std::move(derivedState), ancestralState, (float)mutTime, originalId);
                // Clients can pass a predicate that only includes certain mutations.
                NodeID mutationNodeId = addMutationFromTree(constructionContext, grg, treeSeq, &currentTree, tsMutNode);
                grg->addMutation(theMutation, mutationNodeId);
            }
        }

        // Now get the set of nodes that will be affected by edge additions/deletions between
        // the current tree and the next one. We then take the MRCA of these nodes and mark
        // everything in between as needing-to-be-split.
        const auto modifiedNodesIt = positionToModifiedNodes.find(endOfTree);
        if (modifiedNodesIt != positionToModifiedNodes.end()) {
            if (maintainTopology) {
                constructionContext.markPendingUntilRoot(&currentTree, modifiedNodesIt->second);
            } else {
                constructionContext.markPendingUntilMRCA(&currentTree, modifiedNodesIt->second);
            }
        } else {
            std::cerr << "Missing modification between trees at " << endOfTree << std::endl;
        }
        treeCount++;
    }
    TSKIT_OK_OR_THROW(treeItResult, "Failed while iterating tree-sequence");
    tsk_tree_free(&currentTree);
    return grg;
}

} // namespace grgl
