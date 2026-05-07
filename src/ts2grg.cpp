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
#include "grgl/ts2grg.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/node_data.h"
#include "tskit/core.h"
#include "tskit/trees.h"
#include "tskit_util.h"
#include "util.h"

#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <tskit.h>
#include <vector>

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
    explicit TsToGrgContext(const size_t numNodes, MutableGRGPtr grg)
        : m_argToGrg(numNodes, INVALID_NODE_ID),
          m_pendingSplit(numNodes),
          m_grg(std::move(grg)) {}

    void addNodeMapping(tsk_id_t tsNodeId, NodeID grgNodeId, NodeIDList uncoalesced = {}) {
        // Once a GRG node is no longer on the "frontier", we compact the edge lists.
        if (m_argToGrg.at(tsNodeId) != INVALID_NODE_ID) {
            const NodeID oldGrgNode = m_argToGrg[tsNodeId];
            m_grg->compact(oldGrgNode);
        }
        m_argToGrg[tsNodeId] = grgNodeId;
        if (!uncoalesced.empty()) {
            const auto insertResult = m_nodeToUncoalesced.emplace(grgNodeId, std::move(uncoalesced));
            assert(insertResult.second);
        }
    }

    NodeID getLatestMappedNode(tsk_id_t tsNodeId) { return m_argToGrg.at(tsNodeId); }

    const NodeIDList& getIndividualCoals(const NodeID grgNode) {
        const static NodeIDList empty = {};
        const auto coalIt = m_nodeToUncoalesced.find(grgNode);
        if (coalIt != m_nodeToUncoalesced.end()) {
            return coalIt->second;
        }
        return empty;
    }

    void countCoals(NodeID parent, const NodeIDList& children, bool useExisting) {
        release_assert(!children.empty());
        // Fast paths: just copy, no coalescences can happen.
        if (children.size() == 1 and !useExisting) {
            m_nodeToUncoalesced.emplace(parent, this->getIndividualCoals(children.front()));
            m_grg->setNumIndividualCoals(parent, 0);
            return;
        }

        NodeIDSizeT coalCount = 0;
        NodeIDSet uncoalescedSet;
        if (useExisting) {
            const NodeIDSizeT existing = m_grg->getNumIndividualCoals(parent);
            if (existing != COAL_COUNT_NOT_SET) {
                coalCount += existing;
            }
            const auto coalIt = m_nodeToUncoalesced.find(parent);
            if (coalIt != m_nodeToUncoalesced.end()) {
                for (NodeID indiv : coalIt->second) {
                    uncoalescedSet.emplace(indiv);
                }
            }
        }
        for (NodeID childId : children) {
            for (NodeID indivId : this->getIndividualCoals(childId)) {
                const auto insertPair = uncoalescedSet.emplace(indivId);
                if (!insertPair.second) {
                    coalCount++; // Coalescence!
                    uncoalescedSet.erase(insertPair.first);
                }
            }
        }
        if (!uncoalescedSet.empty()) {
            const auto insertResult = m_nodeToUncoalesced.emplace(parent, NodeIDList());
            for (NodeID indivId : uncoalescedSet) {
                insertResult.first->second.emplace_back(indivId);
            }
        }
        m_grg->setNumIndividualCoals(parent, coalCount);
    }

    inline bool isPendingSplit(const tsk_id_t tsNodeId) {
        if (tsNodeId < m_pendingSplit.size()) {
            return m_pendingSplit[tsNodeId];
        }
        return false;
    }

    inline bool hasExistingNode(const tsk_id_t tsNodeId) { return (m_argToGrg.at(tsNodeId) != INVALID_NODE_ID); }

    inline void removePendingSplit(const tsk_id_t tsNodeId) {
        if (tsNodeId < m_pendingSplit.size()) {
            m_pendingSplit[tsNodeId] = false;
        }
    }

    inline void addPendingSplit(const tsk_id_t tsNodeId) {
        m_pendingSplit.at(tsNodeId) = true;
        // When we mark a node for splitting, the coalescence information at that node (as far as
        // which samples have not yet coalesced) is invalidated.
        if (m_argToGrg.at(tsNodeId) != INVALID_NODE_ID) {
            auto coalIt = m_nodeToUncoalesced.find(m_argToGrg.at(tsNodeId));
            if (coalIt != m_nodeToUncoalesced.end()) {
                m_nodeToUncoalesced.erase(coalIt);
            }
        }
    }

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
        bool isSampleA = true;
        bool isSampleB = true;
        if (nodeA != TSK_NULL) {
            if (isSampleA) {
                isSampleA = (bool)tsk_tree_is_sample(tree, nodeA);
            }
            if (!isSampleA) {
                addPendingSplit(nodeA);
            }
            TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeA, &timeA), "Getting node time");
        }
        if (nodeB != TSK_NULL) {
            if (isSampleB) {
                isSampleB = (bool)tsk_tree_is_sample(tree, nodeB);
            }
            if (!isSampleB) {
                addPendingSplit(nodeB);
            }
            TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeB, &timeB), "Getting node time");
        }
        while (nodeA != nodeB) {
            if (nodeA != TSK_NULL && (nodeB == TSK_NULL || timeA < timeB)) {
                TSKIT_OK_OR_THROW(tsk_tree_get_parent(tree, nodeA, &nodeA), "Getting parent");
                // Can occur if there are multiple roots to the tree.
                if (TSK_NULL != nodeA) {
                    TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeA, &timeA), "Getting node time");
                    if (isSampleA) {
                        isSampleA = (bool)tsk_tree_is_sample(tree, nodeA);
                    }
                    if (!isSampleA) {
                        addPendingSplit(nodeA);
                    }
                }
            } else {
                TSKIT_OK_OR_THROW(tsk_tree_get_parent(tree, nodeB, &nodeB), "Getting parent");
                // Can occur if there are multiple roots to the tree.
                if (TSK_NULL != nodeB) {
                    TSKIT_OK_OR_THROW(tsk_tree_get_time(tree, nodeB, &timeB), "Getting node time");
                    if (isSampleB) {
                        isSampleB = (bool)tsk_tree_is_sample(tree, nodeB);
                    }
                    if (!isSampleB) {
                        addPendingSplit(nodeB);
                    }
                }
            }
        }
        release_assert(nodeA == nodeB);
        return nodeA;
    }

    void markPendingUntilMRCA(const tsk_tree_t* tree, const std::vector<tsk_id_t>& mrcaLeaves) {
        std::vector<tsk_id_t> workList;
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

    void markPendingUntilRoot(const tsk_tree_t* tree, const std::vector<tsk_id_t>& leaves) {
        std::unordered_set<tsk_id_t> seen;
        for (const tsk_id_t node : leaves) {
            bool isSample = true;
            tsk_id_t activeNode = node;
            while (activeNode != TSK_NULL) {
                if (isSample) {
                    isSample = (bool)tsk_tree_is_sample(tree, activeNode);
                }
                if (!isSample) {
                    addPendingSplit(activeNode);
                }
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
    MutableGRGPtr m_grg;
    // Map from tskit ID to GRG IDs, for nodes.
    NodeIDList m_argToGrg;
    // Map from GRG node to a list of uncoalesced individuals at that node.
    std::unordered_map<NodeID, NodeIDList> m_nodeToUncoalesced;
    // Position "i" being true means that tskit node with ID "i" will be duplicated on
    // the next use.
    std::vector<bool> m_pendingSplit;
};

static NodeID duplicateTopoAbove(TsToGrgContext& context,
                                 MutableGRGPtr& grg,
                                 const tsk_tree_t* tree,
                                 const std::unordered_set<tsk_id_t>& tsChildMutNodes,
                                 const tsk_id_t tsParentMutNode,
                                 const bool computeCoals) {
    release_assert(!tsChildMutNodes.empty());
    release_assert(tsParentMutNode != TSK_NULL);

    // We duplicate every node between the children and the parent Mutation. There can be multiple
    // children, and they can share paths to the parent, so we need to make sure we only dup each
    // node once.
    std::unordered_map<tsk_id_t, NodeID> dupedNodes;
    auto getDupedGRGNode = [&](tsk_id_t tsNode, bool& create) {
        auto dupIt = dupedNodes.find(tsNode);
        if (dupIt != dupedNodes.end()) {
            create = false;
            return dupIt->second;
        }
        if (create) {
            const NodeID dupGrgNode = grg->makeNode();
            dupedNodes.emplace(tsNode, dupGrgNode);
            return dupGrgNode;
        }
        return INVALID_NODE_ID;
    };

    /* We duplicate nodes between children and parent, and then connect to the siblings of the nodes
       along this path.  Example:
                                        a  (parent)
                                      /   \
                          (child 1)  b     c
                                         /   \
                                        d     e  (child 2)
       Nodes on the path: c, a
       Siblings to connect: d

       The topology stays exactly the same as shown for the child mutations purposes (and nodes above
       "a"). However, for the parent "a" we create an alternate topology off to the side. In this alternate
       topology, the edges a->b and c->e do not exist, so the only coalescence counting we need to do is
       from "d" upwards. More generally: only collect coalescences from siblings and propagate them upward.
    */

    // Pass 1: Duplicate all nodes between the children and the parent. We have to do this first, because
    // when there are multiple children they can have arbitrary positions on the tree and we don't know how
    // they relate to each other. We only create edges between duplicated nodes, not other nodes.
    for (tsk_id_t tsChildMutNode : tsChildMutNodes) {
        // When the nodes are the same, the parent mutation has no effect and we can just drop it.
        if (tsChildMutNode == tsParentMutNode) {
            continue;
        }
        bool prevJustCreated = false;
        NodeID prevDupNode = INVALID_NODE_ID;
        for (tsk_id_t node = tree->parent[tsChildMutNode]; node != TSK_NULL; node = tree->parent[node]) {
            bool create = true;
            const NodeID dupNode = getDupedGRGNode(node, create);
            release_assert(dupNode != INVALID_NODE_ID);
            // If the previously duplicated node on this path was just created, then add an edge to it.
            if (prevDupNode != INVALID_NODE_ID && prevJustCreated) {
                grg->connect(dupNode, prevDupNode);
            }
            prevDupNode = dupNode;
            prevJustCreated = create;
            if (node == tsParentMutNode) {
                break;
            }
        }
    }

    // Pass 2: Now we have all the necessary duplicated nodes, we can connect them up any non-duplicated
    // nodes that we need (namely: siblings to the child mutations or duplicated nodes)
    NodeIDSet done;
    bool shouldCreate = false;
    for (tsk_id_t tsChildMutNode : tsChildMutNodes) {
        // When the nodes are the same, the parent mutation has no effect and we can just drop it.
        if (tsChildMutNode == tsParentMutNode) {
            continue;
        }

        // We traverse up the tskit tree, instead of the GRG (which is a multi-tree that may reflect multiple
        // trees up to this point). This lets us restrict our search to a single path upwards.
        for (tsk_id_t node = tree->parent[tsChildMutNode]; node != TSK_NULL; node = tree->parent[node]) {
            // This will return the duplicated node if present, or INVALID_NODE_ID otherwise.
            release_assert(shouldCreate == false);
            const NodeID dupNode = getDupedGRGNode(node, shouldCreate);
            release_assert(dupNode != INVALID_NODE_ID);

            // Only add edges once.
            if (done.find(dupNode) == done.end()) {
                // We want to create edges between the duplicated node and any non-duplicated node that does
                // not correspond to a child Mutation. The duplicated nodes represent all paths from the children
                // mutations to the parent mutation.
                for (tsk_id_t nonDupChild = tree->left_child[node]; nonDupChild != TSK_NULL;
                     nonDupChild = tree->right_sib[nonDupChild]) {
                    // Skip child mutation nodes.
                    if (tsChildMutNodes.find(nonDupChild) != tsChildMutNodes.end()) {
                        continue;
                    }
                    // Skip edges to already-duplicated nodes
                    release_assert(shouldCreate == false);
                    const NodeID checkForDup = getDupedGRGNode(nonDupChild, shouldCreate);
                    if (checkForDup != INVALID_NODE_ID) {
                        continue;
                    }
                    const NodeID nonDupChildGRG = context.getLatestMappedNode(nonDupChild);
                    grg->connect(dupNode, nonDupChildGRG);
                }
                done.emplace(dupNode);
            }
            // If we're computing coalescences, we have to do it in an "update" fashion, since each node
            // could be visited >1 time (we don't have a topological order...)
            if (computeCoals) {
                context.countCoals(dupNode, grg->getDownEdges(dupNode), true);
            }
            if (node == tsParentMutNode) {
                break;
            }
        }
    }
    release_assert(shouldCreate == false);
    return getDupedGRGNode(tsParentMutNode, shouldCreate);
}

// Given recurrent mutations, combine the nodes into a single parent node and a single Mutation
static NodeID combineMutations(TsToGrgContext& context,
                               MutableGRGPtr& grg,
                               const NodeIDList& mutNodes,
                               const Mutation& mutation,
                               const bool computeCoals) {
    release_assert(mutNodes.size() >= 2);
    const NodeID newNode = grg->makeNode();
    for (const auto& node : mutNodes) {
        if (node != INVALID_NODE_ID) {
            grg->connect(newNode, node);
        }
    }
    grg->addMutation(mutation, newNode);
    if (computeCoals) {
        context.countCoals(newNode, mutNodes, false);
    }
    return newNode;
}

// TODO make this non-recursive.
static NodeID addMutationFromTree(TsToGrgContext& context,
                                  MutableGRGPtr& grg,
                                  const tsk_treeseq_t* treeSeq,
                                  const tsk_tree_t* tree,
                                  const tsk_id_t tsNodeId,
                                  const bool computeCoals = false) {
    const bool isSample = tsk_treeseq_is_sample(treeSeq, tsNodeId);
    const bool canReuse = context.hasExistingNode(tsNodeId) && (isSample || !context.isPendingSplit(tsNodeId));
    if (canReuse) {
        return context.getLatestMappedNode(tsNodeId);
    }
    context.removePendingSplit(tsNodeId);
    release_assert(!isSample);

    NodeIDSizeT coalCount = 0;
    NodeIDSet uncoalescedSet;
    NodeID newNodeId = grg->makeNode();
    for (tsk_id_t childNodeId = tree->left_child[tsNodeId]; childNodeId != TSK_NULL;
         childNodeId = tree->right_sib[childNodeId]) {
        const NodeID childId = addMutationFromTree(context, grg, treeSeq, tree, childNodeId, computeCoals);
        grg->connect(newNodeId, childId);
        if (computeCoals) {
            for (NodeID indivId : context.getIndividualCoals(childId)) {
                const auto insertPair = uncoalescedSet.emplace(indivId);
                if (!insertPair.second) {
                    coalCount++; // Coalescence!
                    uncoalescedSet.erase(insertPair.first);
                }
            }
        }
    }
    NodeIDList uncoalesced;
    if (computeCoals) {
        for (NodeID indivId : uncoalescedSet) {
            uncoalesced.emplace_back(indivId);
        }
        grg->setNumIndividualCoals(newNodeId, coalCount);
    }
    context.addNodeMapping(tsNodeId, newNodeId, std::move(uncoalesced));
    return newNodeId;
}

using MutAndTSNode = std::pair<Mutation, tsk_id_t>;

std::vector<MutAndTSNode> getMutationsForSite(const tsk_tree_t* currentTree,
                                              const tsk_site_t* site,
                                              const bool binaryMutations,
                                              const bool useNodeTimes) {
    std::vector<MutAndTSNode> result;
    const std::string ancestralState = std::string(site->ancestral_state, site->ancestral_state_length);
    for (tsk_size_t j = 0; j < site->mutations_length; j++) {
        const tsk_id_t tsMutNode = site->mutations[j].node;
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
        }
        double mutTime = site->mutations[j].time;
        if (useNodeTimes) {
            TSKIT_OK_OR_THROW(tsk_tree_get_time(currentTree, tsMutNode, &mutTime), "Failed to get node time");
        }
        Mutation theMutation =
            Mutation((uint64_t)site->position, std::move(derivedState), ancestralState, (uint32_t)mutTime);
        result.emplace_back(std::move(theMutation), tsMutNode);
    }
    return std::move(result);
}

MutableGRGPtr convertTreeSeqToGRG(const tsk_treeseq_t* treeSeq,
                                  bool binaryMutations,
                                  bool useNodeTimes,
                                  bool maintainTopology,
                                  bool computeCoals,
                                  std::pair<size_t, size_t> treeRange) {
    const size_t initialNodeCount = tsk_treeseq_get_num_nodes(treeSeq);
    const size_t numSamples = tsk_treeseq_get_num_samples(treeSeq);
    const size_t numIndividuals = tsk_treeseq_get_num_individuals(treeSeq);
    const size_t ploidy = (numIndividuals == 0) ? 2 : (numSamples / numIndividuals);
    MutableGRGPtr grg = std::make_shared<MutableGRG>(numSamples, static_cast<uint16_t>(ploidy), initialNodeCount);
    tsk_tree_t currentTree;
    TsToGrgContext constructionContext(initialNodeCount, grg);

    if (treeRange.first == treeRange.second) {
        treeRange.second = std::numeric_limits<size_t>::max();
    }

    // Copy over the population data from the samples.
    const tsk_size_t numPopulations = tsk_treeseq_get_num_populations(treeSeq);
    for (size_t i = 0; i < numPopulations; i++) {
        tsk_population_t population;
        TSKIT_OK_OR_THROW(tsk_treeseq_get_population(treeSeq, i, &population), "Failed getting population");
        grg->addPopulation(std::string(population.metadata, population.metadata_length));
    }

    // We have to do some extra sanity checks when computing coalescences.
    std::vector<tsk_id_t> sampleToIndiv(numSamples, TSK_NULL);
    if (computeCoals) {
        if (ploidy != 2) {
            throw ApiMisuseFailure("Can only compute individual coalescences with diploid data.");
        }
        for (tsk_id_t i = 0; i < numIndividuals; i++) {
            tsk_individual_t individual;
            TSKIT_OK_OR_THROW(tsk_treeseq_get_individual(treeSeq, i, &individual), "Failed getting individuals");
            release_assert(individual.nodes_length == 2);
            if ((individual.nodes[0] / 2) != (individual.nodes[1] / 2)) {
                throw ApiMisuseFailure("Cannot convert coalescences if tree-seq does not have ordered individuals");
            }
            sampleToIndiv.at(individual.nodes[0]) = i;
            sampleToIndiv.at(individual.nodes[1]) = i;
        }
    }

    // Setup the GRG sample nodes to have the first block of IDs.
    tsk_id_t maxSample = 0;
    const tsk_id_t* sampleMap = tsk_treeseq_get_samples(treeSeq);
    for (size_t i = 0; i < numSamples; i++) {
        const tsk_id_t tsNodeId = sampleMap[i];
        if (tsNodeId > maxSample) {
            maxSample = tsNodeId;
        }
        tsk_node_t node;
        TSKIT_OK_OR_THROW(tsk_treeseq_get_node(treeSeq, tsNodeId, &node), "Failed getting node");
        NodeIDList individual;
        if (computeCoals) {
            individual.emplace_back(sampleToIndiv.at(tsNodeId));
        }
        constructionContext.addNodeMapping(tsNodeId, i, std::move(individual));

        if (node.population != TSK_NULL) {
            grg->setPopulationId((NodeID)i, node.population);
        }
    }
    // GRG can theoretically store datasets where the last few "generations" of the trees are samples, and not just
    // the leaves, but it is untested. This checks that the K samples are the first K node IDs.
    if (maxSample >= numSamples) {
        throw ApiMisuseFailure("tree-sequences with non-consecutive, non-leafy samples cannot be converted to GRG.");
    }

    // FIXME we could probably save about 4x the RAM if we converted the edge spans into tree indices
    // and used a dense vector from tree index to node information.
    // I think the main RAM usage right now is the MutableGRG itself, because we end up creating a
    // ton of Nodes that are then simplified later. The best way to solve this is probably to do the
    // TS->GRG conversion is sections (by genome range) and then merge them after simplification.
    // This would also allow for parallel conversion.

    // Compute which nodes are involved in edge additions/deletions for each interval.
    std::unordered_map<double, std::vector<tsk_id_t>> positionToModifiedNodes;
    for (size_t i = 0; i < treeSeq->tables->edges.num_rows; i++) {
        const double edgeAddedPosition = treeSeq->tables->edges.left[i];
        const double edgeDeletedPosition = treeSeq->tables->edges.right[i];
        const tsk_id_t parent = treeSeq->tables->edges.parent[i];
        const tsk_id_t child = treeSeq->tables->edges.child[i];
        auto addIt = positionToModifiedNodes.insert(
            std::pair<double, std::vector<tsk_id_t>>(edgeAddedPosition, std::vector<tsk_id_t>()));
        addIt.first->second.push_back(parent);
        addIt.first->second.push_back(child);
        auto deleteIt = positionToModifiedNodes.insert(
            std::pair<double, std::vector<tsk_id_t>>(edgeDeletedPosition, std::vector<tsk_id_t>()));
        deleteIt.first->second.push_back(parent);
        deleteIt.first->second.push_back(child);
    }

    size_t treeCount = 0;
    TSKIT_OK_OR_THROW(tsk_tree_init(&currentTree, treeSeq, 0), "Failed to init tree");
    int treeItResult = 0;
    for (treeItResult = tsk_tree_first(&currentTree); TSK_TREE_OK == treeItResult;
         treeItResult = tsk_tree_next(&currentTree)) {
        if (treeCount >= treeRange.first && treeCount < treeRange.second) {
            const double endOfTree = currentTree.interval.right;
            // Add all the new mutations to the GRG, and the relevant nodes/edges.
            tsk_size_t num_sites = 0;
            const tsk_site_t* sites = nullptr;
            int siteItResult = tsk_tree_get_sites(&currentTree, &sites, &num_sites);
            TSKIT_OK_OR_THROW(siteItResult, "Failed to get sites from tree");
            for (tsk_size_t i = 0; i < num_sites; i++) {
                const tsk_site_t* site = &sites[i];

                // Copy the tskit nodes and tree paths related to Mutations into the GRG. Don't add the Mutations
                // yet (see next step).
                NodeIDList mutNodes;
                std::vector<MutAndTSNode> muts = getMutationsForSite(&currentTree, site, binaryMutations, useNodeTimes);
                for (const auto& mutAndNode : muts) {
                    const NodeID mutationNodeId = addMutationFromTree(
                        constructionContext, grg, treeSeq, &currentTree, mutAndNode.second, computeCoals);
                    mutNodes.emplace_back(mutationNodeId);
                }
                // Duplicate the topology as needed for "nested" (back) Mutations, since the presence of
                // the Mutation beneath another Mutation acts as "set subtraction" for the sample set
                // beneath. The tskit way of storing this makes calculations complex, and GRG requires
                // something simpler.
                std::unordered_map<Mutation, NodeIDList> mutToNodes;
                for (tsk_size_t j = 0; j < site->mutations_length; j++) {
                    NodeID mutationNodeId = mutNodes.at(j);
                    const tsk_id_t tsParentMutNode = site->mutations[j].node;
                    // Find the child (if any) of current mutation
                    std::unordered_set<tsk_id_t> tsChildMutNodes;
                    for (size_t k = 0; k < site->mutations_length; k++) {
                        if (site->mutations[k].parent == site->mutations[j].id) {
                            tsChildMutNodes.emplace(site->mutations[k].node);
                        }
                    }
                    if (!tsChildMutNodes.empty()) {
                        mutationNodeId = duplicateTopoAbove(
                            constructionContext, grg, &currentTree, tsChildMutNodes, tsParentMutNode, computeCoals);
                    }
                    // We skip mutations where the REF == ALT, because we have already accounted for any
                    // "subtractive effect" they had in the original tree.
                    Mutation& mutation = muts.at(j).first;
                    if (mutation.getAllele() != mutation.getRefAllele()) {
                        auto mutNodeIt = mutToNodes.emplace(std::move(mutation), NodeIDList());
                        mutNodeIt.first->second.emplace_back(mutationNodeId);
                    }
                }
                // We do this at the end, because we want the opportunity to create a SINGLE Mutation for potentially
                // multiple separate tskit mutations.
                for (const auto& mutAndNodes : mutToNodes) {
                    const Mutation& mutation = mutAndNodes.first;
                    if (mutAndNodes.second.size() == 1) {
                        grg->addMutation(mutation, mutAndNodes.second.front());
                    } else {
                        combineMutations(constructionContext, grg, mutAndNodes.second, mutation, computeCoals);
                    }
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
                positionToModifiedNodes.erase(modifiedNodesIt);
            } else {
                std::cerr << "Missing modification between trees at " << endOfTree << std::endl;
            }
        }
        treeCount++;
    }
    TSKIT_OK_OR_THROW(treeItResult, "Failed while iterating tree-sequence");
    tsk_tree_free(&currentTree);
    grg->sortMutations();
    return grg;
}

} // namespace grgl
