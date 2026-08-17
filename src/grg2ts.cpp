/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2026 April Wei
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
#include "grgl/grg2ts.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/mutation.h"
#include "tskit/core.h"
#include "tskit/tables.h"
#include "tskit/trees.h"
#include "tskit_util.h"
#include "util.h"

#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <tskit.h>

#define GRG2TS_VALIDATION 0

#define TSKIT_OK_OR_THROW(ok, msg)                                                                                     \
    do {                                                                                                               \
        int tskit_ok_val = (ok);                                                                                       \
        if (tskit_ok_val != 0) {                                                                                       \
            std::stringstream errMsg;                                                                                  \
            errMsg << (msg) << ": err=" << tskit_ok_val << ", " << tsk_strerror(tskit_ok_val);                         \
            throw TskitApiFailure(errMsg.str().c_str());                                                               \
        }                                                                                                              \
    } while (0)

#define TSKIT_ID_OR_THROW(ok, eid, msg)                                                                                \
    do {                                                                                                               \
        int tskit_ok_val = (ok);                                                                                       \
        if (tskit_ok_val < 0) {                                                                                        \
            std::stringstream errMsg;                                                                                  \
            errMsg << (msg) << ": err=" << tskit_ok_val << ", " << tsk_strerror(tskit_ok_val);                         \
            throw TskitApiFailure(msg);                                                                                \
        }                                                                                                              \
        if ((eid) != TSK_NULL && (tskit_ok_val != (eid))) {                                                            \
            std::stringstream errMsg;                                                                                  \
            errMsg << (msg) << ": expect=" << tskit_ok_val << ", got=" << (eid);                                       \
            throw TskitApiFailure(errMsg.str().c_str());                                                               \
        }                                                                                                              \
    } while (0)

namespace grgl {

class GrgToTsContext {
public:
    explicit GrgToTsContext(tsk_table_collection_t* tsTables, GRGPtr& grg)
        : m_tables(tsTables),
          m_grgNodeValid(grg->numNodes(), false),
          m_currentChildren(grg->numNodes(), 0),
          m_numSamples(grg->numSamples()) {
        TSKIT_OK_OR_THROW(tsk_node_table_init(&tsTables->nodes, 0), "Node table init");
        TSKIT_OK_OR_THROW(tsk_edge_table_init(&tsTables->edges, 0), "Edge table init");
        TSKIT_OK_OR_THROW(tsk_individual_table_init(&tsTables->individuals, 0), "Individual table init");
        TSKIT_OK_OR_THROW(tsk_site_table_init(&tsTables->sites, 0), "Site table init");
        TSKIT_OK_OR_THROW(tsk_mutation_table_init(&tsTables->mutations, 0), "Mutation table init");

        // For non-sample nodes, we use the nodeID as the time, since the GRG nodeIDs are in topological
        // order, guaranteeing that we have parent(time) > child(time).
        api_exc_check(grg->nodesAreOrdered(), "GRG->TS requires topologically ordered GRG nodes (e.g., ImmutableGRG)");

        // TODO:
        // 1. Populations
        // 2. Individuals
        // 3. When time is present on Mutations, add to the table.

        // Create all the sample nodes.
        for (m_nextTsId = 0; m_nextTsId < grg->numSamples(); m_nextTsId++) {
            TSKIT_ID_OR_THROW(
                tsk_node_table_add_row(&tsTables->nodes, TSK_NODE_IS_SAMPLE, 0.0, TSK_NULL, TSK_NULL, nullptr, 0),
                m_nextTsId,
                "Failed to add node");
            m_currentRoots.emplace(m_nextTsId);
        }
        // All GRG nodes have the same ID as their TS counterpart, so allocate them now.
        for (m_nextTsId = (tsk_id_t)grg->numSamples(); m_nextTsId < grg->numNodes(); m_nextTsId++) {
            TSKIT_ID_OR_THROW(
                tsk_node_table_add_row(&m_tables->nodes, 0, (double)m_nextTsId, TSK_NULL, TSK_NULL, nullptr, 0),
                m_nextTsId,
                "Failed to add node");
        }
    }

    // Return TSK_NULL if the GRG node is not in the current coalescent tree, otherwise return the
    // tsk_id_t (which will always be identifical to the grgNodeId)
    inline tsk_id_t getCurrentNode(const NodeID grgNodeId) const {
        // Samples exist in EVERY tree
        if (grgNodeId < m_numSamples) {
            return (tsk_id_t)grgNodeId;
        }
        return m_grgNodeValid.at(grgNodeId) ? (tsk_id_t)grgNodeId : TSK_NULL;
    }

    tsk_id_t getTreeParent(tsk_id_t tsNodeId) const {
        auto findIt = m_currentEdges.find(tsNodeId);
        if (findIt == m_currentEdges.end()) {
            return TSK_NULL;
        }
        return findIt->second.parent;
    }

    // Add to the table collection and delete from the current tree.
    void finalizeTreeEdge(const tsk_id_t tsParentId, const tsk_id_t tsChildId, const BpPosition endPos) {
        auto findIt = m_currentEdges.find(tsChildId);
        release_assert(findIt != m_currentEdges.end());
        const TSEdge& edge = findIt->second;
        release_assert(edge.parent == tsParentId);

        TSKIT_ID_OR_THROW(
            tsk_edge_table_add_row(&m_tables->edges, edge.start, (double)endPos, edge.parent, edge.child, nullptr, 0),
            TSK_NULL,
            "Failed to add edge");

        m_currentEdges.erase(findIt);

        NodeIDSizeT& parentsChildren = currentChildren(tsParentId);
        release_assert(parentsChildren > 0);
        parentsChildren--;
    }

    // Clear out the edges above a particular node -- this is only called after orphaning the node, which created
    // the root node, so no additional roots are created.
    void invalidateTreeAbove(const tsk_id_t tsNodeId, const BpPosition position) {
        // We start out by terminating edges, because the first edge is always terminated (it is the one
        // that has started the whole invalidation process).
        bool terminateEdges = true;

        auto findIt = m_currentEdges.find(tsNodeId);
        while (findIt != m_currentEdges.end()) {
            // Deleting the edge above the current node and emit that edge to the tskit table.
            const tsk_id_t parent = findIt->second.parent;
            const tsk_id_t child = findIt->second.child;
            release_assert(m_currentRoots.find(child) == m_currentRoots.end());

            // Move to next edge, so we can reason about it below.
            findIt = m_currentEdges.find(parent);

            // Terminate the edge, if requested.
            if (terminateEdges) {
                finalizeTreeEdge(parent, child, position);

                // If the child node has children, then it is now a root because we deleted its parent.
                // If it has no children, then it is an orphaned node, and is only a root if it is a sample.
                if (currentChildren(child) > 0 || child < m_numSamples) {
                    release_assert(m_currentRoots.emplace(child).second);
                }

                // If we are at the root of the path, and we are terminating edges, then this parent node
                // was a tree root, and it is no longer.
                if (findIt == m_currentEdges.end() && currentChildren(parent) == 0) {
                    release_assert(m_currentRoots.erase(parent) == 1);
                }
            }

            // Invalidate the node: samples beneath no longer match the GRG's samples beneath.
            if (parent < m_grgNodeValid.size()) {
                m_grgNodeValid[parent] = false;
            }

            // If the parent has other children, then we don't want to terminate the edges above the parent,
            // but we need to invalidate all the of the nodes on the path to the root.
            if (terminateEdges && currentChildren(parent) == 1) {
                terminateEdges = false;
            }
        }
    }

    void addTreeParent(tsk_id_t tsChildId, tsk_id_t tsParentId, BpPosition startPos) {
        // Edges are always added to start at the current tree. This can result in some slightly wonky stuff
        // where you get a path new1->new2-> <dangling> in tree i, and then in tree i+1 we finish the path
        // because we had to perform a deletion of an edge in the previous tree before we could add the new edge(s).
        // However, it makes the trees much easier to understand, and it makes it easier for us to add intervals
        // for synthetic edges that we create when adding roots to trees.
        release_assert(startPos >= m_currentTreeStart);
        const auto insertIt = m_currentEdges.emplace(tsChildId, TSEdge({tsParentId, tsChildId, m_currentTreeStart}));
        if (!insertIt.second) {
            release_assert(insertIt.first->second.parent == tsParentId);
            release_assert(m_currentRoots.find(tsChildId) == m_currentRoots.end());
        } else {
            m_currentRoots.erase(tsChildId);
            NodeIDSizeT& numChildren = currentChildren(tsParentId);
            numChildren++;
        }
    }

    tsk_id_t createTsNode(NodeID grgNodeId) {
        const tsk_id_t nodeId = (tsk_id_t)grgNodeId;
        m_grgNodeValid.at(grgNodeId) = true;
        // Add to roots if applicable.
        if (m_currentEdges.find(nodeId) == m_currentEdges.end()) {
            m_currentRoots.emplace(nodeId);
        }
        return nodeId;
    }

    void finalize(BpPosition position) {
        rootTheTree(position);
        for (auto& edgePair : m_currentEdges) {
            const TSEdge& edge = edgePair.second;
            TSKIT_ID_OR_THROW(tsk_edge_table_add_row(
                                  &m_tables->edges, edge.start, (double)position, edge.parent, edge.child, nullptr, 0),
                              TSK_NULL,
                              "Failed to add edge");
        }
        // Clear all the tree metadata so we can't accidentally use it again.
        m_currentEdges.clear();
        m_currentChildren.clear();
        m_currentRoots.clear();
    }

    // If there is more than one root, add a new node that is the root of the tree!
    void rootTheTree(BpPosition nextTreeStart) {
        // If we have multiple edge deletions due to the same mutation, then we don't need to re-root
        // the (previous) tree multiple times.
        if (nextTreeStart == m_currentTreeStart) {
            return;
        }
        if (m_currentRoots.size() > 1) {
            // Create a single new node, add an edge to it from each previous root.
            const tsk_id_t newRoot = m_nextTsId++;
            auto rootListCopy = m_currentRoots;
            for (tsk_id_t oldRoot : rootListCopy) {
                // Add an edge from the new root to the old root, with the same genomic start position as the
                // minimum start position of edges beneath the old root.
                addTreeParent(oldRoot, newRoot, m_currentTreeStart);
            }
            TSKIT_ID_OR_THROW(
                tsk_node_table_add_row(&m_tables->nodes, 0, (double)newRoot, TSK_NULL, TSK_NULL, nullptr, 0),
                newRoot,
                "Failed to add node");

            currentChildren(newRoot) = rootListCopy.size();
            m_currentRoots = {newRoot};
        }
        m_currentTreeStart = nextTreeStart;
    }

    BpPosition currentTreeStart() const { return m_currentTreeStart; }

    size_t numRoots() const { return m_currentRoots.size(); }

#if GRG2TS_VALIDATION
    bool validateRoots() const {
        std::unordered_set<tsk_id_t> realRoots;
        // Every root is either a sample node, or reachable from edges.
        for (NodeIDSizeT s = 0; s < m_numSamples; s++) {
            tsk_id_t p = s;
            auto it = m_currentEdges.find(p);
            while (it != m_currentEdges.end()) {
                p = it->second.parent;
                it = m_currentEdges.find(p);
            }
            // p is a root.
            realRoots.insert(p);
        }
        if (realRoots.size() != m_currentRoots.size()) {
            std::cout << "FAIL: " << realRoots.size() << " vs " << m_currentRoots.size() << "\n";
            for (auto r : realRoots) {
                if (m_currentRoots.find(r) == m_currentRoots.end()) {
                    std::cout << "  Only in COMPUTED roots: " << r << "\n";
                }
            }
            for (auto r : m_currentRoots) {
                if (realRoots.find(r) == realRoots.end()) {
                    std::cout << "  Only in TRACKED roots: " << r << "\n";
                }
            }
            return false;
        }
        return true;
    }
#endif

protected:
    NodeIDSizeT& currentChildren(tsk_id_t tsNode) {
        if (tsNode >= m_currentChildren.size()) {
            m_currentChildren.resize(tsNode + 1, 0);
        }
        return m_currentChildren[tsNode];
    }

    // The table collection representing our TreeSequence.
    tsk_table_collection_t* m_tables;

    // If true, the node n in GRG G is represented in the current coalescent tree T in the same way as
    // in G: i.e., samples below T{S(n)} = G{S(n)} is the same in both structures. When false, the node
    // is either not in T, or T{S(n)} is a subset of G{S(n)}.
    std::vector<bool> m_grgNodeValid;

    struct TSEdge {
        tsk_id_t parent;
        tsk_id_t child;
        BpPosition start;
    };

    // The edges in the current tree (mapped from child -> edge).
    std::unordered_map<tsk_id_t, TSEdge> m_currentEdges;
    // The set of root nodes in the current tree. Note: all other nodes in the tree
    // will have edges associated with them, only the roots have no (up) edge.
    std::unordered_set<tsk_id_t> m_currentRoots;
    // Map from tskit ID to the number of children beneath the node, _in the current tree_
    std::vector<NodeIDSizeT> m_currentChildren;
    // Start position (BP) of the current tree. Defined as the position of the most recent
    // edge deletion, which delineates the previous tree from the current tree.
    BpPosition m_currentTreeStart{};

    // Counter for tskit nodes
    tsk_id_t m_nextTsId{};

    // Number of haplotypes in our dataset.
    NodeIDSizeT m_numSamples;
};

/**
 * Given the context, which includes the current coalescent tree, add the given GRG node to the tree
 * by looking up if we already have a corresponding TS node (or creating one).
 */
static tsk_id_t addMutationToTree(GrgToTsContext& context,
                                  GRGPtr& grg,
                                  const NodeID grgNodeId,
                                  const BpPosition position,
                                  const tsk_id_t tsParentId = TSK_NULL) {
    tsk_id_t tsNodeId = context.getCurrentNode(grgNodeId);
    // If the node is new to this tree (or was invalidated for this tree), then we need to recursively
    // call this function. Otherwise, we can stop after processing this single node.
    const bool recurse = (tsNodeId == TSK_NULL);
    if (recurse) {
        tsNodeId = context.createTsNode(grgNodeId);
    }
    if (tsParentId != TSK_NULL) {
        // First, check for other parents. If we already have a parent, delete it and add the edge
        // to the table (terminating now).
        const tsk_id_t tsOtherParent = context.getTreeParent(tsNodeId);
        if (tsOtherParent != tsParentId) {
            if (tsOtherParent != TSK_NULL) {
                // First, make sure our current tree is properly rooted, because when we delete the
                // edge (next step) it might start a _NEW_ tree.
                context.rootTheTree(position);

                // We also need to invalidate the entire upward path in our tree, because it no longer reaches
                // the list of samples that it did, so cannot reuse any of those nodes.
                context.invalidateTreeAbove(tsNodeId, position);
            }

            // Next, add in our new parent to the tree, and associate the position with the start of the
            // new edge.
            context.addTreeParent(tsNodeId, tsParentId, position);
        }
    }
    if (recurse) {
        for (NodeID child : grg->getDownEdges(grgNodeId)) {
            addMutationToTree(context, grg, child, position, tsNodeId);
        }
    }
    return tsNodeId;
}

using MutAndTSNode = std::pair<Mutation, tsk_id_t>;

void convertGRGToTreeSeq(GRGPtr& grg, tsk_treeseq_t* outTS, std::pair<size_t, size_t> treeRange) {
    api_exc_check(grg->samplesAreOrdered(), "Samples must be ordered from 0...(N-1) for GRG->TS conversion");
    api_exc_check(!grg->hasMissingData(), "GRG has missing data; not supported for GRG->TS conversion");

    tsk_table_collection_t tsTables;
    TSKIT_OK_OR_THROW(tsk_table_collection_init(&tsTables, 0), "Failed allocated table collection");
    const BpPosition lastPosition = std::max(grg->getBPRange().second, grg->getSpecifiedBPRange().second);
    tsTables.sequence_length = (double)lastPosition;

    GrgToTsContext context(&tsTables, grg);

    tsk_id_t lastSiteId = TSK_NULL;
    grgl::BpPosition prevPos = INVALID_POSITION;
    for (auto& mutAndNode : grg->getMutationsToNodeOrdered()) {
        const Mutation& mut = grg->getMutationById(mutAndNode.first);
        api_exc_check(!mut.isMissing(), "GRG has missing data; not supported for GRG->TS conversion");
        const NodeID grgNode = mutAndNode.second;
        // Mutations in a GRG can be floating - no graph association. TS does not allow this.
        if (grgNode == INVALID_NODE_ID) {
            continue;
        }

        // Update the tree topology to reflect this mutation, and return the tskit node that is
        // immediately below the mutation.
        const tsk_id_t tsNode = addMutationToTree(context, grg, grgNode, mut.getPosition(), TSK_NULL);
#if GRG2TS_VALIDATION
        // This is very slow, so we only use it optionally when testing code changes.
        release_assert(context.validateRoots());
#endif

        // Update the sites and mutations tables with this mutation information, and its mapping to the
        // given node.
        if (mut.getPosition() != prevPos) {
            lastSiteId = tsk_site_table_add_row(&tsTables.sites,
                                                (double)mut.getPosition(),
                                                mut.getRefAllele().c_str(),
                                                mut.getRefAllele().size(),
                                                nullptr,
                                                0);
            TSKIT_ID_OR_THROW(lastSiteId, TSK_NULL, "Site add");
            prevPos = mut.getPosition();
        }
        TSKIT_ID_OR_THROW(tsk_mutation_table_add_row(&tsTables.mutations,
                                                     lastSiteId,
                                                     tsNode,
                                                     TSK_NULL,
                                                     TSK_UNKNOWN_TIME,
                                                     mut.getAllele().c_str(),
                                                     mut.getAllele().size(),
                                                     nullptr,
                                                     0),
                          TSK_NULL,
                          "Mutation add");
    }
    release_assert(prevPos == INVALID_POSITION || lastPosition >= prevPos);
    context.finalize(lastPosition);

    TSKIT_OK_OR_THROW(tsk_table_collection_sort(&tsTables, nullptr, 0), "Sort failed");
    TSKIT_OK_OR_THROW(tsk_table_collection_simplify(&tsTables, nullptr, 0, 0, nullptr), "Simplification failed");
    TSKIT_OK_OR_THROW(tsk_table_collection_build_index(&tsTables, 0), "Indexing failed");
    // TSK_TAKE_OWNERSHIP does not work here - we get a segfault during writing of the node table "metadata" column
    // (which should be empty). Instead we free the table collection manually after conversion.
    TSKIT_OK_OR_THROW(tsk_treeseq_init(outTS, &tsTables, 0), "Tree sequence creation failed");
    tsk_table_collection_free(&tsTables);
}

} // namespace grgl
