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

static bool s_debugging = false;
static size_t s_indent = 0;

std::string indentStr(size_t amount) { return std::string(amount, ' '); }

#define DEBUG_OUT(msg)                                                                                                 \
    do {                                                                                                               \
        if (s_debugging) {                                                                                             \
            std::cerr << indentStr(s_indent) << msg << "\n";                                                           \
        }                                                                                                              \
    } while (0)

namespace grgl {

class GrgToTsContext {
public:
    explicit GrgToTsContext(tsk_table_collection_t* tsTables, GRGPtr& grg)
        : m_tables(tsTables),
          m_grgNodeValid(grg->numNodes(), false),
          m_grgNodeUsed(grg->numNodes(), 0),
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
    // tsk_id_t (which will always be identical to the grgNodeId)
    inline tsk_id_t getAndUpdateCurrentNode(const NodeID grgNodeId, const BpPosition position) {
        // Samples exist in EVERY tree
        if (grgNodeId < m_numSamples) {
            return (tsk_id_t)grgNodeId;
        }
        if (!m_grgNodeValid.at(grgNodeId)) {
            return TSK_NULL;
        }
        m_grgNodeUsed[grgNodeId] = position; // Update that we have used this node at the given position.
        return (tsk_id_t)grgNodeId;
    }

    tsk_id_t getTreeParent(tsk_id_t tsNodeId) const {
        auto findIt = m_currentEdges.find(tsNodeId);
        if (findIt == m_currentEdges.end()) {
            return TSK_NULL;
        }
        return findIt->second.parent;
    }

    // Add to the table collection and delete from the current tree.
    void finalizeTreeEdge(const tsk_id_t tsParentId, const tsk_id_t tsChildId) {
        DEBUG_OUT("finalizeTreeEdge(" << tsParentId << ", " << tsChildId << ")");
        auto findIt = m_currentEdges.find(tsChildId);
        release_assert(findIt != m_currentEdges.end());
        const TSEdge& edge = findIt->second;
        release_assert(edge.parent == tsParentId);
        m_edgePendingDeletion.push_back(edge);
        m_currentEdges.erase(findIt);

        NodeIDSizeT& parentsChildren = currentChildren(tsParentId);
        release_assert(parentsChildren > 0);
        parentsChildren--;
    }

    // Clear out the edges above a particular node -- this is only called after orphaning the node, which created
    // the root node, so no additional roots are created.
    // Returns true if a new tree was started.
    void invalidateTreeAbove(const tsk_id_t tsNodeId, const BpPosition position) {
        // We start out by terminating edges, because the first edge is always terminated (it is the one
        // that has started the whole invalidation process).
        bool terminateEdges = true;

        std::vector<std::pair<tsk_id_t, tsk_id_t>> edgesToDelete;

        auto findIt = m_currentEdges.find(tsNodeId);
        while (findIt != m_currentEdges.end()) {
            // Deleting the edge above the current node and emit that edge to the tskit table.
            const tsk_id_t parent = findIt->second.parent;
            const tsk_id_t child = findIt->second.child;
            release_assert(m_currentRoots.find(child) == m_currentRoots.end());

            // Move to next edge, so we can reason about it below.
            findIt = m_currentEdges.find(parent);

            // Edges are [left, right) - so if we terminate an edge at right then the next tree
            // starts at [right, ...). grgNodeUsed[] reflects a value <right if that usage happened
            // within the edge [left, right)

            // XXX there is no guarantee that two mutations at the same position in a GRG will be
            // representable by a single tree. This is obviously problematic for the TS representation,
            // and we need to figure out a way to detect and represent it.
            if (parent < m_grgNodeUsed.size() && m_grgNodeUsed[parent] >= m_currentTreeStart) {
                m_currentTreeStart = m_grgNodeUsed[parent] + 1;
            }

            // Terminate the edge, if requested.
            if (terminateEdges) {
                edgesToDelete.emplace_back(parent, child);
            }

            // Invalidate the node: samples beneath no longer match the GRG's samples beneath.
            if (parent < m_grgNodeValid.size()) {
                m_grgNodeValid[parent] = false;
            }

            // If the parent has other children, then we don't want to terminate the edges above the parent,
            // but we need to invalidate all the of the nodes on the path to the root.
            if (terminateEdges && currentChildren(parent) > 1) {
                terminateEdges = false;
            }
        }

        // FIXME: do we still need this separate, now that we don't decide the right position of the edge
        // until later?
        for (const auto& edge : edgesToDelete) {
            release_assert(m_currentTreeStart >= 1);
            finalizeTreeEdge(edge.first, edge.second);

            // If the child node has children, then it is now a root because we deleted its parent.
            // If it has no children, then it is an orphaned node, and is only a root if it is a sample.
            if (currentChildren(edge.second) > 0 || edge.second < m_numSamples) {
                release_assert(m_currentRoots.emplace(edge.second).second);
            }

            // If we are at the root of the path, and we are terminating edges, then this parent node
            // was a tree root, and it is no longer.
            if (currentChildren(edge.first) == 0) {
                auto findIt = m_currentRoots.find(edge.first);
                if (findIt != m_currentRoots.end()) {
                    m_currentRoots.erase(findIt);
                }
            }
        }
    }

    void addTreeParent(tsk_id_t tsParentId, tsk_id_t tsChildId, BpPosition startPos = INVALID_POSITION) {
        const auto insertIt = m_currentEdges.emplace(tsChildId, TSEdge({tsParentId, tsChildId, startPos}));
        if (!insertIt.second) {
            api_exc_check(false, "TODO: checking whether existing edges are ever duplicated");
            release_assert(insertIt.first->second.parent == tsParentId);
            release_assert(m_currentRoots.find(tsChildId) == m_currentRoots.end());
        } else {
            const bool implicit = (startPos == INVALID_POSITION);
            DEBUG_OUT("addTreeParent(" << tsParentId << "-->" << tsChildId << " @ tree_pos="
                                       << (implicit ? m_currentTreeStart : startPos) << (implicit ? "i" : "") << ")");
            // By default, edges are added eagerly to the current tree, but we won't know where they start until
            // we complete the previous tree, so we defer the setting of the start position.
            if (startPos == INVALID_POSITION) {
                m_pendingChildren.emplace_back(tsChildId);
            }
            NodeIDSizeT& numChildren = currentChildren(tsParentId);
            numChildren++;
        }
    }

    tsk_id_t createTsNode(const NodeID grgNodeId, const BpPosition position) {
        const tsk_id_t nodeId = (tsk_id_t)grgNodeId;
        m_grgNodeValid.at(grgNodeId) = true;
        m_grgNodeUsed.at(grgNodeId) = position;
        return nodeId;
    }

    // Save the current roots in case we have a breaking change.
    void saveRoots() { m_previousRoots = m_currentRoots; }

    // If needed, add this node to the roots.
    void checkAddRoot(const tsk_id_t nodeId) {
        // Add to roots if applicable.
        if (m_currentEdges.find(nodeId) == m_currentEdges.end()) {
            DEBUG_OUT("Adding root " << nodeId);
            m_currentRoots.emplace(nodeId);
        }
    }

    void removeRoot(const tsk_id_t nodeId) {
        DEBUG_OUT("Removing root " << nodeId);
        m_currentRoots.erase(nodeId);
    }

    void finalize(BpPosition position) {
        DEBUG_OUT("FINALIZE: roots are: ");
        for (auto root : m_currentRoots) {
            DEBUG_OUT("   root=" << root);
        }
        m_currentTreeStart = position;
        m_previousRoots = m_currentRoots;
        m_currentRoots.clear();
        flushPendingEdges(/*force=*/true);
        for (auto& edgePair : m_currentEdges) {
            const TSEdge& edge = edgePair.second;
            if (edge.start < position) {
                DEBUG_OUT("Finalizing edge " << edge.parent << "-->" << edge.child << " (" << edge.start << "-"
                                             << position << ")");
                TSKIT_ID_OR_THROW(
                    tsk_edge_table_add_row(
                        &m_tables->edges, edge.start, (double)position, edge.parent, edge.child, nullptr, 0),
                    TSK_NULL,
                    "Failed to add edge");
            }
        }
        // Clear all the tree metadata so we can't accidentally use it again.
        m_currentEdges.clear();
        m_currentChildren.clear();
        m_currentRoots.clear();
    }

    // If there is more than one root, add a new node that is the root of the tree. This is kind of
    // tricky to do as we move from left-to-right, since nodes that are roots in the previous tree
    // may be children in the current tree. We handle both cases:
    // 1. If the root is a root in both trees, add the edge and do not terminate it.
    // 2. If it is only a root in the previous tree, just add the edge spanning (prevStart - curStart)
    void rootTheTree() {
        if (m_previousRoots.size() > 1) {
            DEBUG_OUT("Rooting the tree that starts at " << m_previousTreeStart << " and ends at "
                                                         << m_currentTreeStart);
            // Create a single new node, add an edge to it from each previous root.
            const tsk_id_t newRoot = m_nextTsId++;
            for (tsk_id_t oldRoot : m_previousRoots) {
                // If this root is _not_ a root in the current tree, then we just add the edge to
                // the TreeSequence, not to the current tree.
                if (m_currentRoots.find(oldRoot) == m_currentRoots.end()) {
                    DEBUG_OUT("Root edge: " << newRoot << "-->" << oldRoot << " (" << m_previousTreeStart << "-"
                                            << m_currentTreeStart << ")");
                    TSKIT_ID_OR_THROW(tsk_edge_table_add_row(&m_tables->edges,
                                                             (double)m_previousTreeStart,
                                                             (double)m_currentTreeStart,
                                                             newRoot,
                                                             oldRoot,
                                                             nullptr,
                                                             0),
                                      TSK_NULL,
                                      "Failed to add edge");
                } else {
                    DEBUG_OUT("Continuing root edge: " << newRoot << "-->" << oldRoot << " (" << m_previousTreeStart
                                                       << "-?)");
                    addTreeParent(newRoot, oldRoot, m_previousTreeStart);
                    removeRoot(oldRoot);
                    m_currentRoots.emplace(newRoot);
                }
            }
            TSKIT_ID_OR_THROW(
                tsk_node_table_add_row(&m_tables->nodes, 0, (double)newRoot, TSK_NULL, TSK_NULL, nullptr, 0),
                newRoot,
                "Failed to add node");

            m_previousRoots = {newRoot};
        }
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

    // Flushing the edges produces a new tree, which started at m_previousTreeStart and ended
    // at m_currentTreeStart.
    void flushPendingEdges(bool force = false) {
        if (m_previousTreeStart != m_currentTreeStart || force) {
            rootTheTree();
            m_previousTreeStart = m_currentTreeStart;
        }

        DEBUG_OUT("Flushing deleted edges:");
        for (const auto& edge : m_edgePendingDeletion) {
            release_assert(edge.start <= m_currentTreeStart);
            if (edge.start < m_currentTreeStart) {
                DEBUG_OUT(">" << edge.parent << "-->" << edge.child << " (" << edge.start << "-" << m_currentTreeStart
                              << ")");
                TSKIT_ID_OR_THROW(
                    tsk_edge_table_add_row(
                        &m_tables->edges, edge.start, (double)m_currentTreeStart, edge.parent, edge.child, nullptr, 0),
                    TSK_NULL,
                    "Failed to add edge");
            }
        }
        m_edgePendingDeletion.clear();

        DEBUG_OUT("Flushing added edges:");
        for (const auto& childId : m_pendingChildren) {
            auto& edge = m_currentEdges.at(childId);
            edge.start = m_currentTreeStart;
            DEBUG_OUT(">" << edge.parent << "-->" << edge.child << " starts at " << m_currentTreeStart);
        }
        m_pendingChildren.clear();
    }

protected:
    NodeIDSizeT& currentChildren(tsk_id_t tsNode) {
        if (tsNode >= m_currentChildren.size()) {
            m_currentChildren.resize(tsNode + 1, 0);
        }
        return m_currentChildren[tsNode];
    }

    // The table collection representing our TreeSequence.
    tsk_table_collection_t* m_tables;

    // At what position (left-to-right) was the node last used. "Used" means that the meaning of the GRG
    // node was relied upon w.r.t. the samples beneath it. When edges are deleted from the current tree,
    // the node below which an edge was deleted no longer has the same meaning (samples-beneath changes),
    // and becomes invalid.
    // 1. Only GRG-derived nodes can be "used", by this definition. Synthetic nodes cannot be.
    // 2. A value of INVALID_POSITION means that the node is invalid. The next time a copy of a mutation
    //    from the GRG needs to use this node, the edges will have to be reconstructed.
    std::vector<bool> m_grgNodeValid;
    std::vector<BpPosition> m_grgNodeUsed;

    struct TSEdge {
        tsk_id_t parent;
        tsk_id_t child;
        BpPosition start;
    };

    std::vector<TSEdge> m_edgePendingDeletion;
    std::vector<tsk_id_t> m_pendingChildren;

    // The edges in the current tree (mapped from child -> edge).
    std::unordered_map<tsk_id_t, TSEdge> m_currentEdges;
    // The set of root nodes in the current tree. Note: all other nodes in the tree
    // will have edges associated with them, only the roots have no (up) edge.
    std::unordered_set<tsk_id_t> m_currentRoots;
    std::unordered_set<tsk_id_t> m_previousRoots;
    // Map from tskit ID to the number of children beneath the node, _in the current tree_
    std::vector<NodeIDSizeT> m_currentChildren;
    // Start position (BP) of the current tree. Defined as the position of the most recent
    // edge deletion, which delineates the previous tree from the current tree.
    BpPosition m_currentTreeStart{0};
    BpPosition m_previousTreeStart{0};

    // Counter for tskit nodes
    tsk_id_t m_nextTsId{0};

    // Number of haplotypes in our dataset.
    NodeIDSizeT m_numSamples;
};

// This algorithm is a bit hard (for me, at least) to get your head around.
// Consider the current tree that we are building, while we are processing mutation M at position P.
// After processing M, we do not necessarily have the tree for locus P - we are only gauranteed that
// the subtree below M is exactly how it will be. When processing subsequent mutations (e.g. M2 @ P2)
// we may add or remove edges to the tree that we "backdate" to cover locus P in addition to P2.
//
// When we actually flush the pending deletions and additions, they will be at the same boundary between
// trees:
// 1. After applying those deletions and additions, the TreeSequence tables will match the "current tree"
//    that we have in m_currentEdges (and which we do not yet know the end-range of)
// 2. Therefore, we never really know which "current tree" (m_currentEdges set) will cover the range
//     (m_previousTreeStart, m_currentTreeStart), unless we happened to save a copy of that tree.
//
// To root the tree covering (m_previousTreeStart, m_currentTreeStart), we then need to save a copy of
// the current tree's roots periodically, whenever we get a new "high water mark" for where our deletions
// and additions will be applied.

/**
 * Given the context, which includes the current coalescent tree, add the given GRG node to the tree
 * by looking up if we already have a corresponding TS node (or creating one).
 */
static tsk_id_t addHierarchyToTree(GrgToTsContext& context,
                                   GRGPtr& grg,
                                   const NodeID grgNodeId,
                                   const BpPosition position,
                                   const tsk_id_t tsParentId = TSK_NULL,
                                   size_t indent = 0) {
    s_indent = indent;
    DEBUG_OUT("addHierarchyToTree(" << grgNodeId << ", parent=" << tsParentId << ")");
    tsk_id_t tsNodeId = context.getAndUpdateCurrentNode(grgNodeId, position);
    // If the node is new to this tree (or was invalidated for this tree), then we need to recursively
    // call this function. Otherwise, we can stop after processing this single node.
    const bool recurse = (tsNodeId == TSK_NULL);
    if (recurse) {
        tsNodeId = context.createTsNode(grgNodeId, position);
        DEBUG_OUT("new node = " << tsNodeId);
    }
    if (tsParentId != TSK_NULL) {
        // First, check for other parents. If we already have a parent, delete it and add the edge
        // to the table (terminating now).
        const tsk_id_t tsOtherParent = context.getTreeParent(tsNodeId);
        DEBUG_OUT("new parent = " << tsParentId << ", old parent = " << tsOtherParent);
        if (tsOtherParent != tsParentId) {
            if (tsOtherParent != TSK_NULL) {
                // We need to invalidate the entire upward path in our tree, because it no longer reaches
                // the list of samples that it did, so cannot reuse any of those nodes.
                context.invalidateTreeAbove(tsNodeId, position);
            }

            // Next, add in our new parent to the current tree.
            context.addTreeParent(tsParentId, tsNodeId);
            context.removeRoot(tsNodeId);
        }
    }
    if (recurse) {
        for (NodeID child : grg->getDownEdges(grgNodeId)) {
            addHierarchyToTree(context, grg, child, position, tsNodeId, indent + 2);
        }
    }

    return tsNodeId;
}

/**
 * Given the context, which includes the current coalescent tree, add the given GRG node to the tree
 * by looking up if we already have a corresponding TS node (or creating one). This algorithm has the
 * following properties:
 * 1. Adding a single mutation can both add and remove edges from the tree. All added edges will have
 *    the same start position (the start of the tree) and all removed edges will have the same
 *    end position.
 */
static tsk_id_t
addMutationToTree(GrgToTsContext& context, GRGPtr& grg, const NodeID grgNodeId, const BpPosition position) {
    DEBUG_OUT("\naddMutationToTree(" << grgNodeId << ", position=" << position << ") {");

    context.saveRoots();

    // This constructs a subtree rooted at (the tskit equivalent of) grgNodeId by following
    // all down edges in the GRG, by modifying the current marginal tree in the TreeSequence.
    const tsk_id_t mutNode = addHierarchyToTree(context, grg, grgNodeId, position);
    context.checkAddRoot(mutNode);
    context.flushPendingEdges();

    DEBUG_OUT("}\n");
    return mutNode;
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

    MutationId debugMut = INVALID_MUTATION_ID;
    if (std::getenv("DEBUG_MUT") != nullptr) {
        debugMut = std::atoi(std::getenv("DEBUG_MUT"));
        std::cerr << "Debugging MutationId=" << debugMut << "\n";
    }

    MutationId stopAtMut = INVALID_MUTATION_ID;
    if (std::getenv("STOP_AT_MUT") != nullptr) {
        stopAtMut = std::atoi(std::getenv("STOP_AT_MUT"));
        std::cerr << "Stopping at MutationId=" << stopAtMut << "\n";
    }

    tsk_id_t lastSiteId = TSK_NULL;
    grgl::BpPosition prevPos = INVALID_POSITION;
    for (auto& mutAndNode : grg->getMutationsToNodeOrdered()) {
        if (mutAndNode.first == stopAtMut) {
            break;
        }
        if (mutAndNode.first == debugMut) {
            s_debugging = true;
        } else {
            s_debugging = false;
        }
        const Mutation& mut = grg->getMutationById(mutAndNode.first);
        api_exc_check(!mut.isMissing(), "GRG has missing data; not supported for GRG->TS conversion");
        const NodeID grgNode = mutAndNode.second;
        // Mutations in a GRG can be floating - no graph association. TS does not allow this.
        if (grgNode == INVALID_NODE_ID) {
            continue;
        }

        // Update the tree topology to reflect this mutation, and return the tskit node that is
        // immediately below the mutation.
        const tsk_id_t tsNode = addMutationToTree(context, grg, grgNode, mut.getPosition());
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
