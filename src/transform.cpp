#include "grgl/transform.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/node_data.h"
#include "grgl/visitor.h"
#include "util.h"

#include <memory>
#include <unordered_map>

namespace grgl {

// The general idea is to find sibling node groups that share a lot of common children, and create new
// nodes beneath these siblings to capture this sharing. Doing this optimally is probably somewhat
// intractable, so we: iterate the nodes in topological order, look at all siblings to our current node,
// order them by the one that shares the most with us, construct hierarchy, and continue. We do this
// even if the sharing is identical, because the simplification step (when serializing the graph) will
// clean it up anyways.
//
// This algorithm does not RELY on topological order, but the results are probably better for it. As
// we make more modifications to the graph we are breaking the order (so node order != topo order).
size_t reduceGRG(const MutableGRGPtr& mutGRG) {
    // If we have 2 or fewer shared children it doesn't benefit us to create hierarchy.
    constexpr size_t MIN_SHARED = 2;
    // If we remove 3 edges, we are adding a single node and 2 edges. This ratio is meant to
    // very roughly capture the storage/compute cost of a node vs. an edge.
    constexpr size_t NODE_COST = 3;
    // This is just to help us prioritize work in an expensive optimization problem that we have
    // made very greedy. We modify hierarchy below nodes where at least 60% of their children are shared.
    constexpr double CHILD_PROPORTION = 0.6;

    const NodeIDSizeT origNodes = mutGRG->numNodes();
    size_t removedEdges = 0;
    std::vector<bool> covered(origNodes, false);
    for (NodeID node = mutGRG->numSamples(); node < origNodes; node++) {
        if (covered[node]) {
            continue;
        }
        if (mutGRG->numDownEdges(node) < (MIN_SHARED * NODE_COST)) {
            continue;
        }
        // Go down first
        const auto children = mutGRG->getDownEdges(node);
        const NodeIDSizeT target = children.size();
        // For all our siblings, find the one with the most sharing.
        NodeID bestSibling = INVALID_NODE_ID;
        size_t bestCount = 0;
        std::unordered_map<NodeID, NodeIDList> shared;
        for (const NodeID child : children) {
            for (const NodeID parent : mutGRG->getUpEdges(child)) {
                if (parent >= origNodes || parent == node) {
                    continue;
                }
                auto insertIt = shared.emplace(parent, NodeIDList());
                auto& sharedList = insertIt.first->second;
                sharedList.push_back(child);
                if (sharedList.size() > bestCount) {
                    bestCount = sharedList.size();
                    bestSibling = parent;
                }
            }
        }
        // Nodes have a size cost that is larger than edges, because of NodeData, and because we
        // use 64-bits per node in the CSR representation. So we account for that cost here.
        if (((double)bestCount / (double)target) >= CHILD_PROPORTION && bestCount >= (MIN_SHARED * NODE_COST)) {
            release_assert(bestSibling != INVALID_NODE_ID);
            covered[bestSibling] = true;
            const auto& sharedList = shared.at(bestSibling);
            if (bestCount < target) {
                const NodeID newNode = mutGRG->makeNode();
                for (const NodeID child : sharedList) {
                    release_assert(mutGRG->disconnect(bestSibling, child));
                    release_assert(mutGRG->disconnect(node, child));
                    mutGRG->connect(newNode, child);
                }
                mutGRG->connect(node, newNode);
                mutGRG->connect(bestSibling, newNode);

                // These nodes are either new, or their children have been modified. As such as need to mark
                // them for getting the coalescence counts recalculated later.
                mutGRG->setNumIndividualCoalsGrow(newNode, COAL_COUNT_NOT_SET);
                mutGRG->setNumIndividualCoals(node, COAL_COUNT_NOT_SET);
                mutGRG->setNumIndividualCoals(bestSibling, COAL_COUNT_NOT_SET);

                removedEdges += (sharedList.size() - 2);
            } else {
                release_assert(bestCount == target); // They are equal
                // There are two scenarios:
                // A) Nodes are identical.
                // B) "bestSibling" is a superset of "node"
                // We can handle both cases the same way by placing node beneath bestSibling, in place
                // of all those shared edges.
                for (const NodeID child : sharedList) {
                    release_assert(mutGRG->disconnect(bestSibling, child));
                }
                mutGRG->connect(bestSibling, node);
                // node's coalescence is unchanged, but bestSibling no longer coalesces the children that
                // are beneath node now.
                const NodeIDSizeT nodeCoals = mutGRG->getNumIndividualCoals(node);
                const NodeIDSizeT existingCoals = mutGRG->getNumIndividualCoals(bestSibling);
                if (nodeCoals == COAL_COUNT_NOT_SET || existingCoals == COAL_COUNT_NOT_SET) {
                    mutGRG->setNumIndividualCoals(bestSibling, COAL_COUNT_NOT_SET);
                } else {
                    release_assert(existingCoals >= nodeCoals);
                    mutGRG->setNumIndividualCoals(bestSibling, existingCoals - nodeCoals);
                }
                // Since we didn't create a node, we need to tell the GRG explicitly that we broke node ordering.
                mutGRG->nodesAreOrdered() = false;
                removedEdges += sharedList.size();
            }
        }
    }
    return removedEdges;
}

size_t reduceGRGUntil(const grgl::MutableGRGPtr& mutGRG,
                      const size_t iterations,
                      const size_t minDropped,
                      const double fractionDropped,
                      const bool verbose) {
    const size_t origEdges = mutGRG->numEdges();
    size_t totalRemoved = 0;
    size_t lastRemoved = 0;
    size_t it = 0;
    do {
        lastRemoved = grgl::reduceGRG(mutGRG);
        totalRemoved += lastRemoved;
        if (verbose) {
            std::cout << "REDUCE: Removed " << lastRemoved << " edges in iteration " << it << "\n";
        }
        it++;
    } while ((it < iterations) && (totalRemoved < ((size_t)(fractionDropped * (double)origEdges))) &&
             (lastRemoved > minDropped));

    // Calculate the coalescence counts for new nodes and nodes that we modified, above.
    calculateMissingCoals(mutGRG);

    return it;
}

/**
 * Visitor that marks nodes which need to propagate individual information to their parents.
 */
class GRGCoalMarkVisitor : public GRGVisitor {
public:
    GRGCoalMarkVisitor() = default;

    bool visit(const GRGPtr& grg, NodeID node, TraversalDirection direction, DfsPass dfsPass = DFS_PASS_NONE) override {
        if (m_marked.empty()) {
            release_assert(direction == TraversalDirection::DIRECTION_UP);
            m_marked.resize(grg->numNodes(), false);
        }
        if (dfsPass == DFS_PASS_BACK_AGAIN) {
            bool mark = (grg->getNumIndividualCoals(node) == COAL_COUNT_NOT_SET);
            if (!mark) {
                for (const NodeID parent : grg->getUpEdges(node)) {
                    if (m_marked[parent]) {
                        mark = true;
                        break;
                    }
                }
            }
            m_marked[node] = mark;
        }
        return true;
    }

    std::vector<bool> m_marked;
};

/**
 * Visitor that recalculates coalescence information for particular nodes.
 */
class GRGCoalCalcVisitor : public GRGVisitor {
public:
    explicit GRGCoalCalcVisitor(const std::vector<bool>& marked)
        : m_marked(marked) {}

    bool visit(const GRGPtr& grg, NodeID node, TraversalDirection direction, DfsPass dfsPass = DFS_PASS_NONE) override {
        if (m_nodeToIndivs.uninitialized()) {
            release_assert(direction == TraversalDirection::DIRECTION_DOWN);
            m_nodeToIndivs = RefCountedNodeData<NodeIDList>(grg->numNodes(), flipDir(direction));
        }
        if (dfsPass == DFS_PASS_BACK_AGAIN) {
            NodeIDListUPtr uncoalesced = NodeIDListUPtr(new NodeIDList());
            const NodeIDSizeT coals = getCoalsForParent(
                grg, m_nodeToIndivs, node, grg->getDownEdges(node), *uncoalesced, /*implicitSamples=*/true);
            if (!uncoalesced->empty()) {
                m_nodeToIndivs.add(grg, node, std::move(uncoalesced));
            }
            const NodeIDSizeT existingCoals = grg->getNumIndividualCoals(node);
            if (existingCoals == COAL_COUNT_NOT_SET) {
                grg->setNumIndividualCoalsGrow(node, coals);
                m_updated++;
            } else {
                release_assert(existingCoals == coals);
            }
            return m_marked[node];
        }
        return true;
    }

    size_t m_updated{};

private:
    const std::vector<bool>& m_marked;
    RefCountedNodeData<NodeIDList> m_nodeToIndivs;
};

size_t calculateMissingCoals(const GRGPtr& grg) {
    if (grg->getPloidy() == 2) {
        // Mark the relevant nodes.
        GRGCoalMarkVisitor markVisitor;
        grg->visitDfs(markVisitor, TraversalDirection::DIRECTION_UP, grg->getSampleNodes());

        // Traverse them and recalculate coalescence information.
        GRGCoalCalcVisitor visitor(markVisitor.m_marked);
        fastCompleteDFS(grg, visitor);
        return visitor.m_updated;
    }
    return 0;
}

} // namespace grgl
