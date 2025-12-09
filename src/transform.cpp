#include "grgl/transform.h"
#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "util.h"

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
        std::unordered_map<NodeID, NodeIDList> shared;
        if (mutGRG->numDownEdges(node) < (MIN_SHARED * NODE_COST)) {
            continue;
        }
        // Go down first
        const auto children = mutGRG->getDownEdges(node);
        const NodeIDSizeT target = children.size();
        if (target <= MIN_SHARED) {
            continue;
        }
        // For all our siblings, find the one with the most sharing.
        NodeID bestSibling = INVALID_NODE_ID;
        size_t bestCount = 0;
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
        // We don't want to trade off nodes for edges, so we have bestCount >= 6
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
    return it;
}

} // namespace grgl
