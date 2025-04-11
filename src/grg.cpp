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
#include "grgl/grg.h"
#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/visitor.h"
#include "util.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <unordered_map>

namespace grgl {

// This number needs to be more principled in how it is chosen. I somewhat randomly chose this,
// as I want the pairwise case to definitely use the sparse topological search, but I don't really
// know when the switch from dense to sparse is optimal. This only matters for large graphs.
static constexpr size_t TOPO_DENSE_THRESHOLD = 10;

using bitvect = std::vector<bool>;

MutationId GRG::addMutation(const Mutation& mutation, const NodeID nodeId) {
    const MutationId mutId = m_mutations.size();
    m_mutations.push_back(mutation);
    if (nodeId != INVALID_NODE_ID) {
        release_assert(nodeId < this->numNodes());
    }
    this->m_mutIdsByNodeId.emplace_back(nodeId, mutId);
    this->m_mutIdsByNodeIdSorted = false;
    if (m_mutsAreOrdered) {
        m_mutsAreOrdered = false;
    }
    return mutId;
}

using MutIdAndNode = std::pair<MutationId, NodeID>;

struct MutIdAndNodeLt {
    explicit MutIdAndNodeLt(GRGPtr& theGRG)
        : grg(theGRG) {}

    // Sort order: position, ref allele, alt allele.
    bool operator()(const MutIdAndNode& lhs, const MutIdAndNode& rhs) {
        const auto& mlhs = grg->getMutationById(lhs.first);
        const auto& mrhs = grg->getMutationById(rhs.first);
        if (mlhs.getPosition() == mrhs.getPosition()) {
            const auto leftAlleles = mlhs.getBothAlleles();
            const auto rightAlleles = mrhs.getBothAlleles();
            if (leftAlleles.first == rightAlleles.first) {
                return leftAlleles.second < rightAlleles.second;
            }
            return leftAlleles.first < rightAlleles.first;
        }
        return mlhs.getPosition() < mrhs.getPosition();
    }

    GRGPtr& grg;
};

struct MutIdAndNodeLtFast {
    bool operator()(const MutIdAndNode& lhs, const MutIdAndNode& rhs) const { return lhs.first < rhs.first; }
};

std::vector<MutIdAndNode> GRG::getMutationsToNodeOrdered() {
    std::vector<MutIdAndNode> result;
    for (const auto& nodeAndMutId : m_mutIdsByNodeId) {
        result.emplace_back(nodeAndMutId.second, nodeAndMutId.first);
    }
    GRGPtr sharedThis = shared_from_this();
    if (mutationsAreOrdered()) {
        std::sort(result.begin(), result.end(), MutIdAndNodeLtFast());
    } else {
        std::sort(result.begin(), result.end(), MutIdAndNodeLt(sharedThis));
    }
    return std::move(result);
}

void GRG::visitBfs(GRGVisitor& visitor,
                   TraversalDirection direction,
                   const NodeIDList& seedList,
                   ssize_t maxQueueWidth) {
    GRGPtr sharedThis = shared_from_this();
    std::deque<NodeID> fifo;
    for (const auto& nodeId : seedList) {
        fifo.push_back(nodeId);
    }
    while (!fifo.empty()) {
        auto nodeId = fifo.front();
        fifo.pop_front();
        bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_NONE);
        if (keepGoing) {
            auto successors = (direction == DIRECTION_UP) ? this->getUpEdges(nodeId) : this->getDownEdges(nodeId);
            for (const auto& succId : successors) {
                if ((maxQueueWidth < 0) || fifo.size() < maxQueueWidth) {
                    fifo.push_back(succId);
                }
            }
        }
    }
}

void GRG::visitDfs(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList, bool forwardOnly) {
    GRGPtr sharedThis = shared_from_this();
    const NodeMark is2ndPass = NODE_MARK_1;
    // Most STL implementations implement this as a packed bitvector.
    std::unique_ptr<bitvect> alreadySeen(new bitvect(this->numNodes()));
    std::list<NodeID> lifo;
    for (const auto& nodeId : seedList) {
        assert(!alreadySeen->at(nodeId));
        lifo.push_back(nodeId);
    }
    while (!lifo.empty()) {
        const auto markedNodeId = lifo.back();
        const auto nodeId = removeMarks(markedNodeId);
        assert(nodeId < this->numNodes());
        if (!hasMark(markedNodeId, is2ndPass)) {
            // First (forward) pass.
            if (alreadySeen->at(nodeId)) {
                lifo.pop_back();
                continue;
            }
            if (forwardOnly) {
                lifo.pop_back();
            } else {
                lifo.back() = markNodeId(nodeId, is2ndPass, true);
            }
            bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_THERE);
            if (keepGoing) {
                auto successors = (direction == DIRECTION_UP) ? this->getUpEdges(nodeId) : this->getDownEdges(nodeId);
                for (const auto& succId : successors) {
                    lifo.push_back(succId);
                }
            }
        } else {
            // Second (back up) pass.
            lifo.pop_back();
            visitor.visit(sharedThis, nodeId, direction, DFS_PASS_BACK_AGAIN);
            alreadySeen->at(nodeId) = true;
        }
    }
}

void MutableGRG::connect(const NodeID srcId, const NodeID tgtId) {
    m_nodes.at(srcId)->addDownEdge(tgtId);
    m_nodes.at(tgtId)->addUpEdge(srcId);
}

void MutableGRG::disconnect(const NodeID srcId, const NodeID tgtId) {
    m_nodes.at(srcId)->deleteDownEdge(tgtId);
    m_nodes.at(tgtId)->deleteUpEdge(srcId);
}

class TopoOrderVisitor : public grgl::GRGVisitor {
public:
    TopoOrderVisitor() = default;

    bool visit(const GRGPtr& grg,
               const NodeID nodeId,
               const TraversalDirection direction,
               const DfsPass dfsPass = grgl::DfsPass::DFS_PASS_NONE) override {
        (void)direction;
        if (m_result.empty()) {
            m_result.resize(grg->numNodes());
        }
        if (dfsPass == grgl::DfsPass::DFS_PASS_BACK_AGAIN) {
            this->m_result[nodeId] = m_counter++;
        }
        return true;
    }

    std::vector<NodeIDSizeT> m_result;

private:
    NodeIDSizeT m_counter{};
};

std::vector<NodeIDSizeT> MutableGRG::topologicalSort(TraversalDirection direction) {
    TopoOrderVisitor visitor;
    NodeIDList seeds =
        (direction == TraversalDirection::DIRECTION_DOWN) ? this->getRootNodes() : this->getSampleNodes();
    this->visitDfs(visitor, direction, seeds);
    return visitor.m_result;
}

class HeapNode {
public:
    HeapNode(const NodeID nodeId, const NodeIDSizeT order)
        : m_nodeId(nodeId),
          m_order(order) {}

    NodeID m_nodeId;
    NodeIDSizeT m_order;
};

// We invert the comparison so we can get a min-heap instead of the default max-heap
static bool cmpHeapNode(const HeapNode& hnode1, const HeapNode& hnode2) { return hnode1.m_order > hnode2.m_order; }

void MutableGRG::visitTopo(GRGVisitor& visitor,
                           TraversalDirection direction,
                           const NodeIDList& seedList,
                           const std::vector<NodeIDSizeT>* sortOrder) {
    // Short-cut for much faster access, when nodes are ordered (they usually are).
    if (this->nodesAreOrdered()) {
        if (sortOrder != nullptr) {
            throw ApiMisuseFailure("This GRG has its node topologically ordered; sortOrder unneeded");
        }
        if (seedList.size() > TOPO_DENSE_THRESHOLD) {
            this->visitTopoNodeOrderedDense(visitor, direction, seedList);
        } else {
            this->visitTopoNodeOrderedSparse(visitor, direction, seedList);
        }
        return;
    }

    GRGPtr sharedThis = shared_from_this();
    std::vector<NodeIDSizeT> _localSortOrder;
    if (nullptr == sortOrder) {
        _localSortOrder = topologicalSort((direction == DIRECTION_UP) ? DIRECTION_DOWN : DIRECTION_UP);
        sortOrder = &_localSortOrder;
    }
    NodeIDSet alreadySeen;
    std::vector<HeapNode> heap;
    heap.reserve(seedList.size() * 2);
    for (const auto& seedId : seedList) {
        heap.emplace_back(seedId, (*sortOrder).at(seedId));
    }
    NodeIDSizeT lastPopped = 0;
    std::make_heap(heap.begin(), heap.end(), cmpHeapNode);
    while (!heap.empty()) {
        std::pop_heap(heap.begin(), heap.end(), cmpHeapNode);
        const NodeID nodeId = heap.back().m_nodeId;
        release_assert(lastPopped <= heap.back().m_order);
        lastPopped = heap.back().m_order;
        heap.pop_back();

        if (alreadySeen.find(nodeId) != alreadySeen.end()) {
            continue;
        }
        alreadySeen.insert(nodeId);

        const bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_NONE);
        if (keepGoing) {
            auto successors = (direction == DIRECTION_UP) ? this->getUpEdges(nodeId) : this->getDownEdges(nodeId);
            for (const auto& succId : successors) {
                heap.emplace_back(succId, (*sortOrder).at(succId));
                std::push_heap(heap.begin(), heap.end(), cmpHeapNode);
            }
        }
    }
}

std::vector<NodeIDSizeT> CSRGRG::topologicalSort(TraversalDirection direction) {
    std::vector<NodeIDSizeT> result;
    if (direction == TraversalDirection::DIRECTION_DOWN) {
        for (size_t i = 0; i < this->numNodes(); i++) {
            result.emplace_back(i);
        }
    } else {
        for (size_t i = this->numNodes(); i > 0; i--) {
            result.emplace_back(i - 1);
        }
    }
    return result;
}

// MinHeap for the "up" direction (smaller nodeIDs are below larger ones)
static bool cmpNodeIdGt(const NodeID& node1, const NodeID& node2) { return node1 > node2; }
// MaxHeap for the "down" direction (larger nodeIDs are above smaller ones)
static bool cmpNodeIdLt(const NodeID& node1, const NodeID& node2) { return node1 < node2; }

// Visit the nodes in topological order, assuming that the nodeIDs are already in topological order.
// This function takes a set of seeds though, so only a subset of nodes may be visited. Using this
// linear bit-vector approach is actually significantly faster than using a heap of nodes, which is
// what we did previously.
void GRG::visitTopoNodeOrderedDense(GRGVisitor& visitor,
                                    const TraversalDirection direction,
                                    const NodeIDList& seedList) {
    GRGPtr sharedThis = shared_from_this();
    static std::vector<bool> yesVisit;
    yesVisit.clear();
    yesVisit.resize(numNodes());
    for (const auto& seedId : seedList) {
        yesVisit.at(seedId) = true;
    }
    ssize_t toVisit = static_cast<ssize_t>(seedList.size());
    if (direction == TraversalDirection::DIRECTION_DOWN) {
        for (NodeID nodeIdP1 = numNodes(); nodeIdP1 > 0; nodeIdP1--) {
            const NodeID nodeId = nodeIdP1 - 1;
            if (yesVisit[nodeId]) {
                assert(toVisit > 0);
                toVisit--;
                const bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_NONE);
                if (keepGoing) {
                    for (const NodeID succId : this->getDownEdges(nodeId)) {
                        if (!yesVisit[succId]) {
                            yesVisit[succId] = true;
                            toVisit++;
                        }
                    }
                }
                if (toVisit == 0) {
                    break;
                }
            }
        }
    } else {
        for (NodeID nodeId = 0; nodeId < numNodes(); nodeId++) {
            if (yesVisit[nodeId]) {
                assert(toVisit > 0);
                toVisit--;
                const bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_NONE);
                if (keepGoing) {
                    for (const NodeID succId : this->getUpEdges(nodeId)) {
                        if (!yesVisit[succId]) {
                            yesVisit[succId] = true;
                            toVisit++;
                        }
                    }
                }
                if (toVisit == 0) {
                    break;
                }
            }
        }
    }
}

void GRG::visitTopoNodeOrderedSparse(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList) {
    GRGPtr sharedThis = shared_from_this();
    NodeIDSet alreadySeen;
    const NodeIDSizeT numNodes = this->numNodes();

#define HEAP_ID(nodeId) (direction == DIRECTION_DOWN) ? (nodeId) : (numNodes - (nodeId))

    std::vector<NodeID> heap;
    for (const NodeID& seedId : seedList) {
        heap.emplace_back(HEAP_ID(seedId));
    }

#ifndef NDEBUG
    NodeIDSizeT prevPopped = (direction == DIRECTION_UP) ? 0 : numNodes;
#endif
    std::make_heap(heap.begin(), heap.end());
    while (!heap.empty()) {
        std::pop_heap(heap.begin(), heap.end());
        const NodeID nodeId = HEAP_ID(heap.back());
#ifndef NDEBUG
        if (direction == DIRECTION_UP) {
            assert(prevPopped <= nodeId);
        } else {
            assert(prevPopped >= nodeId);
        }
        prevPopped = nodeId;
#endif
        heap.pop_back();

        auto insertIt = alreadySeen.insert(nodeId);
        if (!insertIt.second) { // Was already inserted? Then skip this node.
            continue;
        }

        const bool keepGoing = visitor.visit(sharedThis, nodeId, direction, DFS_PASS_NONE);
        if (keepGoing) {
            const auto& successors =
                (direction == DIRECTION_UP) ? this->getUpEdges(nodeId) : this->getDownEdges(nodeId);
            for (const auto& succId : successors) {
                heap.emplace_back(HEAP_ID(succId));
                std::push_heap(heap.begin(), heap.end());
            }
        }
    }
}

void CSRGRG::visitTopo(GRGVisitor& visitor,
                       TraversalDirection direction,
                       const NodeIDList& seedList,
                       const std::vector<NodeIDSizeT>* sortOrder) {
    if (sortOrder != nullptr) {
        throw ApiMisuseFailure("CSRGRG has its node topologically ordered; sortOrder unneeded");
    }
    if (seedList.size() > TOPO_DENSE_THRESHOLD) {
        this->visitTopoNodeOrderedDense(visitor, direction, seedList);
    } else {
        this->visitTopoNodeOrderedSparse(visitor, direction, seedList);
    }
}

void MutableGRG::compact(NodeID nodeId) {
    if (nodeId == INVALID_NODE_ID) {
        constexpr size_t CONST_FACTOR = 2;
        for (auto& node : m_nodes) {
            if (node->m_downEdges.capacity() > node->m_downEdges.size() + CONST_FACTOR) {
                node->m_downEdges.shrink_to_fit();
            }
            if (node->m_upEdges.capacity() > node->m_upEdges.size() + CONST_FACTOR) {
                node->m_upEdges.shrink_to_fit();
            }
        }
    } else {
        m_nodes.at(nodeId)->m_downEdges.shrink_to_fit();
    }
}

} // namespace grgl
