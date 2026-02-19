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
#ifndef GRG_NODE_H
#define GRG_NODE_H

#include <cstdint>
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>

#include "mutation.h"

#define PERFORM_DUP_EDGE_CHECKS 0

#if PERFORM_DUP_EDGE_CHECKS
#include <algorithm>
#include <iostream>

#define CHECK_DUP_EDGES(edge_list)                                                                                     \
    do {                                                                                                               \
        std::sort((edge_list).begin(), (edge_list.end()));                                                             \
        const size_t origSize = (edge_list).size();                                                                    \
        (edge_list).erase(std::unique((edge_list).begin(), (edge_list).end()), (edge_list).end());                     \
        release_assert(origSize == (edge_list).size());                                                                \
    } while (0)
#else
#define CHECK_DUP_EDGES(edge_list)
#endif

namespace grgl {

using EdgeSizeT = uint64_t;

#ifdef COMPACT_NODE_IDS
using NodeID = uint32_t;
using NodeIDSizeT = uint32_t;
using SignedNodeID = int32_t;

// We support about 2 billion nodes per graph.
static constexpr NodeIDSizeT MAX_GRG_NODES = 0x7fffffffU - 1;
static constexpr NodeID INVALID_NODE_ID = 0x7fffffffU;

static constexpr NodeID GRG_NODE_IS_NEGATIVE = 0x80000000U;
static constexpr NodeID GRG_NODE_NEGATIVE_MASK = 0x7fffffffU;
#else
using NodeID = uint64_t;
using NodeIDSizeT = uint64_t;
using SignedNodeID = int64_t;

// We support about trillion nodes.
static constexpr NodeIDSizeT MAX_GRG_NODES = 0x000000ffffffffff - 1;
static constexpr NodeID INVALID_NODE_ID = 0x000000ffffffffff;

static constexpr NodeID GRG_NODE_IS_NEGATIVE = 0x8000000000000000;
static constexpr NodeID GRG_NODE_NEGATIVE_MASK = 0xffffff0000000000;
#endif

static_assert((INVALID_NODE_ID & ~GRG_NODE_NEGATIVE_MASK) == 0, "Invalid masks");

using NodeIDList = std::vector<NodeID>;
using NodeIDSetOrdered = std::set<NodeID>;
#ifdef USE_TREE_SETS
using NodeIDSet = NodeIDSetOrdered;
#else
using NodeIDSet = std::unordered_set<NodeID>;
#endif

using NodeIDListUPtr = std::unique_ptr<NodeIDList>;

// Just on the off case that someone was using this externally (highly, highly unlikely).
// We used to have a bit that could be stolen generically for marking nodes.
#define NODE_MARK_1 "THIS IS NO LONGER VALID! The node marks are not part of the public API"
#define markNodeId  NODE_MARK_1
#define hasMark     NODE_MARK_1
#define removeMarks NODE_MARK_1

/**
 * Negative nodes, or "prepended nodes" are topologically beneath all the positive,
 * or "regular" nodes. This is useful for when you need to build out a graph beneath
 * the sample nodes.
 */
inline bool nodeIsNegative(const NodeID nodeId) { return (bool)(nodeId & GRG_NODE_IS_NEGATIVE); }

/**
 * Strip off the negative flag on a NodeID;
 */
inline NodeID nodeStripNegative(const NodeID nodeId) {
    if (nodeIsNegative(nodeId)) {
        return -nodeId;
    }
    return nodeId;
}

class MutableGRG;

class GRGNode;
using GRGNodePtr = std::shared_ptr<GRGNode>;

/**
 * A node in the Genomic Representation Graph (GRG). Typically only accessed via the GRG
 * class, and not independently.
 */
class GRGNode {
public:
    GRGNode() = default;

    virtual ~GRGNode() = default;
    GRGNode(const GRGNode&) = delete;
    GRGNode(GRGNode&&) = delete;
    GRGNode& operator=(const GRGNode&) = delete;
    GRGNode& operator=(GRGNode&&) = delete;

    const NodeIDList& getDownEdges() const { return this->m_downEdges; }

    const NodeIDList& getUpEdges() const { return this->m_upEdges; }

    bool operator==(const GRGNode& rhs) const = delete;

protected:
    static bool edgeDelete(NodeIDList& edges, const NodeID nodeId) {
        size_t laggingCounter = 0;
        size_t leadingCounter = 0;
        for (leadingCounter = 0; leadingCounter < edges.size(); leadingCounter++) {
            if (edges[leadingCounter] != nodeId) {
                if (laggingCounter < leadingCounter) {
                    edges[laggingCounter] = edges[leadingCounter];
                }
                laggingCounter++;
            }
        }
        if (leadingCounter != laggingCounter) {
            edges.resize(laggingCounter);
            return true;
        }
        return false;
    }

    void addDownEdge(const NodeID target) {
        this->m_downEdges.push_back(target);
        CHECK_DUP_EDGES(this->m_downEdges);
    }

    bool deleteDownEdge(const NodeID target) { return edgeDelete(this->m_downEdges, target); }

    void addUpEdge(const NodeID source) { this->m_upEdges.push_back(source); }

    bool deleteUpEdge(const NodeID source) { return edgeDelete(this->m_upEdges, source); }

    NodeIDList m_downEdges;
    NodeIDList m_upEdges;

    friend MutableGRG;
};

} // namespace grgl

#endif /* GRG_NODE_H */
