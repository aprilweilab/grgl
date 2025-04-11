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
#include "util.h"
#include <algorithm>

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

#ifdef COMPACT_NODE_IDS
using NodeID = uint32_t;
using NodeIDSizeT = uint32_t;

// We support about 2 billion nodes per graph.
static constexpr NodeIDSizeT MAX_GRG_NODES = 0x7fffffffU - 1;
static constexpr NodeID INVALID_NODE_ID = 0x7fffffffU;
// The upper bit is available for flags.
static constexpr NodeID GRG_NODE_FLAG_MASK = 0x80000000U;

using NodeMark = NodeID;
constexpr NodeMark NODE_MARK_1 = 0x80000000U;
#else
using NodeID = uint64_t;
using NodeIDSizeT = uint64_t;

// We support about trillion nodes.
#define MAX_GRG_NODES      (0x000000ffffffffff - 1)
#define INVALID_NODE_ID    (0x000000ffffffffff)
// The upper 24 bits are available for flags.
#define GRG_NODE_FLAG_MASK (0xffffff0000000000)

enum NodeMark {
    NODE_MARK_1 = 0x0000010000000000,
    NODE_MARK_2 = 0x0000020000000000,
    NODE_MARK_3 = 0x0000040000000000,
    NODE_MARK_4 = 0x0000080000000000,
};
#endif

static_assert((INVALID_NODE_ID & GRG_NODE_FLAG_MASK) == 0, "Invalid masks");

using NodeIDList = std::vector<NodeID>;
using NodeIDSetOrdered = std::set<NodeID>;
#ifdef USE_TREE_SETS
using NodeIDSet = NodeIDSetOrdered;
#else
using NodeIDSet = std::unordered_set<NodeID>;
#endif

inline NodeID markNodeId(const NodeID nodeId, NodeMark markNum, bool value) {
    if (value) {
        return nodeId | markNum;
    }
    return nodeId & ~static_cast<NodeID>(markNum);
}

inline bool hasMark(const NodeID nodeId, NodeMark markNum) { return (nodeId & markNum) > 0; }

inline NodeID removeMarks(const NodeID nodeId) { return nodeId & (~GRG_NODE_FLAG_MASK); }

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
