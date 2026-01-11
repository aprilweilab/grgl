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
#ifndef GRG_VISITOR_H
#define GRG_VISITOR_H

#include <memory>

#include "grgnode.h"

namespace grgl {

enum TraversalDirection { DIRECTION_DOWN = 1, DIRECTION_UP = 2 };

enum DfsPass { DFS_PASS_NONE = 0, DFS_PASS_THERE = 1, DFS_PASS_BACK_AGAIN = 2 };

class GRG;
using GRGPtr = std::shared_ptr<GRG>;

/**
 * Visitor abstract base class for traversing the nodes of a GRG.
 *
 * To traverse a GRG, create a class that inherits from `GRGVisitor` and overrides the `visit()`
 * function. Then use the GRG::visitDfs(), GRG::visitBfs(), or GRG::visitTopo() functions and
 * pass your visitor as an argument.
 */
class GRGVisitor {
public:
    GRGVisitor() = default;
    virtual ~GRGVisitor() = default;

    GRGVisitor(const GRGVisitor&) = delete;
    GRGVisitor(GRGVisitor&&) = delete;
    GRGVisitor& operator=(const GRGVisitor&) = delete;
    GRGVisitor& operator=(GRGVisitor&&) = delete;

    /**
     * Visit a node in the GRG.
     * @param[in] grg The graph.
     * @param[in] node The current node we are visiting.
     * @param[in] direction Whether we're traversing bottom-up (DIRECTION_UP) or
     *      top-down (DIRECTION_DOWN).
     * @param[in] dfsPass For DFS traversals, is this the forward (DFS_PASS_THERE)
     *      pass or the backward (DFS_PASS_BACK_AGAIN) pass after we have computed
     *      values for all successors.
     */
    virtual bool
    visit(const GRGPtr& grg, NodeID node, TraversalDirection direction, DfsPass dfsPass = DFS_PASS_NONE) = 0;
};

} // namespace grgl

#endif /* GRG_VISITOR_H */
