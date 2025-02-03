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
#ifndef GRG_HAP_INDEX_H
#define GRG_HAP_INDEX_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <list>
#include <unordered_map>
#include <vector>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "hap_helpers.h"
#include "lean_bk_tree.h"
#include "util.h"

namespace grgl {

using NodeToHapVect = std::vector<HaplotypeVector>;

/**
 * Index a dataset based on a vector representing each haplotype.
 */
class HaplotypeIndex {
public:
    explicit HaplotypeIndex(std::function<size_t(const NodeID&, const NodeID&)> distFunc)
        : m_bkTree(std::move(distFunc)) {}

    virtual ~HaplotypeIndex() = default;

    HaplotypeIndex(HaplotypeIndex&) = delete;
    HaplotypeIndex(HaplotypeIndex&&) = default;
    HaplotypeIndex& operator=(HaplotypeIndex&) = delete;
    HaplotypeIndex& operator=(HaplotypeIndex&&) = default;

    /**
     * Add a new (hash, node) pair to the index.
     */
    void add(NodeID nodeId);

    /**
     * Remove a (hash, node) pair from the index.
     */
    void remove(NodeID nodeId);

    /**
     * Find the nearest neighbor to targetHash (distance=D) and then find all other neighbors
     * at the same distance D and return the list.
     */
    NodeIDList getMostSimilarNodes(NodeID nodeId, bool collectAll);

    void emitStats() const {
        std::cout << " -- Index Stats --" << std::endl;
        std::cout << "  -> Comparisons: " << m_comparisons << std::endl;
    }

private:
    LeanBKTree<NodeID> m_bkTree;
    // Keep track of how many comparisons we do.
    size_t m_comparisons{};
};

} // namespace grgl

#endif /* GRG_HAP_INDEX_H */
