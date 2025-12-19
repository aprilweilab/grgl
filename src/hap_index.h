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
#include <limits>
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
    explicit HaplotypeIndex(std::function<size_t(const NodeID&, const NodeID&)> distFunc,
                            const double rebuildProportion = 2.0)
        : m_bkTree(std::move(distFunc)),
          m_rebuildProportion(rebuildProportion),
          m_maxComparisons(std::numeric_limits<size_t>::max()),
          m_comparisonBudget(std::numeric_limits<size_t>::max()) {}

    virtual ~HaplotypeIndex() = default;

    HaplotypeIndex(HaplotypeIndex&) = delete;
    HaplotypeIndex(HaplotypeIndex&&) = default;
    HaplotypeIndex& operator=(HaplotypeIndex&) = delete;
    HaplotypeIndex& operator=(HaplotypeIndex&&) = delete;

    void setMaxComparePerQuery(const size_t maxCompare) {
        m_maxComparisons = maxCompare;
        // We give an initial budget of 5 times the average, so the first couple queries can
        // work harder if needed.
        m_comparisonBudget = std::max(maxCompare, maxCompare * 5);
    }

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

    void emitStats(std::ostream& stream) const {
        stream << " -- Index Stats --" << std::endl;
        stream << "  -> Comparisons: " << m_comparisons << std::endl;
        stream << "  -> Queries: " << m_queries << std::endl;
        stream << "  -> Truncated Queries: " << m_trunc << std::endl;
        m_bkTree.dumpStats(stream);
    }

private:
    std::vector<std::shared_ptr<LeanBKTreeNode<NodeID>>> m_nodeToBKNode;
    LeanBKTree<NodeID> m_bkTree;
    const double m_rebuildProportion;
    // Keep track of how many comparisons we do.
    size_t m_comparisons{};
    // How many comparison we're allowed, while maintaining the per-query average.
    size_t m_comparisonBudget = std::numeric_limits<size_t>::max();
    // Max comparisons allowed per query
    size_t m_maxComparisons;

    size_t m_queries{};
    size_t m_trunc{};
};

} // namespace grgl

#endif /* GRG_HAP_INDEX_H */
