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
#include "hap_index.h"

#include <limits>
#include <unordered_map>

#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "picovcf.hpp"
#include "util.h"

namespace grgl {

void dumpHash(const HaplotypeVector& hash) {
    for (auto part : hash) {
        std::cout << part << ", ";
    }
    std::cout << std::endl;
}

void HaplotypeIndex::add(const NodeID nodeId) { m_bkTree.insert(nodeId, m_comparisons); }

void HaplotypeIndex::remove(const NodeID nodeId) {}

NodeIDList HaplotypeIndex::getMostSimilarNodes(const NodeID nodeId, const bool collectAll) {
    /* There are a lot of considerations/parameters here:
     * - Only get nodes at the distance of the nearest neighbor, or allow a larger radius?
     * - Limit the number of neighbors? How do we choose the limit (is it per-level)?
     * - There can be hash collisions; improving the hash function would probably help, but
     *   also we filter the results checking the hamming distance. We could try different hash
     *   functions and also _not_ checking the hamming distance (time overhead seems low though)
     *
     * One thing I've learned is that you _have_ to tune against datasets with large samplesets,
     * as you get fairly different results for what is optimal, and that is the case where size
     * of the GRG really matters.
     */
    NodeIDList result;

    size_t nearestDistance = 0;
    auto bkTreeNodes = m_bkTree.lookup(nodeId, nearestDistance, m_comparisons, collectAll);
    for (auto& node : bkTreeNodes) {
        node->moveElements(result);
    }
    for (auto it = result.begin(); it != result.end();) {
        if (*it == nodeId) {
            it = result.erase(it);
        } else {
            ++it;
        }
    }
    return result;
}

} // namespace grgl
