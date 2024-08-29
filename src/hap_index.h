#ifndef GRG_HAP_INDEX_H
#define GRG_HAP_INDEX_H

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <list>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "util.h"
#include "lean_bk_tree.h"
#include "hap_helpers.h"

namespace grgl {

using NodeToHapVect = std::vector<HaplotypeVector>;

/**
 * Index a dataset based on a vector representing each haplotype.
 */
class HaplotypeIndex {
public:
    explicit HaplotypeIndex(std::function<size_t(const NodeID&, const NodeID&)> distFunc)
        : m_bkTree(std::move(distFunc)) {
    }

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

}

#endif /* GRG_HAP_INDEX_H */