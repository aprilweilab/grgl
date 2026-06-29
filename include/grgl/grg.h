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
#ifndef GRG_H
#define GRG_H

#include <algorithm>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/file_vector.h"
#include "grgnode.h"
#include "mutation.h"
#include "node_data.h"
#include "visitor.h"

#if USE_AVX
#include <immintrin.h>
#endif

// To avoid dependency on gtest.
#ifndef FRIEND_TEST
#define FRIEND_TEST(x, y)
#endif

namespace grgl {

// Returned by the edge-count functionality when there are no edges loaded/stored for the up direction, as opposed to
// actually having no edges. This is safer than returning 0 for the number of edges, as anything that requires up
// will fail hard trying to iterate this many edges.
constexpr size_t NO_UP_EDGES = std::numeric_limits<size_t>::max();

/**
 * Simple class to encapsulate the filtering of a GRG to shrink it.
 */
class GRGOutputFilter {
public:
    GRGOutputFilter()
        : direction(TraversalDirection::DIRECTION_UP) {}

    GRGOutputFilter(TraversalDirection dir, NodeIDList seeds)
        : direction(dir),
          seedList(std::move(seeds)) {}

    bool isSpecified() const { return !seedList.empty(); }

    // Traverse UP or DOWN from the seeds.
    TraversalDirection direction;
    // These are the nodes to start traversing from.
    NodeIDList seedList;
    std::pair<BpPosition, BpPosition> bpRange;
};

/**
 * Abstract GRG base class.
 */
class GRG : public std::enable_shared_from_this<GRG> {
public:
    static constexpr size_t DEFAULT_NODE_CAPACITY = 1024;

    /**
     * A node and a mutation attached to that node.
     */
    using NodeAndMut = std::pair<NodeID, MutationId>;

    /**
     * A node and a mutation attached to that node, and the missingness node associated with the mutation.
     */
    using NodeMutMiss = std::tuple<NodeID, MutationId, NodeID>;

    /**
     * A mutation and node, where the mutation comes first.
     */
    using MutAndNode = std::pair<MutationId, NodeID>;

    /**
     * A mutation, node, and missingness node, where the mutation comes first.
     */
    using MutNodeMiss = std::tuple<MutationId, NodeID, NodeID>;

    explicit GRG(size_t numSamples, uint16_t ploidy, const bool phased = true)
        : m_numSamples(numSamples),
          m_ploidy(ploidy),
          m_phased(phased) {
        release_assert(numSamples % ploidy == 0);
    }

    virtual ~GRG() = default;
    GRG(const GRG&) = delete;
    GRG(GRG&&) = delete;
    GRG& operator=(const GRG&) = delete;
    GRG& operator=(GRG&&) = delete;

    // Don't compare GRGs. Use equivalentGRGs if needed.
    bool operator==(const GRG& rhs) const = delete;

    /**
     * Is the given nodeID associated with a sample node?
     *
     * In most GRGs, the first numSamples node IDs are reserved for samples,
     * so you can also just use `nodeId < grg.numSamples()` _if_ you are using an
     * immutable GRG (CSRGRG) or more generally if samplesAreOrdered() returns true.
     *
     * @param[in] nodeId The node ID to check.
     */
    virtual bool isSample(const NodeID nodeId) const { return nodeId < this->m_numSamples; }

    /**
     * The number of samples (haplotypes) in the graph.
     */
    virtual NodeIDSizeT numSamples() const { return m_numSamples; }

    /**
     * The number of individuals (haplotypes/ploidy) in the graph.
     */
    NodeIDSizeT numIndividuals() const { return numSamples() / m_ploidy; }

    /**
     * How many haploid samples are there per individual?
     *
     * @return The ploidy, usually 1 or 2. Individual coalescence support only works when ploidy==2.
     */
    uint16_t getPloidy() const { return m_ploidy; }

    /**
     * Is the dataset phased?
     *
     * If not, GRG still represents the data as a mapping between mutations and haploid samples, but
     * you cannot depend on any haplotype-based analyses representing the true haplotypes. Individual-based
     * analyses are still completely accurate.
     */
    bool isPhased() const { return m_phased; }

    /**
     * Returns true if nodes are ordered by ID in bottom-up topological order.
     *
     * This means you can just iterate the node IDs from 0...(numNodes - 1) to
     * perform a topological traversal.
     */
    virtual bool nodesAreOrdered() const = 0;

    /**
     * Returns true if nodes are ordered in bottom-up topological order, If this is true
     * and nodesAreOrdered is false, then you cannot iterate by NodeID (which are
     * always positive), but you can iterate by using getOrderedNodes().
     *
     * Most users will not need this method, unless they are heavily modifying the GRG
     * (e.g., by adding negative nodes to the graph).
     */
    virtual bool nodesAreTopo() const = 0;

    /**
     * Number of nodes (including sample nodes) in the graph.
     */
    virtual size_t numNodes() const = 0;

    /**
     * Number of edges (total) in the graph.
     */
    virtual size_t numEdges() const = 0;

    /**
     * Number of down edges for a particular node in the graph.
     */
    virtual size_t numDownEdges(NodeID nodeId) = 0;

    /**
     * Number of up edges for a particular node in the graph.
     */
    virtual size_t numUpEdges(NodeID nodeId) = 0;

    /**
     * Get list of down edges for a particular node in the graph.
     */
    virtual NodeIDList getDownEdges(NodeID nodeId) = 0;

    /**
     * Get list of up edges for a particular node in the graph.
     */
    virtual NodeIDList getUpEdges(NodeID nodeId) = 0;

    /**
     * Get all of the node IDs in topological order for the given direction (UP means do
     * the leaves first, DOWN means do the roots first).
     *
     * @param[in] direction Direction for the order.
     * @param[in] allowSort Set to true to allow performing a topological sort to get the result,
     *      if the nodes are not already ordered. This can impose a substantial performance hit
     *      and is only necessary when the graph has been modified in ways that did not preserve
     *      the topological order! It is almost always possible to modify the graph so that the
     *      topological order is preserved, so this option should rarely (if ever) be used.
     *      Default: false.
     */
    virtual NodeIDList getOrderedNodes(TraversalDirection direction, bool allowSort = false) = 0;

    /**
     * True if this graph stores up edges, false if only down edges are stored.
     */
    virtual bool hasUpEdges() const = 0;

    /**
     * True if the edges in this graph are always in ascending NodeID order.
     */
    virtual bool edgesAreOrdered() const = 0;

    /**
     * True iff the mutations are ordered according to position and the alternate allele (in that order).
     */
    bool mutationsAreOrdered() const { return m_mutsAreOrdered; }

    /**
     * Get the base-pair range this GRG covers as a pair {first, last+1} where first is
     * the base-pair position of the first mutation, and last is similarly defined.
     * In general we treat all ranges as left-inclusive and right-exlusive, hence the
     * +1 on the last.
     */
    std::pair<BpPosition, BpPosition> getBPRange() {
        if (mutationsAreOrdered() && !m_mutations.empty()) {
            const size_t lastMutId = m_mutations.size() - 1;
            return {m_mutations[0].getPosition(), m_mutations[lastMutId].getPosition() + 1};
        }
        BpPosition firstPos = std::numeric_limits<BpPosition>::max();
        BpPosition lastPos = 0;
        for (MutationId i = 0; i < numMutations(); i++) {
            const auto position = m_mutations[i].getPosition();
            if (position < firstPos) {
                firstPos = position;
            }
            if (position > lastPos) {
                lastPos = position;
            }
        }
        return {firstPos, lastPos + 1};
    }

    /**
     * Get the specified base-pair range this GRG covers as a pair {first, last+1}. Unlike getBPRange,
     * this function just returns a pair of values that were specified by a tool/user at some point.
     */
    std::pair<BpPosition, BpPosition> getSpecifiedBPRange() const { return m_specifiedRange; }

    /**
     * Set the specified base-pair range this GRG covers as a pair {first, last+1}, i.e. the lower value
     * is inclusive and the higher value is exclusive.
     */
    void setSpecifiedBPRange(std::pair<BpPosition, BpPosition> bpRange) { m_specifiedRange = bpRange; }

    /**
     * Get the NodeIDList for all samples in the graph, in order of sample ID (from the 0th sample
     * to the (N-1)th sample, where N is the number of haplotypes).
     */
    virtual NodeIDList getSampleNodes() const {
        NodeIDList result(this->numSamples());
        std::iota(result.begin(), result.end(), 0);
        return std::move(result);
    }

    /**
     * Get the NodeIDList for all roots in the graph.
     * Note: this is an O(numNodes) operation.
     */
    NodeIDList getRootNodes() {
        NodeIDList result;
        if (hasUpEdges()) {
            for (grgl::NodeID nodeId = 0; nodeId < this->numNodes(); nodeId++) {
                if (this->numUpEdges(nodeId) == 0) {
                    result.push_back(nodeId);
                }
            }
        } else {
            std::vector<bool> isTarget(numNodes(), false);
            for (grgl::NodeID nodeId = 0; nodeId < this->numNodes(); nodeId++) {
                for (grgl::NodeID targetId : getDownEdges(nodeId)) {
                    isTarget[targetId] = true;
                }
            }
            for (grgl::NodeID nodeId = 0; nodeId < this->numNodes(); nodeId++) {
                if (!isTarget[nodeId]) {
                    result.push_back(nodeId);
                }
            }
        }
        return std::move(result);
    }

    const std::vector<Mutation>& getMutations() const { return m_mutations.vector(); }

    size_t numMutations() const { return m_mutations.size(); }

    /**
     * Does this GRG have Mutations which have missing data represented?
     */
    bool hasMissingData() const { return !m_mutIdsByNodeIdAndMiss.empty(); }

    /**
     * Iterator for traversing the NodeID<-->MutationID mapping, along with the (sometimes present) NodeID
     * for the missingness node associated with the Mutation.
     *
     * Each Mutation can have two different "kinds" of nodes associated with it:
     * 1. The Mutation Node, which defines the mapping between a Mutation and all of the samples that have
     *    that Mutation.
     * 2. The Missingness Node, which defines the mapping between a Mutation and all of the samples that
     *    are missing an allele at the site associated with the Mutation.
     */
    template <typename T> class NodeAndMutIterator {
    public:
        class VectIter3 {
        public:
            VectIter3(const std::vector<NodeAndMut>* pair, size_t pos)
                : m_pair(pair),
                  m_position(pos) {}
            VectIter3(const std::vector<NodeMutMiss>* triple, size_t pos)
                : m_triple(triple),
                  m_position(pos) {}
            VectIter3& operator++() {
                m_position++;
                return *this;
            }
            VectIter3 operator++(int) {
                VectIter3 retval = *this;
                ++(*this);
                return retval;
            }
            bool operator==(VectIter3 other) const {
                return m_pair == other.m_pair && m_triple == other.m_triple && m_position == other.m_position;
            }
            bool operator!=(VectIter3 other) const { return !(*this == other); }
            T operator*();
            size_t size() const { return (nullptr != m_pair) ? m_pair->size() : m_triple->size(); }

            using difference_type = size_t;
            using value_type = T;
            using pointer = T*;
            using reference = const T&;
            using iterator_category = std::forward_iterator_tag;

        private:
            const std::vector<NodeAndMut>* m_pair{};
            const std::vector<NodeMutMiss>* m_triple{};
            size_t m_position{};
        };

        explicit NodeAndMutIterator(const std::vector<NodeAndMut>* pair)
            : m_begin(pair, 0),
              m_end(pair, pair->size()) {}

        explicit NodeAndMutIterator(const std::vector<NodeMutMiss>* triple)
            : m_begin(triple, 0),
              m_end(triple, triple->size()) {}

        /**
         * The number of items in the vector underlying the iterator.
         */
        size_t size() const { return m_begin.size(); }

        /**
         * The beginning of the vector.
         */
        VectIter3 begin() const { return m_begin; }
        /**
         * The end of the vector.
         */
        VectIter3 end() const { return m_end; }

    private:
        VectIter3 m_begin;
        VectIter3 m_end;
    };

    // DEPRECATED! Use getNodesAndMutations() instead.
    NodeAndMutIterator<NodeAndMut> getNodeMutationPairs() GRGL_DEPRECATED { return getNodesAndMutations<NodeAndMut>(); }

    /**
     * Get the list of all corresponding nodes and mutations, for iterating over the mapping
     * between graph nodes and Mutations. By default, just returns the (NodeID, MutationID) pairs,
     * but using the template argument of GRG::NodeMutMiss will return a tuple (NodeID, MutationID, NodeID)
     * where the additional NodeID is the missingness node (or INVALID_NODE_ID) associated with this
     * Mutation.
     *
     * @param[in] allowSort Allow the mutations-to-node mapping to be resorted; if false and the mapping is
     *      not sorted then an exception will be thrown. Sorting repeatedly is slow (obviously), so this flag
     *      can help avoid/debug performance issues.
     *
     * @return An iterator over (NodeID, MutationId) pairs or (NodeID, MutationId, NodeID) tuples,
     *      depending on the function template argument.
     */
    template <typename T = NodeAndMut> NodeAndMutIterator<T> getNodesAndMutations(const bool allowSort = true) {
        if (!m_mutIdsByNodeIdSorted) {
            api_exc_check(allowSort, "Mutation-to-node mapping needs to be sorted, but allowSort is false.");
            sortMutIdsByNodeID();
        }
        if (!m_mutIdsByNodeId.empty()) {
            return NodeAndMutIterator<T>(&m_mutIdsByNodeId);
        }
        return NodeAndMutIterator<T>(&m_mutIdsByNodeIdAndMiss);
    }

    const Mutation& getMutationById(MutationId mutId) {
        api_exc_check(mutId < m_mutations.size(), "Invalid MutationID: " << mutId);
        return m_mutations.cref(mutId);
    }

    void setMutationById(MutationId mutId, Mutation mutation) {
        api_exc_check(mutId < m_mutations.size(), "Invalid MutationID: " << mutId);
        m_mutations.ref(mutId) = std::move(mutation);
        m_mutsAreOrdered = false;
    }

    /**
     * Get the mutations associated with no nodes. By default, the template argument is MutationId,
     * so only the vector of mutations is returned. If you set the template argument to
     * GRG::MutAndNode then the mutation and it's corresponding missingness node will be returned.
     *
     * @return A vector of either MutationId or GRG::MutAndNode (a std::pair), depending on the template
     *      parameter.
     */
    template <typename T = MutationId> std::vector<T> getUnmappedMutations() {
        return getMutationsForNode<T>(INVALID_NODE_ID);
    }

    /**
     * Get the mutations associated with the node. By default, the template argument is MutationId,
     * so only the vector of mutations is returned. If you set the template argument to
     * GRG::MutAndNode then the mutation and it's corresponding missingness node will be returned.
     *
     * @param[in] nodeId The node to lookup.
     * @param[in] allowSort Allow the mutations-to-node mapping to be resorted; if false and the mapping is
     *      not sorted then an exception will be thrown. Sorting repeatedly is slow (obviously), so this flag
     *      can help avoid/debug performance issues.
     * @return A vector of either MutationId or GRG::MutAndNode (a std::pair), depending on the template
     *      parameter.
     */
    template <typename T = MutationId>
    std::vector<T> getMutationsForNode(const NodeID nodeId, const bool allowSort = true) {
        if (!m_mutIdsByNodeIdSorted) {
            api_exc_check(allowSort, "Mutation-to-node mapping needs to be sorted, but allowSort is false.");
            sortMutIdsByNodeID();
        }
        std::vector<T> result;
        if (!m_mutIdsByNodeId.empty()) {
            const NodeAndMut query = {nodeId, 0};
            auto mutIdIt = std::lower_bound(m_mutIdsByNodeId.begin(), m_mutIdsByNodeId.end(), query, NodeAndMutLt());
            while (mutIdIt != m_mutIdsByNodeId.end() && mutIdIt->first == nodeId) {
                if (mutIdIt->second != INVALID_MUTATION_ID) {
                    result.push_back(std::move(nodeAndMutToType<T>(*mutIdIt)));
                }
                mutIdIt++;
            }
        } else {
            const NodeMutMiss query = {nodeId, 0, 0};
            auto mutIdIt = std::lower_bound(
                m_mutIdsByNodeIdAndMiss.begin(), m_mutIdsByNodeIdAndMiss.end(), query, NodeMutMissLt());
            while (mutIdIt != m_mutIdsByNodeIdAndMiss.end() && std::get<0>(*mutIdIt) == nodeId) {
                if (std::get<1>(*mutIdIt) != INVALID_MUTATION_ID) {
                    result.push_back(nodeAndMutAndMissToType<T>(*mutIdIt));
                }
                mutIdIt++;
            }
        }
        return std::move(result);
    }

    /**
     * Does the node have at least one Mutation associated with it?
     *
     * @param[in] nodeId The node to lookup.
     * @param[in] allowSort Allow the mutations-to-node mapping to be resorted; if false and the mapping is
     *      not sorted then an exception will be thrown. Sorting repeatedly is slow (obviously), so this flag
     *      can help avoid/debug performance issues.
     * @return true if the node has one or more Mutation.
     */
    bool nodeHasMutations(const NodeID nodeId, const bool allowSort = true) {
        if (!m_mutIdsByNodeIdSorted) {
            api_exc_check(allowSort, "Mutation-to-node mapping needs to be sorted, but allowSort is false.");
            sortMutIdsByNodeID();
        }
        if (!m_mutIdsByNodeId.empty()) {
            const std::pair<NodeID, MutationId> query = {nodeId, 0};
            auto mutIdIt = std::lower_bound(m_mutIdsByNodeId.begin(), m_mutIdsByNodeId.end(), query, NodeAndMutLt());
            while (mutIdIt != m_mutIdsByNodeId.end() && mutIdIt->first == nodeId) {
                if (mutIdIt->second != INVALID_MUTATION_ID) {
                    return true;
                }
                mutIdIt++;
            }
            return false;
        }
        const std::tuple<NodeID, MutationId, NodeID> query = {nodeId, 0, 0};
        auto mutIdIt =
            std::lower_bound(m_mutIdsByNodeIdAndMiss.begin(), m_mutIdsByNodeIdAndMiss.end(), query, NodeMutMissLt());
        while (mutIdIt != m_mutIdsByNodeIdAndMiss.end() && std::get<0>(*mutIdIt) == nodeId) {
            if (std::get<1>(*mutIdIt) != INVALID_MUTATION_ID) {
                return true;
            }
            mutIdIt++;
        }
        return false;
    }

    /**
     * Get corresponding mutation IDs and node IDs, ordered by the mutation position + allele (ascending).
     * By default, a MutAndNode is returned (std::pair). If you set the template argument to MutIdNodeMiss
     * (std::tuple) then the third item in the returned tuple will be the missingness node associated with
     * the Mutation.
     *
     * @return A vector of pairs, MutationID and NodeID (in that order), or a vector of tuples of
     *      (MutationID, NodeID, NodeID), depending on function template argument.
     */
    template <typename T = MutAndNode> std::vector<T> getMutationsToNodeOrdered() {
        std::vector<T> result;
        if (!m_mutIdsByNodeId.empty()) {
            for (const auto& nodeAndMutId : m_mutIdsByNodeId) {
                if (nodeAndMutId.second != INVALID_MUTATION_ID) {
                    result.emplace_back(std::move(nodeAndMutToRevType<T>(nodeAndMutId)));
                }
            }
        } else {
            for (const auto& triple : m_mutIdsByNodeIdAndMiss) {
                if (std::get<1>(triple) != INVALID_MUTATION_ID) {
                    result.emplace_back(nodeAndMutAndMissToRevType<T>(triple));
                }
            }
        }
        sortByMutation(result);
        return std::move(result);
    }

    /**
     * Visit nodes breadth-first, starting at the given nodes and following up or
     * down edges.
     *
     * @param[in] visitor The visitor that will be called back. Typically owns any state
     *      that is being computed.
     * @param[in] direction Whether to go up or down the graph.
     * @param[in] seedList The nodes to start traversal from.
     * @param[in] maxQueueWidth The maximum width of the queue; restricts the number of
     *      end-to-end paths that will be visited.
     */
    void
    visitBfs(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList, ssize_t maxQueueWidth = -1);

    /**
     * Visit nodes depth-first, starting at the given nodes and following up or
     * down edges. A node is never visited more than one time in either pass (forward
     * or backward), unless forwardOnly is true.
     *
     * @param[in] visitor The visitor that will be called back. Typically owns any state
     *      that is being computed.
     * @param[in] direction Whether to go up or down the graph.
     * @param[in] seedList The nodes to start traversal from.
     * @param[in] maxPathDepth The maximum depth of any path that will be explored.
     *      Traversal stops when this is reached. Typically only useful when forwardOnly
     *      is true.
     * @param[in] forwardOnly Typically depth-first search reaches a given node, computes
     *      (recursively) the values for all successors, and then computes the current
     *      node value. This causes each node to be visited twice. However, setting
     *      forwardOnly will only visit nodes in the forward direction. It also causes
     *      nodes to be visited an arbitrary number of times.
     */
    void
    visitDfs(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList, bool forwardOnly = false);

    virtual std::vector<NodeIDSizeT> topologicalSort(TraversalDirection direction) = 0;

    virtual void visitTopo(GRGVisitor& visitor,
                           TraversalDirection direction,
                           const NodeIDList& seedList,
                           const std::vector<NodeIDSizeT>* sortOrder = nullptr) = 0;

    /**
     * Add a population to the GRG.
     */
    size_t addPopulation(std::string populationDescription) {
        const auto popId = m_populations.size();
        m_populations.push_back(std::move(populationDescription));
        return popId;
    }

    /**
     * The populations (their descriptions) represented by this GRG.
     */
    const std::vector<std::string>& getPopulations() const { return m_populations; }

    /**
     * Add the given mutation to the GRG, associated with the given node. If there is no node associated
     * with the mutation, use INVALID_NODE_ID.
     *
     * @param[in] mutation The Mutation to add.
     * @param[in] nodeId The nodeId to associated it with, or INVALID_NODE_ID.
     * @param[in] missNodeId The nodeId for the missingness node of this Mutation, or INVALID_NODE_ID if none.
     * @return The MutationId for the newly added Mutation.
     */
    MutationId addMutation(const Mutation& mutation, NodeID nodeId, NodeID missNodeId = INVALID_NODE_ID);

    /**
     * Remove the Mutation with the given ID and nodeId from the graph. If the Mutation has a node
     * associated with it, the mapping of that Mutation to its node will also be removed. This is a memory
     * efficient operation; it leaves an empty Mutation at the MutationID. Therefore, if you want the
     * MutationID ordering to match position/allele ordering you must call sortMutations().
     *
     * @param[in] mutId The MutationId.
     * @param[in] nodeId The NodeId.
     */
    void removeMutation(MutationId mutId, NodeID nodeId);

    /**
     * If needed, sort the mutations by their (position, allele) and renumber them so
     * that the MutationId ascending order matches this order. Can be a memory-intensive
     * operation. Alternatively, you can just write the GRG to disk (saveGrg()) and then
     * reload it (loadImmutableGRG()) and it will use less memory (but be slower).
     */
    void sortMutations();

    /**
     * Check whether the Mutations in GRG are unique: there is a single Mutation (with MutationId)
     * for each unique (position, ref allele, alt allele) combination.
     * Since the mapping between MutationId and nodes is many-to-1 (not many-to-many), having unique
     * mutations means that each Mutation is associated with a single node in the graph. When a GRG
     * is constructed via `grg construct` this is always true of the resulting GRG. When a GRG
     * is constructed arbitrarily from the API, it may not be true.
     *
     * This operation is O(M), where M is the number of mutations.
     */
    bool mutationsAreUnique() const;

    // Internal enum used for determining the node value initialization behavior of matrix multiplication.
    enum NodeInitEnum {
        NIE_ZERO = 0,   // Zero values.
        NIE_VECTOR = 1, // A vector of values, one for each row of the input matrix.
        NIE_MATRIX = 2, // A matrix of values, of dimensions ROWS x NODES
        NIE_XTX = 3,    // The number of coalescences * 2.
    };

    // Internal method, used by Python API and the main API below.
    // The outputMatrix is added to. So if you want the pure matrix multiplication, it must
    // be zeroed out.
    template <typename IOType, typename NodeValueType, bool useBitVector>
    void matrixMultiplication(const IOType* inputMatrix,
                              size_t inputCols,
                              size_t inputRows,
                              TraversalDirection direction,
                              IOType* outputMatrix,
                              size_t outputSize,
                              bool emitAllNodes = false,
                              bool byIndividual = false,
                              const IOType* initMatrix = nullptr,
                              NodeInitEnum nodeInit = NIE_ZERO,
                              IOType* missMatrix = nullptr);

    /**
     * Compute one of two possible matrix multiplications across the entire
     * graph. The input matrix \f$V\f$ can be either \f$K \times N\f$ (\f$N\f$
     * is number of samples) or \f$K \times M\f$ (\f$M\f$ is number of
     * mutations). The given direction determines which input matrix is
     * expected. Let \f$X\f$ be the \f$N \times M\f$ genotype matrix. For an
     * \f$K \times N\f$ input \f$V\f$, the product performed is \f$V \times X\f$
     * which gives a \f$K \times M\f$ result. I.e., the input matrix is a column
     * per sample and the output matrix is a column per mutation. For an \f$K
     * \times M\f$ input \f$V\f$, the product performed is \f$V \times X^T\f$
     * which gives a \f$K \times N\f$ result. I.e., the input matrix is a column
     * per mutation and the output matrix is a column per sample.
     *
     * The simplest case to consider is a vector input (e.g., a \f$1 \times N\f$
     * matrix). This vector-matrix product in the graph works by seeding the
     * input nodes (samples in this example) with the corresponding values from
     * the input vector and then traversing the graph in the relevant direction
     * (up or down). The ancestor/descendant values are summed at each node,
     * until the terminal nodes (mutations in this example) are reached. The
     * values at the terminal nodes are then the output vector.
     * When a \f$K\f$-row matrix is input, instead of a vector, the only
     * difference is that each node stores \f$K\f$ values instead of 1.
     *
     * Note: the RAM used will be \f$O(K * nodes)\f$ where \f$nodes\f$ is the
     * total number of nodes in the graph.
     *
     * @param[in] inputMatrix The \f$K \times N\f$ or \f$K \times M\f$ input
     *      matrix in row-major order. For samples, the \f$i\f$th column are
     *      the values for sample \f$i\f$. Similarly for mutations, the ordering
     *      follows the MutationId ordering.
     * @param[in] numRows The number of rows in the matrix. The size of the input
     *      must be divisible by numRows, and determines the number of columns.
     * @param[in] direction Use TraversalDirection::DIRECTION_DOWN for mutations
     *      as input, and TraversalDirection::DIRECTION_UP for samples as input.
     * @return The resulting matrix as described above, in row-major order.
     */
    template <typename T>
    std::vector<T> matMul(const std::vector<T>& inputMatrix, const size_t numRows, TraversalDirection direction) {
        if (numRows == 0 || (inputMatrix.size() % numRows != 0)) {
            throw ApiMisuseFailure("inputMatrix must be divisible by numRows");
        }
        const size_t numCols = inputMatrix.size() / numRows;
        const size_t outSize =
            (numRows * ((direction == TraversalDirection::DIRECTION_DOWN) ? numSamples() : numMutations()));
        std::vector<T> result(outSize);
        matrixMultiplication<T, T, false>(inputMatrix.data(), numCols, numRows, direction, result.data(), outSize);
        return std::move(result);
    }

    /**
     * Get the ID of the population associated with the node. Will be POPULATION_UNSPECIFIED if the
     * node is not a sample, or if there is no population associated.
     * @param[in] nodeId The node to retrieve.
     */
    PopulationID getPopulationId(NodeID nodeId) { return m_nodeData.getPopId(nodeId); }

    /**
     * Set the ID of the population associated with the node.
     * @param[in] nodeId The node to modify.
     * @param[in] popId The population ID, which is already associated with a population description
     *      via addPopulation().
     */
    void setPopulationId(NodeID nodeId, PopulationID popId) {
        api_exc_check(popId < m_populations.size(),
                      "Invalid population ID " << popId << " (only have " << m_populations.size()
                                               << " populations configured)");
        m_nodeData.setPopId(nodeId, popId);
    }

    /**
     * Get the number of individuals that coalescence _at_ this node (not below it, but exactly at it).
     * @param[in] nodeId The node to retrieve.
     * @return The number of individuals. For diploid data this is a number between 0...numSamples()/2.
     */
    NodeIDSizeT getNumIndividualCoals(NodeID nodeId) {
        const auto coals = m_nodeData.getNumCoals(numSamples(), nodeId);
        assert(coals == COAL_COUNT_NOT_SET || coals <= numIndividuals());
        return coals;
    }

    /**
     * Set the number of individuals that coalesce at this node, growing the underlying
     * data size if needed.
     *
     * NOT THREADSAFE.
     *
     * @param[in] nodeId The node to modify.
     * @param[in] coals The number of individuals that coalesce, between 0...numSamples()/ploidy.
     */
    void setNumIndividualCoalsGrow(NodeID nodeId, NodeIDSizeT coals) {
        assert(coals == COAL_COUNT_NOT_SET || coals <= numIndividuals());
        const NodeIDSizeT samples = numSamples();
        m_nodeData.allocNumCoals(numNodes() - samples);
        m_nodeData.setNumCoals(samples, nodeId, coals);
    }

    /**
     * Set the number of individuals that coalesce at this node.
     *
     * The value is updated atomically, so this can be used from threaded code (as long
     * as the rest of the GRG, e.g. number of nodes, is not changing).
     *
     * @param[in] nodeId The node to modify.
     * @param[in] coals The number of individuals that coalesce, between 0...numSamples()/ploidy.
     */
    void setNumIndividualCoals(NodeID nodeId, NodeIDSizeT coals) {
        assert(coals == COAL_COUNT_NOT_SET || coals <= numIndividuals());
        m_nodeData.setNumCoals(numSamples(), nodeId, coals);
    }

    /**
     * Does this dataset have any individual coalescence information? If not, you won't be
     * able to quickly get the exact sample variance for each mutation (diploid genotypes only).
     */
    bool hasIndividualCoals() const { return m_nodeData.hasCoalCounts(); }

    /**
     * Add the next individual's identifier. Must be called in order, for individuals
     * 0...(N-1).
     *
     * @param identifier The identifier to add.
     */
    void addIndividualId(const std::string& identifier) {
        if (m_individualIds.numNodes() == 0) {
            m_individualIds = CSRStringTable(this->numIndividuals());
        }
        m_individualIds.appendData(reinterpret_cast<const uint8_t*>(identifier.c_str()), identifier.size());
    }

    /**
     * Clear all individual identifiers.
     */
    void clearIndividualIds() { m_individualIds = CSRStringTable(); }

    /**
     * Are there any string identifiers for the individuals in this GRG?
     *
     * @return true if individual identifiers are available.
     */
    bool hasIndividualIds() const { return m_individualIds.numNodes() > 0; }

    /**
     * Get the string identifier for the k'th individual. Returns empty string if there is no identifier
     * for the given individual.
     *
     * @return String identifier for the given individual.
     */
    std::string getIndividualId(NodeIDSizeT individualIndex) {
        if (!hasIndividualIds()) {
            throw ApiMisuseFailure("No individual IDs on this GRG; check hasIndividualIds()");
        }
        std::vector<uint8_t> characters;
        m_individualIds.getData(individualIndex, characters);
        return {reinterpret_cast<const char*>(characters.data()), characters.size()};
    }

    // When nodeIds change (e.g., during merging) we need to make a post-processing
    // pass over the missinginess nodes to adjust their IDs.
    void adjustMissingnessNodeIds(const size_t afterIndex, const std::vector<NodeID>& replaceMap) {
        if (!m_mutIdsByNodeIdAndMiss.empty()) {
            for (size_t i = afterIndex; i < m_mutIdsByNodeIdAndMiss.size(); i++) {
                NodeID& missingnessNode = std::get<2>(m_mutIdsByNodeIdAndMiss[i]);
                if (missingnessNode != INVALID_NODE_ID) {
                    missingnessNode = replaceMap.at(missingnessNode);
                }
            }
        }
    }

    virtual bool isMutable() const { return false; }
    virtual bool samplesAreOrdered() const { return true; }

protected:
    struct NodeAndMutLt {
        bool operator()(const GRG::NodeAndMut& lhs, const GRG::NodeAndMut& rhs) const {
            if (lhs.first == INVALID_NODE_ID && rhs.first != INVALID_NODE_ID) {
                return true;
            }
            if (lhs.first != INVALID_NODE_ID && rhs.first == INVALID_NODE_ID) {
                return false;
            }
            return lhs < rhs;
        }
    };

    struct NodeMutMissLt {
        bool operator()(const GRG::NodeMutMiss& lhs, const GRG::NodeMutMiss& rhs) const {
            if (std::get<0>(lhs) == INVALID_NODE_ID && std::get<0>(rhs) != INVALID_NODE_ID) {
                return true;
            }
            if (std::get<0>(lhs) != INVALID_NODE_ID && std::get<0>(rhs) == INVALID_NODE_ID) {
                return false;
            }
            return lhs < rhs;
        }
    };

    void visitTopoNodeOrderedDense(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList);
    void visitTopoNodeOrderedSparse(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList);
    void populateMutMissInfo();

    // m_mutIdsByNodeId gets appended to and then sorted periodically or as needed, using this
    // function. We rarely need to lookup Mutations by NodeID while modifying a GRG, so this ends
    // up being about 3x smaller than the previous multimap<> and many times faster.
    void sortMutIdsByNodeID();

    // Templated helpers for transparently handling the cases where we do or don't have
    // missing data nodes.
    template <typename T> T nodeAndMutToType(const NodeAndMut& pair);
    template <typename T> T nodeAndMutAndMissToType(const NodeMutMiss& triple);
    template <typename T> void sortByMutation(std::vector<T>&);
    template <typename T> T nodeAndMutToRevType(const GRG::NodeAndMut& pair);
    template <typename T> T nodeAndMutAndMissToRevType(const GRG::NodeMutMiss& triple);

    // The position is the mutation ID for both of these vectors.
    EagerFileVector<Mutation> m_mutations;
    // These two vectors are equivalent, except one has more info (missingness). Only one of these will
    // be non-empty. The first two items are MutationIDs ordered by their associated NodeID, ascending.
    // The optional third item is the missingness node for the mutation.
    std::vector<std::pair<NodeID, MutationId>> m_mutIdsByNodeId;
    std::vector<std::tuple<NodeID, MutationId, NodeID>> m_mutIdsByNodeIdAndMiss;
    bool m_mutIdsByNodeIdSorted{false};

    // (Optional) list of population descriptions. The position corresponds to the population
    // ID, which can be used to tag nodes.
    std::vector<std::string> m_populations;

    // The range of base-pair positions covered by this GRG, according to the user/tool that created it.
    std::pair<BpPosition, BpPosition> m_specifiedRange{};

    // Node data.
    NodeDataContainer m_nodeData;

    // Sample (individual) identifiers
    CSRStringTable m_individualIds;

    const NodeIDSizeT m_numSamples;
    const uint16_t m_ploidy;
    const bool m_phased;

    // True if the mutationId order matches the (position, allele) order.
    bool m_mutsAreOrdered{false};

    friend void readGrgCommon(const GRGFileHeader& header, const GRGPtr& grg, IFSPointer& inStream);
    friend std::pair<NodeIDSizeT, EdgeSizeT>
    simplifyAndSerialize(const GRGPtr& grg, std::ostream& outStream, const GRGOutputFilter& filter, bool allowSimplify);
    friend class RenumberAndWriteVisitor;

    // Google-test unit tests that need private/protected access.
    FRIEND_TEST(GRG, TestTopoVisit);
    FRIEND_TEST(GRG, TestFrontier);
};

template <> inline MutationId GRG::nodeAndMutToType<MutationId>(const NodeAndMut& pair) { return pair.second; }

template <> inline GRG::MutAndNode GRG::nodeAndMutToType<GRG::MutAndNode>(const NodeAndMut& pair) {
    return {pair.second, INVALID_NODE_ID};
}

template <> inline MutationId GRG::nodeAndMutAndMissToType<MutationId>(const NodeMutMiss& triple) {
    return std::get<1>(triple);
}

template <> inline GRG::MutAndNode GRG::nodeAndMutAndMissToType<GRG::MutAndNode>(const NodeMutMiss& triple) {
    return {std::get<1>(triple), std::get<2>(triple)};
}

using GRGPtr = std::shared_ptr<GRG>;
using ConstGRGPtr = std::shared_ptr<const GRG>;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;
using ConstMutableGRGPtr = std::shared_ptr<const MutableGRG>;

template <>
inline std::tuple<NodeID, MutationId, NodeID>
GRG::NodeAndMutIterator<std::tuple<NodeID, MutationId, NodeID>>::VectIter3::operator*() {
    if (nullptr != m_pair) {
        const auto& item = (*m_pair).at(m_position);
        return {item.first, item.second, INVALID_NODE_ID};
    }
    return (*m_triple).at(m_position);
}

template <>
inline std::pair<NodeID, MutationId> GRG::NodeAndMutIterator<std::pair<NodeID, MutationId>>::VectIter3::operator*() {
    if (nullptr != m_triple) {
        const auto& item = (*m_triple).at(m_position);
        return {std::get<0>(item), std::get<1>(item)};
    }
    return (*m_pair).at(m_position);
}

/**
 * A MutableGRG can be changed by adding/removing nodes and edges.
 */
class MutableGRG : public GRG {
public:
    /**
     * Construct a GRG for a given number of samples.
     *
     * @param[in] The number of samples that will be used to construct the graph.
     * @param[in] (Optional) the initial capacity of the node vector. If you know in advance roughly
     *      how many nodes will be created this can improve performance.
     */
    explicit MutableGRG(NodeIDSizeT numSamples,
                        uint16_t ploidy,
                        bool phased = true,
                        NodeIDSizeT initialNodeCapacity = DEFAULT_NODE_CAPACITY,
                        bool useUpEdges = true)
        : GRG(numSamples, ploidy, phased),
          m_hasUpEdges(useUpEdges) {
        api_exc_check(numSamples > 0, "Must have at least one sample.");
        if (initialNodeCapacity < numSamples) {
            initialNodeCapacity = numSamples * 2;
        }
        this->m_nodes.reserve(initialNodeCapacity);
        this->makeNode(numSamples);
    }

    NodeIDSizeT numSamples() const override {
        if (m_orderedSampleList.empty()) {
            return GRG::numSamples();
        }
        return m_orderedSampleList.size();
    }

    bool isSample(const NodeID nodeId) const override {
        if (m_orderedSampleList.empty()) {
            return GRG::isSample(nodeId);
        }
        return m_sampleSet.find(nodeId) != m_sampleSet.end();
    }

    NodeIDList getSampleNodes() const override {
        if (m_orderedSampleList.empty()) {
            return std::move(GRG::getSampleNodes());
        }
        return m_orderedSampleList;
    }

    bool nodesAreOrdered() const override { return m_nodesAreOrdered && m_negativeNodes.empty(); }

    bool nodesAreTopo() const override { return m_nodesAreOrdered; }

    size_t numNodes() const override { return m_nodes.size(); }

    size_t numEdges() const override {
        size_t edgeCount = 0;
        for (const auto& node : m_nodes) {
            edgeCount += node->getDownEdges().size();
        }
        return edgeCount;
    }

    size_t numDownEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getDownEdges().size(); }

    size_t numUpEdges(NodeID nodeId) override {
        if (!m_hasUpEdges) {
            release_assert(m_nodes.at(nodeId)->getUpEdges().empty());
            return NO_UP_EDGES;
        }
        return m_nodes.at(nodeId)->getUpEdges().size();
    }

    NodeIDList getDownEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getDownEdges(); }

    NodeIDList getUpEdges(NodeID nodeId) override {
        api_exc_check(m_hasUpEdges, "No up edges were loaded/enabled");
        return m_nodes.at(nodeId)->getUpEdges();
    }

    bool hasUpEdges() const override { return m_hasUpEdges; }

    bool edgesAreOrdered() const override { return false; }

    NodeIDList getOrderedNodes(TraversalDirection direction, bool allowSort = false) override;

    /**
     * Create a new node in the GRG.
     *
     * This is the only valid way to construct a GRG node. When you create a GRG you specify the number
     * of sample nodes, and those nodes are created right away. Each call to `makeNode()` after that will
     * generate sequential-ID nodes. Use `connect()` to connect the newly created node to other nodes.
     *
     * @param[in] count The number of nodes to create.
     * @param[in] negative For "normal" nodes, the topological order is based on counting up with NodeIDs
     *      but for negative nodes this is not the case, they count (and build) down. If you set this true,
     *      the node will be added to a special list of prepended nodes so that the topological order can
     *      be maintained (via getOrderedNodes()... the integer iteration will no longer work). This does not
     *      return a negative value for NodeID (it is unsigned), and the NodeID increments just like a
     *      normal node. It is up to the user to use the NodeID appropriately - which means passing it as
     *      a negative value to connect() whenever edges are created with it.
     */
    NodeID makeNode(const size_t count = 1, bool negative = false) {
        const auto nextId = this->m_nodes.size();
        if (nextId + count > MAX_GRG_NODES) {
            std::stringstream ssErr;
            ssErr << "Cannot create more than " << MAX_GRG_NODES << " nodes in a GRG";
            throw ApiMisuseFailure(ssErr.str().c_str());
        }
        for (size_t i = 0; i < count; i++) {
            this->m_nodes.push_back(std::unique_ptr<GRGNode>(new GRGNode()));
            if (negative) {
                const NodeID nextNode = nextId + i;
                m_negativeNodes.emplace_back(nextNode);
#ifdef GRGL_CHECK_NEGATIVE
                m_isNegative.emplace(nextNode);
#endif
            }
        }
        if (this->m_nodes.size() > this->m_numSamples) {
            this->m_nodeData.allocNumCoals((this->m_nodes.size() - this->m_numSamples) + 1);
        }
        return nextId;
    }

    /**
     * Create a graph edge from one node (the source) to another node (the target).
     *
     * In GRG, "source node" always refers to the node "above" the other node. I.e., this function constructs
     * a down edge from source to target, and an up edge from target to source. In ARG parlance, the source
     * node is the parent and the target is the child.
     *
     * @param[in] srcId The ID of the source node.
     * @param[in] tgtId The ID of the target node.
     */
    void connect(SignedNodeID srcId, SignedNodeID tgtId);

    /**
     * Remove a graph edge between two nodes.
     *
     * The source and target nodes have the same meaning as in connect().
     *
     * @param[in] srcId The ID of the source node.
     * @param[in] tgtId The ID of the target node.
     */
    bool disconnect(NodeID srcId, NodeID tgtId);

    /**
     * Ensure that every Mutation in the GRG is unique: there is a single Mutation (with MutationId)
     * for each unique (position, ref allele, alt allele) combination. It does so by finding duplicate
     * Mutations, adding a new node that becomes a parent to them, and creating a new (single) mutation
     * representing them all. At the end of this method, the mutations will be unordered, and either
     * serializing the GRG to disk or calling sortMutations() will be necessary to get the new, correct
     * MutationIds.
     *
     * This operation is O(M), where M is the number of mutations.
     *
     * @return The number of nodes that were added to the graph.
     */
    NodeIDSizeT ensureUniqueMutations();

    /**
     * Merge one or more GRGs into this one. Only succeeds if all GRGs have the same number of
     * samples.
     *
     * This assumes that the GRGs were constructed from the same sampleset -- e.g., they
     * could be constructed from two subsets of the same sampleset (as long as both were
     * constructed with the same sample node numbering) or from a subset of mutations against
     * the same sampleset.
     *
     * The specified range of the resulting GRG will be (min(range of any input), max(range of any input)),
     * even if the provided GRGs do not span a contiguous region. It is up to the caller of this
     * method to ensure that either (1) the span is contiguous or (2) they adjust the specified range
     * appropriately afterwards.
     *
     * @param[in] otherGrgFiles The list of GRG filenames to load and merge.
     * @param[in] combineNodes Set to false to never combine nodes from the graphs.
     * @param[in] useSampleSets Set to true to use the slower, more RAM intensive version that tracks
     *      sets of samples beneath each graph node, and combines nodes that have the same sample set.
     *      The default algorithm just uses the (mapped) children to determine if two nodes can be
     *      combined, which combines fewer nodes overall, but also retains more hierarchy in the final
     *      graph.
     * @param[in] verbose Emit more information to stdout.
     * @param[in] ignoreRangeViolations Merge graphs that might have overlapping positions.
     * @param[in] positionAdjust Adjust each otherGrgFiles GRG by the given basepair position. I.e.,
     *      for each Mutation in otherGrgFiles[i], add positionAdjust[i] to its mutation position.
     */
    void merge(const std::list<std::string>& otherGrgFiles,
               bool combineNodes = true,
               bool useSampleSets = false,
               bool verbose = false,
               bool ignoreRangeViolations = false,
               std::vector<BpPosition> positionAdjust = {});

    /**
     * Change which nodes in the graph are the sample nodes. These nodes do not need to be leaf nodes,
     * but they should not be reachable from each other (this is not checked: the user must ensure).
     *
     * The order that the NodeIDs are passed in the input list is the order that the samples will be
     * numbered. So for example, sampleNodes[0] is sample 0, samplesNodes[1] is 1, etc. All methods
     * that deal with samples will treat them in this order, and when the GRG is written to disk the
     * nodes will be renumbered such that sampleNodes[0] becomes NodeID=0, and so on.
     *
     * @param[in] sampleNodes The ordered list of new sample nodes. This list length must be evenly
     *  divisible by getPloidy().
     */
    void setSamples(const NodeIDList& sampleNodes) {
        api_exc_check(sampleNodes.size() < numNodes(), "Too many sample nodes");
        api_exc_check(!sampleNodes.empty(), "Must have at least one sample");
        m_orderedSampleList.clear();
        m_sampleSet.clear();
        m_orderedSampleList.shrink_to_fit();
        m_orderedSampleList.reserve(sampleNodes.size());
        for (const NodeID node : sampleNodes) {
            api_exc_check(node < numNodes(), "Invalid sample node: " << node);
            m_orderedSampleList.emplace_back(node);
            m_sampleSet.insert(node);
        }
        api_exc_check(m_orderedSampleList.size() == m_sampleSet.size(), "Sample nodes must be unique");
    }

    /**
     * Retrieve a node by ID.
     *
     * @param[in] nodeId The ID of the node in question.
     */
    GRGNode& getNode(const NodeID nodeId) const {
        api_exc_check(nodeId < this->m_nodes.size(), "Invalid NodeID: " << nodeId);
        return *this->m_nodes[nodeId];
    }

    std::vector<NodeIDSizeT> topologicalSort(TraversalDirection direction) override;

    void visitTopo(GRGVisitor& visitor,
                   TraversalDirection direction,
                   const NodeIDList& seedList,
                   const std::vector<NodeIDSizeT>* sortOrder = nullptr) override;

    /**
     * Compact the edge vectors in the GRG. If done infrequently won't affect the amortized edge addition
     * cost but will reduce overall RAM usage by a non-trivial amount.
     */
    void compact(NodeID nodeId = INVALID_NODE_ID);

    bool isMutable() const override { return true; }

    bool samplesAreOrdered() const override { return m_orderedSampleList.empty(); }

private:
    // Empty, or the map from NodeIDs that correspond to the samples to the "sample ID"
    // (e.g., the NodeID that it will get when/if the GRG is written to disk).
    NodeIDSet m_sampleSet;
    NodeIDList m_orderedSampleList;

    // The list of nodes. The node's position in this vector must match its ID.
    std::vector<std::unique_ptr<GRGNode>> m_nodes;

    friend MutableGRGPtr readMutableGrg(IFSPointer&, bool);

    // Negative nodes in the order they were added to the GRG. See makeNode() for details.
    std::vector<NodeID> m_negativeNodes;
    // In debug builds, sanity check the use of negative nodes. At large scale this can have
    // non-trivial overhead, so only in debug builds.
#ifdef GRGL_CHECK_NEGATIVE
    NodeIDSet m_isNegative;
#endif

    // True if the nodeId order matches the bottom-up topological order.
    bool m_nodesAreOrdered{true};
    // True if we are tracking up edges -- not all analyses need them.
    bool m_hasUpEdges;
};

class CSRGRG : public GRG {
public:
    explicit CSRGRG(size_t numSamples, size_t nodeCount, uint16_t ploidy, bool loadUpEdges = true, bool phased = true)
        : GRG(numSamples, ploidy, phased),
          m_downEdges(nodeCount),
          m_upEdges(loadUpEdges ? nodeCount : 0) {
        api_exc_check(numSamples > 0, "Must have at least one sample.");
    }

    bool nodesAreOrdered() const override { return true; }

    bool nodesAreTopo() const override { return true; }

    size_t numNodes() const override { return m_downEdges.numNodes(); }

    size_t numEdges() const override { return m_downEdges.numValues(); }

    size_t numDownEdges(NodeID nodeId) override { return m_downEdges.numValuesAt(nodeId); }

    bool hasUpEdges() const override { return m_upEdges.numNodes() > 0; }

    bool edgesAreOrdered() const override { return true; }

    size_t numUpEdges(NodeID nodeId) override {
        if (!hasUpEdges()) {
            return NO_UP_EDGES;
        }
        return m_upEdges.numValuesAt(nodeId);
    }

    NodeIDList getDownEdges(NodeID nodeId) override {
        NodeIDList result;
        m_downEdges.getData(nodeId, result);
        return std::move(result);
    }

    NodeIDList getUpEdges(NodeID nodeId) override {
        if (!hasUpEdges()) {
            throw ApiMisuseFailure("No up edges were loaded/created");
        }
        NodeIDList result;
        m_upEdges.getData(nodeId, result);
        return std::move(result);
    }

    // For CSRGRG, you can just iterate the nodes in order from 0...(numNodes-1). This function
    // only exists for the Python API, and also to have an agnostic function if you don't
    // know (or care) what kind of GRG you have.
    NodeIDList getOrderedNodes(TraversalDirection direction, bool allowSort = false) override {
        (void)allowSort;
        NodeIDList result(numNodes());
        if (direction == grgl::TraversalDirection::DIRECTION_DOWN) {
            std::iota(result.rbegin(), result.rend(), 0);
        } else {
            std::iota(result.begin(), result.end(), 0);
        }
        return std::move(result);
    }

    std::vector<NodeIDSizeT> topologicalSort(TraversalDirection direction) override;

    void visitTopo(GRGVisitor& visitor,
                   TraversalDirection direction,
                   const NodeIDList& seedList,
                   const std::vector<NodeIDSizeT>* sortOrder = nullptr) override;

    EagerCSREdges32 m_downEdges;
    EagerCSREdges32 m_upEdges;

private:
    friend GRGPtr readImmutableGrg(IFSPointer& inStream, bool loadUpEdges, bool loadDownEdges);
};

using CSRGRGPtr = std::shared_ptr<CSRGRG>;
using ConstCSRGRGPtr = std::shared_ptr<const CSRGRG>;

// This header is just separated out to try to keep grg.h well organized. It contains templated
// functions/classes that are needed by matrix multiplication.
#include "internal_mult.h"

template <typename NodeValueType, bool useBitVector>
inline void matmulSumChildrenDown(GRG* grg,
                                  const NodeID nodeId,
                                  std::vector<NodeValueType>& nodeValues,
                                  const size_t effectiveInputRows) {
    const size_t base = nodeId * effectiveInputRows;
    for (NodeID childId : grg->getDownEdges(nodeId)) {
        const size_t cbase = childId * effectiveInputRows;
        vectorAdd<NodeValueType, useBitVector>(nodeValues.data(), nodeValues.data(), cbase, base, effectiveInputRows);
    }
}

template <typename NodeValueType, bool useBitVector>
inline void matmulSumChildrenUp(GRG* grg,
                                const NodeID nodeId,
                                std::vector<NodeValueType>& nodeValues,
                                const size_t effectiveInputRows) {
    const size_t base = nodeId * effectiveInputRows;
    for (NodeID childId : grg->getDownEdges(nodeId)) {
        const size_t cbase = childId * effectiveInputRows;
        vectorAdd<NodeValueType, useBitVector>(nodeValues.data(), nodeValues.data(), base, cbase, effectiveInputRows);
    }
}

// Notes on parameter constraints/assumptions (_CALLERS_ must ensure these):
// * initMatrix and nodeInit match -- the dimensions of the former are defined by the latter
// * missMatrix, if provided, must be a vector of length inputCols (C). Unlike inputMatrix which
//   is (k x C), missMatrix is (1 x C). I.e., you can't provide different missing values for each
//   row (it is assumed to be constant).
template <typename IOType, typename NodeValueType, bool useBitVector>
void GRG::matrixMultiplication(const IOType* inputMatrix,
                               size_t inputCols,
                               size_t inputRows,
                               TraversalDirection direction,
                               IOType* outputMatrix,
                               size_t outputSize,
                               bool emitAllNodes,
                               bool byIndividual,
                               const IOType* initMatrix,
                               NodeInitEnum nodeInit,
                               IOType* missMatrix) {
    release_assert(inputCols > 0);
    release_assert(inputRows > 0);
    const size_t outputCols = outputSize / inputRows;
    validateMatMulInputs(this, inputCols, inputRows, direction, outputSize, emitAllNodes, byIndividual, outputCols);
    // When we do bitvector calculations, we must have the number of input rows be a multiple of the
    // element size the bitvector is using.
    const size_t effectiveInputRows =
        useBitVector ? roundUpToMultiple(inputRows, sizeof(NodeValueType) * 8) : inputRows;

    // The node value storage stores the different row values consecutively (so
    // similar to column-major order), in contrast to the input/output matrices
    // which are row-major.
    const size_t nodesPerElem = useBitVector ? sizeof(NodeValueType) * 8 : 1;
    release_assert(effectiveInputRows % nodesPerElem == 0);
    std::vector<NodeValueType> nodeValues(numNodes() * effectiveInputRows / nodesPerElem);

    api_exc_check(nodeInit == NIE_ZERO || !useBitVector, "Node initialization and bit vector type cannot be mixed");
    switch (nodeInit) {
    case NIE_XTX:
        for (NodeID i = 0; i < numNodes(); i++) {
            const size_t base = i * effectiveInputRows;
            const NodeIDSizeT coalCount = this->getNumIndividualCoals(i);
            api_exc_check(coalCount != COAL_COUNT_NOT_SET,
                          "Coalescent counts not available for this GRG (init=XTX won't work)");
            matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                nodeValues.data(), base, coalCount * 2, effectiveInputRows);
        }
        break;
    case NIE_VECTOR:
        for (NodeID i = 0; i < numNodes(); i++) {
            const size_t base = i * effectiveInputRows;
            for (size_t row = 0; row < inputRows; row++) {
                matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                    nodeValues.data(), base + row, initMatrix, row);
            }
        }
        break;
    case NIE_MATRIX:
        for (NodeID i = 0; i < numNodes(); i++) {
            const size_t base = i * effectiveInputRows;
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * numNodes();
                matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                    nodeValues.data(), base + row, initMatrix, rowStart + i);
            }
        }
        break;
    case NIE_ZERO:
    default: break;
    }

    // Downward, we are calculating "how do the mutations impact the samples?"
    if (direction == DIRECTION_DOWN) {
        for (const auto& tuple : this->getNodesAndMutations<GRG::NodeMutMiss>()) {
            const NodeID& nodeId = std::get<0>(tuple);
            const MutationId& mutId = std::get<1>(tuple);
            const NodeID& missingnessNode = std::get<2>(tuple);
            assert(mutId < inputCols);
            if (nodeId != INVALID_NODE_ID) {
                const size_t base = nodeId * effectiveInputRows;
                for (size_t row = 0; row < inputRows; row++) {
                    const size_t rowStart = row * inputCols;
                    matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                        nodeValues.data(), base + row, inputMatrix, rowStart + mutId);
                }
            }
            if (missMatrix != nullptr && missingnessNode != INVALID_NODE_ID) {
                const size_t base = missingnessNode * effectiveInputRows;
                for (size_t row = 0; row < inputRows; row++) {
                    const size_t rowStart = row * inputCols;
                    matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                        nodeValues.data(), base + row, missMatrix, rowStart + mutId);
                }
            }
        }
        if (this->nodesAreOrdered()) {
            for (NodeID i = numNodes(); i > 0; i--) {
                matmulSumChildrenDown<NodeValueType, useBitVector>(this, i - 1, nodeValues, effectiveInputRows);
            }
        } else if (this->nodesAreTopo()) {
            for (NodeID nodeId : this->getOrderedNodes(TraversalDirection::DIRECTION_DOWN)) {
                matmulSumChildrenDown<NodeValueType, useBitVector>(this, nodeId, nodeValues, effectiveInputRows);
            }
        } else {
            ValueSumVisitor<NodeValueType> valueSumVisitor(nodeValues, effectiveInputRows);
            this->visitDfs(valueSumVisitor, DIRECTION_UP, getSampleNodes());
        }
        if (!emitAllNodes) {
            const NodeIDList& samples = getSampleNodes();
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * outputCols;
                for (NodeID sampleIdx = 0; sampleIdx < samples.size(); sampleIdx++) {
                    const NodeID sampleNodeId = samples[sampleIdx];
                    const size_t srcBase = sampleNodeId * effectiveInputRows; // Node space
                    assert(rowStart + sampleIdx < outputSize);
                    const size_t generalSample = byIndividual ? (sampleIdx / m_ploidy) : sampleIdx; // IO space
                    matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                        outputMatrix, rowStart + generalSample, nodeValues.data(), srcBase + row);
                }
            }
        }
        // Upward, we are calculating "how do the samples impact the mutations?"
    } else {
        // Copy from the input matrix to the sample node values. The destination is per-sample-node,
        // and the source is per-sample-index (ID of 0...(N-1) where N is number of haplotypes).
        const NodeIDList& samples = getSampleNodes();
        for (NodeID sampleIdx = 0; sampleIdx < samples.size(); sampleIdx++) {
            const NodeID sampleNodeId = samples[sampleIdx];
            assert(sampleIdx < inputCols);
            const size_t destBase = sampleNodeId * effectiveInputRows;
            const size_t generalSample = byIndividual ? (sampleIdx / m_ploidy) : sampleIdx;
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * inputCols;
                matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                    nodeValues.data(), destBase + row, inputMatrix, rowStart + generalSample);
            }
        }
        if (this->nodesAreOrdered()) {
            const NodeID start = this->samplesAreOrdered() ? numSamples() : 0;
            for (NodeID nodeId = start; nodeId < numNodes(); nodeId++) {
                matmulSumChildrenUp<NodeValueType, useBitVector>(this, nodeId, nodeValues, effectiveInputRows);
            }
        } else if (this->nodesAreTopo()) {
            for (NodeID nodeId : this->getOrderedNodes(TraversalDirection::DIRECTION_UP)) {
                matmulSumChildrenUp<NodeValueType, useBitVector>(this, nodeId, nodeValues, effectiveInputRows);
            }
        } else {
            ValueSumVisitor<NodeValueType> valueSumVisitor(nodeValues, effectiveInputRows);
            this->visitDfs(valueSumVisitor, DIRECTION_DOWN, getRootNodes());
        }
        if (!emitAllNodes) {
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * outputCols;
                for (const auto& triple : this->getNodesAndMutations<GRG::MutNodeMiss>()) {
                    const NodeID& nodeId = std::get<0>(triple);
                    const MutationId& mutId = std::get<1>(triple);
                    assert(rowStart + mutId < outputSize);
                    if (nodeId != INVALID_NODE_ID) {
                        const size_t base = nodeId * effectiveInputRows;
                        matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                            outputMatrix, rowStart + mutId, nodeValues.data(), base + row);
                    }
                    const NodeID& missingnessNode = std::get<2>(triple);
                    if (missMatrix != nullptr && missingnessNode != INVALID_NODE_ID) {
                        const size_t base = missingnessNode * effectiveInputRows;
                        // Add the missing-data value for this mutation to the output vector position
                        // associated with the mutation.
                        matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                            missMatrix, rowStart + mutId, nodeValues.data(), base + row);
                    }
                }
            }
        }
    }
    if (emitAllNodes) {
        for (size_t row = 0; row < inputRows; row++) {
            const size_t rowStart = row * outputCols;
            for (NodeID nodeId = 0; nodeId < numNodes(); nodeId++) {
                const size_t base = nodeId * effectiveInputRows;
                matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                    outputMatrix, rowStart + nodeId, nodeValues.data(), base + row);
            }
        }
    }
}

} // namespace grgl

#endif /* GRG_H */
