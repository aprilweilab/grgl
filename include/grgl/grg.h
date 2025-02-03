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
#include <vector>

#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/file_vector.h"
#include "grgnode.h"
#include "mutation.h"
#include "node_data.h"
#include "visitor.h"

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

    TraversalDirection direction;
    NodeIDList seedList;
    std::pair<BpPosition, BpPosition> bpRange;
};

/**
 * Abstract GRG base class.
 */
class GRG : public std::enable_shared_from_this<GRG> {
public:
    enum {
        DEFAULT_NODE_CAPACITY = 1024,
    };

    explicit GRG(size_t numSamples, uint16_t ploidy)
        : m_numSamples(numSamples),
          m_ploidy(ploidy) {
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
     * The first numSamples node IDs are reserved for samples, so you can also just use
     * `nodeId < grg.numSamples()`.
     *
     * @param[in] nodeId The node ID to check.
     */
    bool isSample(const NodeID nodeId) const { return nodeId < this->m_numSamples; }

    size_t numSamples() const { return m_numSamples; }

    size_t numIndividuals() const { return m_numSamples / m_ploidy; }

    /**
     * How many haploid samples are there per individual?
     *
     * @return The ploidy, usually 1 or 2. Individual coalescence support only works when ploidy==2.
     */
    uint16_t getPloidy() const { return m_ploidy; }

    /**
     * Returns true if nodes are ordered in bottom-up topological order.
     */
    virtual bool nodesAreOrdered() const = 0;

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
     * True if this graph stores up edges, false if only down edges are stored.
     */
    virtual bool hasUpEdges() const = 0;

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
            return {m_mutations[0].getPosition(), m_mutations[lastMutId].getPosition()};
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
     * Get the NodeIDList for all samples in the graph.
     */
    NodeIDList getSampleNodes() const {
        NodeIDList result;
        for (grgl::NodeID i = 0; i < this->numSamples(); i++) {
            result.push_back(i);
        }
        return std::move(result);
    }

    /**
     * Get the NodeIDList for all roots in the graph.
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

    const std::vector<std::pair<NodeID, MutationId>>& getNodeMutationPairs() {
        if (!m_mutIdsByNodeIdSorted) {
            sortMutIdsByNodeID();
        }
        return m_mutIdsByNodeId;
    }

    const Mutation& getMutationById(MutationId mutId) { return m_mutations.cref(mutId); }

    void setMutationById(MutationId mutId, Mutation mutation) { m_mutations.atRef(mutId) = std::move(mutation); }

    std::vector<MutationId> getUnmappedMutations() { return getMutationsForNode(INVALID_NODE_ID); }

    std::vector<MutationId> getMutationsForNode(const NodeID nodeId) {
        if (!m_mutIdsByNodeIdSorted) {
            sortMutIdsByNodeID();
        }
        std::vector<MutationId> result;
        std::pair<NodeID, MutationId> query = {nodeId, 0};
        auto mutIdIt = std::lower_bound(m_mutIdsByNodeId.begin(), m_mutIdsByNodeId.end(), query);
        while (mutIdIt != m_mutIdsByNodeId.end() && mutIdIt->first == nodeId) {
            result.push_back(mutIdIt->second);
            mutIdIt++;
        }
        return std::move(result);
    }

    bool nodeHasMutations(const NodeID nodeId) {
        if (!m_mutIdsByNodeIdSorted) {
            sortMutIdsByNodeID();
        }
        std::pair<NodeID, MutationId> query = {nodeId, 0};
        auto mutIdIt = std::lower_bound(m_mutIdsByNodeId.begin(), m_mutIdsByNodeId.end(), query);
        return mutIdIt != m_mutIdsByNodeId.end() && (mutIdIt->first == nodeId);
    }

    /**
     * Get pairs of mutation IDs and node IDs, ordered by the mutation position + allele (ascending).
     *
     * @return A vector of pairs, MutationID and NodeID (in that order).
     */
    std::vector<std::pair<MutationId, NodeID>> getMutationsToNodeOrdered();

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
     * down edges.
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
     * @return The MutationId for the newly added Mutation.
     */
    MutationId addMutation(const Mutation& mutation, NodeID nodeId);

    // Internal method, used by Python API and the main API below.
    void dotProduct(const double* inputData,
                    size_t inputLength,
                    TraversalDirection direction,
                    double* outputData,
                    size_t outputLength);

    /**
     * Compute one of two possible dot products across the entire graph. The input vector \f$V\f$ can be
     * either \f$1 \times N\f$ (\f$N\f$ is number of samples) or \f$1 \times M\f$ (\f$M\f$ is number of
     * mutations). The given direction determines which input vector is expected. Let \f$X\f$ be the
     * \f$N \times M\f$ genotype matrix.
     * For an \f$1 \times N\f$ input \f$V\f$, the product performed is \f$V \cdot X\f$ which gives a
     * \f$1 \times M\f$ result. I.e., the input vector is a value per sample and the output vector is
     * a value per mutation.
     * For an \f$1 \times M\f$ input \f$V\f$, the product performed is \f$V \cdot X^T\f$ which gives a
     * \f$1 \times N\f$ result. I.e., the input vector is a value per mutation and the output vector is
     * a value per sample.
     *
     * Dot product in the graph works by seeding the input nodes (samples or mutations) with the corresponding
     * values from the input vector and then traversing the graph in the relevant direction (up or down). The
     * ancestor/descendant values are summed at each node, until the terminal nodes (mutations or samples) are
     * reached.
     *
     * @param[in] inputData The \f$1 \times N\f$ or \f$1 \times M\f$ input vector. For samples, the \f$i\f$th
     *      entry in the vector is the value for sample \f$i\f$. Similarly for mutations, the ordering follows
     *      the MutationId ordering.
     * @param[in] direction Use TraversalDirection::DIRECTION_DOWN for mutations as input, and
     *      TraversalDirection::DIRECTION_UP for samples as input.
     * @return The resulting vector as described above.
     */
    std::vector<double> dotProduct(const std::vector<double>& inputData, TraversalDirection direction) {
        const size_t outSize = (direction == TraversalDirection::DIRECTION_DOWN) ? numSamples() : numMutations();
        std::vector<double> result(outSize);
        dotProduct(inputData.data(), inputData.size(), direction, result.data(), outSize);
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
     * Set the number of individuals that coalesce at this node.
     * @param[in] nodeId The node to modify.
     * @param[in] coals The number of individuals that coalesce, between 0...numSamples()/ploidy.
     */
    void setNumIndividualCoals(NodeID nodeId, NodeIDSizeT coals) {
        assert(coals == COAL_COUNT_NOT_SET || coals <= numIndividuals());
        m_nodeData.setNumCoals(numSamples(), nodeId, coals);
    }

protected:
    void visitTopoNodeOrdered(GRGVisitor& visitor, TraversalDirection direction, const NodeIDList& seedList);

    // m_mutIdsByNodeId gets appended to and then sorted periodically or as needed, using this
    // function. We rarely need to lookup Mutations by NodeID while modifying a GRG, so this ends
    // up being about 3x smaller than the previous multimap<> and many times faster.
    void sortMutIdsByNodeID() {
        std::sort(m_mutIdsByNodeId.begin(), m_mutIdsByNodeId.end());
        m_mutIdsByNodeIdSorted = true;
    }

    // The position is the mutation ID for both of these vectors.
    EagerFileVector<Mutation> m_mutations;
    // This vector is MutationIDs ordered by their associated NodeID, ascending.
    std::vector<std::pair<NodeID, MutationId>> m_mutIdsByNodeId;
    bool m_mutIdsByNodeIdSorted{false};

    // (Optional) list of population descriptions. The position corresponds to the population
    // ID, which can be used to tag nodes.
    std::vector<std::string> m_populations;

    // The range of base-pair positions covered by this GRG, according to the user/tool that created it.
    std::pair<BpPosition, BpPosition> m_specifiedRange{};

    // Node data.
    NodeDataContainer m_nodeData;

    const size_t m_numSamples;
    const uint16_t m_ploidy;

    // True if the mutationId order matches the (position, allele) order.
    bool m_mutsAreOrdered{false};

    friend void readGrgCommon(const GRGFileHeader& header, const GRGPtr& grg, IFSPointer& inStream);
};

using GRGPtr = std::shared_ptr<GRG>;
using ConstGRGPtr = std::shared_ptr<const GRG>;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;
using ConstMutableGRGPtr = std::shared_ptr<const MutableGRG>;

class MutableGRG : public GRG {
public:
    /**
     * Construct a GRG for a given number of samples.
     *
     * @param[in] The number of samples that will be used to construct the graph.
     * @param[in] (Optional) the initial capacity of the node vector. If you know in advance roughly
     *      how many nodes will be created this can improve performance.
     */
    explicit MutableGRG(size_t numSamples, uint16_t ploidy, size_t initialNodeCapacity = DEFAULT_NODE_CAPACITY)
        : GRG(numSamples, ploidy) {
        if (initialNodeCapacity < numSamples) {
            initialNodeCapacity = numSamples * 2;
        }
        this->m_nodes.reserve(initialNodeCapacity);
        this->makeNode(numSamples, true);
    }

    bool nodesAreOrdered() const override { return m_nodesAreOrdered; }

    size_t numNodes() const override { return m_nodes.size(); }

    size_t numEdges() const override {
        size_t edgeCount = 0;
        for (const auto& node : m_nodes) {
            edgeCount += node->getDownEdges().size();
        }
        return edgeCount;
    }

    size_t numDownEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getDownEdges().size(); }

    size_t numUpEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getUpEdges().size(); }

    NodeIDList getDownEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getDownEdges(); }

    NodeIDList getUpEdges(NodeID nodeId) override { return m_nodes.at(nodeId)->getUpEdges(); }

    bool hasUpEdges() const override { return true; }

    /**
     * Create a new node in the GRG.
     *
     * This is the only valid way to construct a GRG node. When you create a GRG you specify the number
     * of sample nodes, and those nodes are created right away. Each call to `makeNode()` after that will
     * generate sequential-ID nodes. Use `connect()` to connect the newly created node to other nodes.
     *
     * @param[in] count The number of nodes to create.
     * @param[in] forceOrdered If true, do not "break" the topological ordering property of the graph if
     *      it already exists. Use only when you are certain that newly added nodes (and their edges) will
     *      main topological order.
     */
    NodeID makeNode(const size_t count = 1, bool forceOrdered = false) {
        const auto nextId = this->m_nodes.size();
        for (size_t i = 0; i < count; i++) {
            this->m_nodes.push_back(std::make_shared<GRGNode>());
        }
        if (!forceOrdered) {
            this->m_nodesAreOrdered = false;
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
    void connect(NodeID srcId, NodeID tgtId);

    /**
     * Remove a graph edge between two nodes.
     *
     * The source and target nodes have the same meaning as in connect().
     *
     * @param[in] srcId The ID of the source node.
     * @param[in] tgtId The ID of the target node.
     */
    void disconnect(NodeID srcId, NodeID tgtId);

    /**
     * Merge another GRG into this one. Only succeeds if both GRGs have the same number of
     * samples.
     *
     * This assumes that the GRGs were constructed from the same sampleset -- e.g., they
     * could be constructed from two subsets of the same sampleset (as long as both were
     * constructed with the same sample node numbering) or from a subset of mutations against
     * the same sampleset.
     *
     * @param otherGrgFiles The list of GRG filenames to load and merge.
     */
    void merge(const std::list<std::string>& otherGrgFiles, bool combineNodes = true);

    /**
     * Retrieve a node by ID.
     *
     * @param[in] nodeId The ID of the node in question.
     */
    GRGNodePtr getNode(const NodeID nodeId) const {
        const NodeID checkId = removeMarks(nodeId);
        return this->m_nodes.at(checkId);
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

private:
    // The list of nodes. The node's position in this vector must match its ID.
    std::vector<GRGNodePtr> m_nodes;

    friend MutableGRGPtr readMutableGrg(IFSPointer& inStream);

    // True if the nodeId order matches the bottom-up topological order.
    bool m_nodesAreOrdered{true};
};

// Eager CSR sorted 32-bit edges, encoded.
using EagerCSREdges32 = CSRStorageImm<EagerFileVector, uint32_t, true, true>;
// Lazy CSR sorted 32-bit edges, encoded.
using LazyCSREdges32 = CSRStorageImm<LazyFileVector, uint32_t, true, true>;

class CSRGRG : public GRG {
public:
    explicit CSRGRG(size_t numSamples, size_t nodeCount, uint16_t ploidy, bool loadUpEdges = true)
        : GRG(numSamples, ploidy),
          m_downEdges(nodeCount),
          m_upEdges(loadUpEdges ? nodeCount : 0) {}

    bool nodesAreOrdered() const override { return true; }

    size_t numNodes() const override { return m_downEdges.numNodes(); }

    size_t numEdges() const override { return m_downEdges.numValues(); }

    size_t numDownEdges(NodeID nodeId) override { return m_downEdges.numValuesAt(nodeId); }

    bool hasUpEdges() const override { return m_upEdges.numNodes() > 0; }

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

} // namespace grgl

#endif /* GRG_H */
