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
#ifndef GRG_SERIALIZE_H
#define GRG_SERIALIZE_H

#include "grgl/csr_storage.h"
#include "grgl/file_vector.h"
#include "grgl/grgnode.h"
#include "grgl/node_data.h"
#include "grgl/visitor.h"
#include <iosfwd>
#include <memory>
#include <stdexcept>

namespace grgl {

/**
 * Exception thrown when a GRG fails to serialize/deserialize.
 */
class SerializationFailure : public std::runtime_error {
public:
    explicit SerializationFailure(char const* const message)
        : std::runtime_error(message) {}
};

class GRG;
class MutableGRG;
using GRGPtr = std::shared_ptr<GRG>;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;
using IFSPointer = std::shared_ptr<std::istream>;

/**
 * Serialize the GRG to the given outstream.
 *
 * @param[in] grg The GRG to be serialized.
 * @param[in] out The (binary) output stream.
 */
std::pair<NodeIDSizeT, EdgeSizeT> writeGrg(const GRGPtr& grg, std::ostream& out, bool allowSimplify = true);

class GRGOutputFilter;

std::pair<NodeIDSizeT, EdgeSizeT> simplifyAndSerialize(const GRGPtr& grg,
                                                       std::ostream& outStream,
                                                       const GRGOutputFilter& filter,
                                                       bool allowSimplify = true);

/**
 * Deserialize the GRG from the given input stream.
 *
 * @param[in] inStream The (binary) input stream.
 */
MutableGRGPtr readMutableGrg(IFSPointer& inStream, bool loadUpEdges = true);
GRGPtr readImmutableGrg(IFSPointer& inStream, bool loadUpEdges = true);

// This class is only exposed for unit testing. It is not part of the GRGL API.
class RenumberAndWriteVisitor : public GRGVisitor {
public:
    explicit RenumberAndWriteVisitor(std::ostream& outStream, NodeIDSizeT numNodes, bool allowSimplify);
    bool getChildren(const grgl::GRGPtr& grg, const grgl::NodeID nodeId, NodeIDList& result, NodeIDSizeT& parentCoals);
    bool shouldKeep(const grgl::GRGPtr& grg, const grgl::NodeID nodeId) const;
    NodeIDSizeT setKeepSamples(const grgl::GRGPtr& grg, const grgl::NodeIDList& sampleIdList, bool warn);
    NodeIDList setKeepMutations(const grgl::GRGPtr& grg, const grgl::NodeIDList& mutationIDList);
    bool keepMutation(const MutationId mutId);
    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override;
    void finalize();
    NodeIDSizeT getNumNodes() const;
    EdgeSizeT getNumEdges() const;
    NodeID getNewID(const NodeID nodeId) const;
    NodeID getOldID(const NodeID nodeId) const;

    // Maps new NodeID to old NodeID (so the order is new NodeID acscending)
    std::vector<NodeIDSizeT> m_revIdMap;
    std::vector<bool> m_keepBeneath;
    std::vector<bool> m_keepMutations;
    // True if the node is a missingness node, false otherwise.
    std::vector<bool> m_isMissingnessNode;

    // Updated node data
    NodeDataContainer m_newNodeData;

    EagerCSREdges32 m_edgeCSR;
    size_t m_bytesWritten{};

    bool hasIndividualIds(const GRGPtr& grg);
    CSRStringTable& getIndividualIds(const GRGPtr& grg);
    size_t getPloidy(const GRGPtr& grg) const;

private:
    NodeIDSizeT getAdjustedCoalCount(const GRGPtr& grg, NodeID nodeId);

    std::ostream& m_outStream;
    // Maps old NodeID to new NodeID
    std::vector<NodeIDSizeT> m_nodeIdMap;
    // Adjustments for coalescences for each node. Only used if we are downsampling (by individual)
    // when serializing the GRG. When m_coalAdjustIsSubtract then these values represent the coalescences
    // of the removed samples, and we subtract them from the existing values. Otherwise these are the
    // new node coalescences that replace existing values.
    std::vector<NodeIDSizeT> m_coalAdjust;
    bool m_coalAdjustIsSubtract{};

    NodeIDSizeT m_newSampleCount{};
    NodeIDSizeT m_nodeCounter{};
    std::unique_ptr<CSRStringTable> m_filteredIndivIds;
    EdgeSizeT m_edgeCounter{};
    size_t m_ploidy{};
    bool m_allowSimplify;
};

/**
 * Load the GRG so that it can be modified, by adding/removing nodes and edges.
 *
 * @param[in] filename The file to load.
 * @param[in] loadUpEdges Whether to load the "up" edges (in addition to the down edges). Default: true.
 * @return A shared_ptr to the MutableGRG object.
 */
MutableGRGPtr loadMutableGRG(const std::string& filename, const bool loadUpEdges = true);

/**
 * Load the GRG read-only. Edges and nodes cannot be changed, but the Mutations and other auxiliary
 * information still can be.
 *
 * @param[in] filename The file to load.
 * @param[in] loadUpEdges Whether to load the "up" edges (in addition to the down edges). Default: false.
 * @return A shared_ptr to the GRG object.
 */
GRGPtr loadImmutableGRG(const std::string& filename, bool loadUpEdges = false);

std::pair<NodeIDSizeT, EdgeSizeT> saveGRG(const GRGPtr& theGRG, const std::string& filename, bool allowSimplify = true);

/**
 * Save a subset of the GRG, either by Mutations (downward traversal) or Samples (upward traversal).
 *
 * @param[in] theGRG The GRG object to save.
 * @param[in] filename The output filename to use.
 * @param[in] direction Subset by mutations (TraversalDirection::DIRECTION_DOWN) or samples
 * (TraversalDirection::DIRECTION_UP)
 * @param[in] seedList The list of MutationIDs or SampleIDs (NodeIDs of samples) to keep. For samples, the order of this
 * list matters - it defines the order of the samples in the new GRG. For example, [0, 1, 2, 3] keeps the first four
 * samples in their original order, but [3, 2, 1, 0] keeps the same samples but reverses their order in the new GRG.
 * This means that the SampleIDs will change as 0->3, 1->2, 2->1, 3->0. The IndividualIDs mapped to individuals will be
 * properly maintained if present.
 * @param[in] byRange Just metadata about the range of base-pairs that this subset will cover. It does not have any
 * effect on the actual filtering option, IT ONLY CHANGES METADATA.
 */
bool saveGRGSubset(const GRGPtr& theGRG,
                   const std::string& filename,
                   const TraversalDirection direction,
                   const NodeIDList& seedList,
                   std::pair<BpPosition, BpPosition> bpRange);

}; // namespace grgl

#endif /* GRG_SERIALIZE_H */
