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
    NodeIDSizeT setKeepSamples(const grgl::GRGPtr& grg, grgl::NodeIDList sampleIDList, bool warn);
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
    std::ostream& m_outStream;
    // Maps old NodeID to new NodeID
    std::vector<NodeIDSizeT> m_nodeIdMap;
    NodeIDSizeT m_newSampleCount{};
    NodeIDSizeT m_nodeCounter{};
    std::unique_ptr<CSRStringTable> m_filteredIndivIds;
    EdgeSizeT m_edgeCounter{};
    size_t m_ploidy{};
    bool m_allowSimplify;
};

}; // namespace grgl

#endif /* GRG_SERIALIZE_H */
