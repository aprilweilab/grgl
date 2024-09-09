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
 * should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "grgl/serialize.h"
#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "util.h"
#include "varint.h"

#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>

constexpr uint64_t GRG_MAGIC = 0xE9366C64DDC8C5B0;
constexpr uint32_t GRG_MAJOR_VERSION = 4;
constexpr uint32_t GRG_MINOR_VERSION = 1;

// The file header stores the nodeID size. This special value means the nodeID size
// should use the variable-sized int encoding (see varint.h)
#define USE_VARINT (0xFF)

#define assert_deserialization(condition, msg)                                                                         \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::stringstream errMsg;                                                                                  \
            errMsg << "Failed (de)serialization: " << (msg);                                                           \
            throw SerializationFailure(msg);                                                                           \
        }                                                                                                              \
    } while (0)

namespace grgl {

constexpr uint64_t HAS_INDIVIDUAL_COALS = 0x2;
constexpr uint64_t MUTATION_IDS_ARE_ORDERED = 0x4;

#pragma pack(push, 1)
struct GRGFileHeader {
    uint64_t magic;
    uint32_t versionMajor;
    uint32_t versionMinor;
    uint64_t idSize;
    uint64_t sampleCount;
    uint64_t mutationCount;
    uint64_t nodeCount;
    uint64_t edgeCount;
    uint16_t populationCount;
    uint16_t ploidy;
    uint32_t unused32;
    uint64_t flags;
    uint64_t unused[7];
};
#pragma pack(pop)

static size_t getIdBytes(const GRGPtr& grg) {
    if (grg->numNodes() < std::numeric_limits<uint8_t>::max()) {
        return sizeof(uint8_t);
    } else if (grg->numNodes() < std::numeric_limits<uint16_t>::max()) {
        return sizeof(uint16_t);
    } else if (grg->numNodes() < std::numeric_limits<uint32_t>::max()) {
        return sizeof(uint32_t);
    }
    return sizeof(uint64_t);
}

template <typename T> static inline void writeScalar(T intValue, std::ostream& outStream) {
    outStream.write(reinterpret_cast<const char*>(&intValue), sizeof(intValue));
}

static inline void writeString(bool useVarInt, const std::string& value, std::ostream& outStream) {
    if (useVarInt) {
        writeVarInt(value.size(), outStream);
    } else {
        writeScalar<uint64_t>(value.size(), outStream);
    }
    outStream.write(value.c_str(), static_cast<std::streamsize>(value.size()));
}

static void writeNodeID(const NodeID nodeId, size_t idSize, std::ostream& outStream) {
    switch (idSize) {
    case sizeof(uint8_t): writeScalar<uint8_t>(nodeId, outStream); break;
    case sizeof(uint16_t): writeScalar<uint16_t>(nodeId, outStream); break;
    case sizeof(uint32_t): writeScalar<uint32_t>(nodeId, outStream); break;
    case USE_VARINT: writeVarInt(nodeId, outStream); break;
    default: writeScalar<uint64_t>(nodeId, outStream); break;
    }
}

// Renumber and simplify
// DFS visit each node, and assign it a new number. Since we only serialize down edges, we only
// need to renumber our down edges (which we can do).
class RenumberAndWriteVisitor : public GRGVisitor {
public:
    explicit RenumberAndWriteVisitor(std::ostream& outStream, size_t idSize, bool allowSimplify)
        : m_outStream(outStream),
          m_idSize(idSize),
          m_allowSimplify(allowSimplify) {}

    // Note: none of the (TS->GRG) simplification comes from a parent node having only one child, it
    // all comes from a child only having one parent. E.g.:
    //      P1  ---P---  P2
    //        \ |  |  | /
    //          A  C  B
    // Only "C" can be simplified in ths above, where C's children get moved up to P. Is this really a
    // simplification that we want?
    bool getChildren(const grgl::GRGPtr& grg, const grgl::NodeID nodeId, NodeIDList& result, NodeData& parentData) {
        bool hasChildren = false;
        const auto& children = grg->getDownEdges(nodeId);
        for (const auto childId : children) {
            if (m_skipNode[childId]) {
                // Coalescences may have occurred at the nodes that we are deleting. These now move to the parent.
                parentData.numIndividualCoals += grg->getNodeData(childId).numIndividualCoals;
                hasChildren |= getChildren(grg, childId, result, parentData);
            } else {
                hasChildren = true;
                result.push_back(childId);
            }
        }
        return hasChildren;
    }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        release_assert(direction == TraversalDirection::DIRECTION_DOWN);
        if (m_nodeVector.empty()) {
            m_nodeVector.resize(grg->numNodes() + 1);
            m_nodeIdMap.resize(grg->numNodes() + 1);
            m_revIdMap.resize(grg->numNodes() + 1);
            m_skipNode.resize(grg->numNodes());
            m_nodeCounter = grg->numSamples();
        }
        if (dfsPass == DfsPass::DFS_PASS_BACK_AGAIN) {
            if (grg->isSample(nodeId)) {
                m_nodeIdMap[nodeId] = nodeId;
                m_revIdMap[nodeId] = nodeId;
                m_nodeVector[nodeId + 1] = 0;
            } else {
                NodeIDList children;
                getChildren(grg, nodeId, children, grg->getNodeData(nodeId));
                const auto numParents = grg->numUpEdges(nodeId);
                const bool extraneousNode =
                    (children.size() <= 1 || numParents <= 1 || (children.size() == 2 && numParents == 2));
                if (!m_allowSimplify || !extraneousNode || grg->nodeHasMutations(nodeId)) {
                    // Save the position of the start of our edges.
                    const NodeID newNodeId = m_nodeCounter++;
                    m_nodeVector[newNodeId + 1] = m_edgeCounter;
                    m_nodeIdMap[nodeId] = newNodeId;
                    m_revIdMap[newNodeId] = nodeId;

                    for (const auto childId : children) {
                        const auto newChildId = m_nodeIdMap.at(childId);
                        writeNodeID(newChildId, m_idSize, m_outStream);
                        m_edgeCounter++;
                    }
                    assert_deserialization(m_outStream.good(), "Writing GRG failed");
                } else {
                    m_skipNode[nodeId] = true;
                    m_nodeIdMap[nodeId] = std::numeric_limits<NodeIDSizeT>::max();
                }
            }
        }
        return true;
    }

    void writeNodeVector() {
        if (m_nodeVector.empty()) {
            m_nodeVector.resize(1);
        }
        release_assert(m_nodeVector.size() >= m_nodeCounter + 1);
        m_nodeVector.resize(m_nodeCounter + 1);
        for (const auto position : m_nodeVector) {
            writeNodeID(position, m_idSize, m_outStream);
        }
    }

    NodeIDSizeT getNumNodes() const { return m_nodeCounter; }

    NodeIDSizeT getNumEdges() const { return m_edgeCounter; }

    NodeID getNewID(const NodeID nodeId) const {
        if (nodeId == INVALID_NODE_ID) {
            return nodeId;
        }
        return m_nodeIdMap[nodeId];
    }

    NodeID getOldID(const NodeID nodeId) const {
        if (nodeId == INVALID_NODE_ID) {
            return nodeId;
        }
        return m_revIdMap[nodeId];
    }

private:
    std::ostream& m_outStream;
    std::vector<NodeIDSizeT> m_nodeVector;
    // Maps old NodeID to new NodeID
    std::vector<NodeIDSizeT> m_nodeIdMap;
    // Maps new NodeID to old NodeID (so the order is new NodeID acscending)
    std::vector<NodeIDSizeT> m_revIdMap;
    std::vector<bool> m_skipNode;
    size_t m_idSize;
    NodeIDSizeT m_nodeCounter{};
    NodeIDSizeT m_edgeCounter{};
    bool m_allowSimplify;
};

template <typename Container>
static void writeNodeList(const Container& list,
                          size_t idSize,
                          std::ostream& outStream,
                          const RenumberAndWriteVisitor* mapper = nullptr) {
    if (USE_VARINT == idSize) {
        writeVarInt(list.size(), outStream);
    } else {
        writeScalar<uint64_t>(list.size(), outStream);
    }
    for (NodeID nodeId : list) {
        if (mapper != nullptr) {
            nodeId = mapper->getNewID(nodeId);
        }
        writeNodeID(nodeId, idSize, outStream);
    }
}

static void writeMutation(const Mutation& mutation, size_t idSize, std::ostream& outStream) {
    writeScalar<uint64_t>(mutation.getPosition(), outStream);
    writeString(USE_VARINT == idSize, mutation.getAllele(), outStream);
    writeString(USE_VARINT == idSize, mutation.getRefAllele(), outStream);
    writeScalar<EXTERNAL_ID>(mutation.getOriginalId(), outStream);
    writeScalar<float>(mutation.getTime(), outStream);
}

// 0...   8...     16...      24...        32...    40...
// <Magic><Id Size><Sample Ct><Mutation Ct><Node Ct><Edge Ct>
void writeGrg(const GRGPtr& grg, std::ostream& outStream, bool useVarInt, bool allowSimplify) {
    assert_deserialization(outStream.good(), "Bad output stream");
    size_t numMutations = grg->getMutations().size();
    release_assert(grg->getNodeMutationPairs().size() == numMutations);
    GRGFileHeader header = {
        GRG_MAGIC,
        GRG_MAJOR_VERSION,
        GRG_MINOR_VERSION,
        static_cast<uint64_t>(useVarInt ? USE_VARINT : getIdBytes(grg)),
        static_cast<uint64_t>(grg->numSamples()),
        static_cast<uint64_t>(numMutations),
        static_cast<uint64_t>(grg->numNodes()),
        static_cast<uint64_t>(grg->numEdges()),
        static_cast<uint16_t>(grg->getPopulations().size()),
        grg->getPloidy(),
        0, /* Unused */
        HAS_INDIVIDUAL_COALS | MUTATION_IDS_ARE_ORDERED,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    };
    // Header
    outStream.write(reinterpret_cast<const char*>(&header), sizeof(header));
    assert_deserialization(outStream.good(), "Writing GRG failed");

    // Write Edges and Nodes
    RenumberAndWriteVisitor visitor(outStream, header.idSize, allowSimplify);
    fastCompleteDFS(grg, visitor);
    visitor.writeNodeVector();
    assert_deserialization(outStream.good(), "Writing GRG failed");
    header.nodeCount = visitor.getNumNodes();
    header.edgeCount = visitor.getNumEdges();
    const auto savedPos = outStream.tellp();
    outStream.seekp(0);
    outStream.write(reinterpret_cast<const char*>(&header), sizeof(header));
    assert_deserialization(outStream.good(), "Writing GRG failed");
    outStream.seekp(savedPos);
    std::cout << "Wrote simplified GRG with:" << std::endl;
    std::cout << "  Nodes: " << header.nodeCount << std::endl;
    std::cout << "  Edges: " << header.edgeCount << std::endl;

    // Mutations
    for (const auto& nodeAndMutId : grg->getMutationsToNodeOrdered()) {
        writeMutation(grg->getMutationById(nodeAndMutId.first), header.idSize, outStream);
        assert_deserialization(outStream.good(), "Writing GRG failed");
        writeNodeID(visitor.getNewID(nodeAndMutId.second), header.idSize, outStream);
        assert_deserialization(outStream.good(), "Writing GRG failed");
    }
    // Populations
    if (!grg->getPopulations().empty()) {
        for (const auto& popDescription : grg->getPopulations()) {
            writeString(useVarInt, popDescription, outStream);
        }
        for (size_t i = 0; i < grg->numSamples(); i++) {
            auto nodeData = grg->getNodeData(i);
            writeScalar<uint16_t>(nodeData.populationId, outStream);
        }
    }
    for (NodeID nodeId = grg->numSamples(); nodeId < grg->numNodes(); nodeId++) {
        writeNodeID(grg->getNodeData(visitor.getOldID(nodeId)).numIndividualCoals, header.idSize, outStream);
    }
    assert_deserialization(outStream.good(), "Writing GRG failed");
}

template <typename T> static inline T readScalar(std::istream& inStream) {
    T simpleValue = 0;
    inStream.read(reinterpret_cast<char*>(&simpleValue), sizeof(simpleValue));
    return simpleValue;
}

static inline std::string readString(bool useVarInt, std::istream& inStream) {
    size_t length = 0;
    if (useVarInt) {
        length = readVarInt(inStream);
    } else {
        length = readScalar<uint64_t>(inStream);
    }
    std::string result;
    result.resize(length);
    // This is slightly sketchy. I should check the C++ standard to see if this
    // behavior is defined.
    inStream.read(const_cast<char*>(result.c_str()), static_cast<std::streamsize>(length));
    return result;
}

static NodeID readNodeID(size_t idSize, std::istream& inStream) {
    switch (idSize) {
    case sizeof(uint8_t): return readScalar<uint8_t>(inStream);
    case sizeof(uint16_t): return readScalar<uint16_t>(inStream);
    case sizeof(uint32_t): return readScalar<uint32_t>(inStream);
    case USE_VARINT: return readVarInt(inStream);
    default: return readScalar<uint64_t>(inStream);
    }
}

// TODO:
// - Switch to using FILE*, and use posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
// - Read the list 16kb at a time, _OR_ just read it directly into the CSR vectors
//   and store things as varint in RAM (as well as disk)
static NodeIDList readNodeIDList(size_t idSize, std::istream& inStream) {
    NodeIDList result;
    size_t length = 0;
    if (USE_VARINT == idSize) {
        length = static_cast<size_t>(readVarInt(inStream));
    } else {
        length = static_cast<size_t>(readScalar<uint64_t>(inStream));
    }
    for (size_t i = 0; i < length; i++) {
        const NodeID nodeId = readNodeID(idSize, inStream);
        assert_deserialization(inStream.good(), "Malformed GRG");
        result.push_back(nodeId);
    }
    return result;
}

static Mutation readMutation(size_t idSize, std::istream& inStream) {
    const auto position = readScalar<uint64_t>(inStream);
    const std::string allele = readString(USE_VARINT == idSize, inStream);
    const std::string refAllele = readString(USE_VARINT == idSize, inStream);
    const auto originalId = readScalar<EXTERNAL_ID>(inStream);
    const auto time = readScalar<float>(inStream);
    return {position, allele, refAllele, time, originalId};
}

void readGrgCommon(const GRGFileHeader& header, const GRGPtr& grg, std::istream& inStream) {
    // Add all the mutations to the nodes.
    for (size_t i = 0; i < header.mutationCount; i++) {
        Mutation mut = readMutation(header.idSize, inStream);
        assert_deserialization(inStream.good(), "Malformed GRG file");
        NodeID nodeId = readNodeID(header.idSize, inStream);
        grg->addMutation(mut, nodeId);
    }
    // Populations
    if (header.populationCount > 0) {
        for (size_t i = 0; i < header.populationCount; i++) {
            grg->addPopulation(std::move(readString(header.idSize == USE_VARINT, inStream)));
            assert_deserialization(inStream.good(), "Malformed GRG file");
        }
        for (size_t i = 0; i < header.sampleCount; i++) {
            const auto popId = readScalar<uint16_t>(inStream);
            grg->getNodeData(i).populationId = popId;
        }
        assert_deserialization(inStream.good(), "Malformed GRG file");
    }

    if ((bool)(header.flags & HAS_INDIVIDUAL_COALS)) {
        for (NodeID nodeId = header.sampleCount; nodeId < header.nodeCount; nodeId++) {
            const NodeIDSizeT coalCount = readNodeID(header.idSize, inStream);
            grg->getNodeData(nodeId).numIndividualCoals = coalCount;
        }
        assert_deserialization(inStream.good(), "Malformed GRG file");
    }
}

MutableGRGPtr readMutableGrg(std::istream& inStream) {
    assert_deserialization(inStream.good(), "Bad input stream");
    // Read header.
    GRGFileHeader header = {};
    inStream.read(reinterpret_cast<char*>(&header), sizeof(header));
    assert_deserialization(GRG_MAGIC == header.magic, "Invalid file header");
    assert_deserialization(GRG_MAJOR_VERSION == header.versionMajor, "Incompatible file major version");
    assert_deserialization(header.nodeCount <= MAX_GRG_NODES, "Malformed GRG file");
    // For backwards compatibility to versions prior to 4.1
    if (header.ploidy == 0 && header.versionMajor == 4 && header.versionMinor == 0) {
        header.ploidy = 2;
    }
    assert_deserialization(header.ploidy != 0, "Malformed GRG file: ploidy was 0");

    // Construct GRG and allocate all the nodes.
    MutableGRGPtr grg = std::make_shared<MutableGRG>(header.sampleCount, header.ploidy, header.nodeCount);
    grg->makeNode(header.nodeCount - header.sampleCount);
    release_assert(grg->numNodes() == header.nodeCount);
    std::vector<NodeID> edges;
    edges.reserve(header.edgeCount);
    for (size_t i = 0; i < header.edgeCount; i++) {
        edges.push_back(readNodeID(header.idSize, inStream));
    }
    assert_deserialization(inStream.good(), "Malformed GRG file");

    NodeID currentNodeID = 0;
    NodeIDSizeT prevEdgePos = 0;
    for (size_t i = 0; i < header.nodeCount + 1; i++) {
        const NodeIDSizeT edgeVectorPos = readNodeID(header.idSize, inStream);
        if (i > 1) {
            currentNodeID = i - 2;
            // If we have any edges...
            if (edgeVectorPos - prevEdgePos > 0) {
                for (size_t j = prevEdgePos; j < edgeVectorPos; j++) {
                    grg->connect(currentNodeID, edges[j]);
                }
                grg->compact(currentNodeID);
            }
        }
        prevEdgePos = edgeVectorPos;
    }
    currentNodeID++;
    assert_deserialization(inStream.good(), "Malformed GRG file");
    assert_deserialization(currentNodeID == (header.nodeCount - 1), "Malformed GRG file");
    for (size_t j = prevEdgePos; j < header.edgeCount; j++) {
        grg->connect(currentNodeID, edges[j]);
    }
    edges.clear();
    grg->compact();

    readGrgCommon(header, grg, inStream);
    if ((bool)(header.flags & MUTATION_IDS_ARE_ORDERED)) {
        grg->m_mutsAreOrdered = true;
    }
    return grg;
}

GRGPtr readImmutableGrg(std::istream& inStream, bool loadUpEdges, bool loadDownEdges) {
    assert_deserialization(inStream.good(), "Bad input stream");
    // Read header.
    GRGFileHeader header = {};
    inStream.read(reinterpret_cast<char*>(&header), sizeof(header));
    assert_deserialization(GRG_MAGIC == header.magic, "Invalid file header");
    assert_deserialization(GRG_MAJOR_VERSION == header.versionMajor, "Incompatible file major version");
    // For backwards compatibility to versions prior to 4.1
    if (header.ploidy == 0 && header.versionMajor == 4 && header.versionMinor == 0) {
        header.ploidy = 2;
    }
    assert_deserialization(header.ploidy != 0, "Malformed GRG file: ploidy was 0");

    // Construct GRG and allocate all the nodes.
    CSRGRGPtr grg =
        std::make_shared<CSRGRG>(header.sampleCount, header.edgeCount, header.nodeCount, header.ploidy, loadUpEdges);

    // The GRG serialization only encodes the down edges, we deserialize them here.
    for (size_t i = 0; i < header.edgeCount; i++) {
        const NodeID targetNode = readNodeID(header.idSize, inStream);
        grg->m_downEdges[i] = targetNode;
        if (loadUpEdges) {
            // We initially set the upPositions to just be the in-degree of each node, which tells us
            // how many "up edges" they will have.
            grg->m_upPositions.at(targetNode + 2)++;
        }
    }
    assert_deserialization(inStream.good(), "Malformed GRG file");
    for (size_t i = 0; i < header.nodeCount + 1; i++) {
        const NodeIDSizeT edgeVectorPos = readNodeID(header.idSize, inStream);
        grg->m_downPositions[i] = edgeVectorPos;
    }

    if (loadUpEdges) {
        // We reconstruct the up edges from the down edges. First we convert the in-degrees
        // to a proper CSR node position vector.
        for (size_t i = 1; i < grg->m_upPositions.size(); i++) {
            grg->m_upPositions[i] += grg->m_upPositions[i - 1];
        }
        release_assert(grg->m_upPositions[grg->m_upPositions.size() - 1] == header.edgeCount);
    }
    grg->finalize();
    if (loadUpEdges) {
        // Next we take a pass over all the down edges and invert them, putting them into the up edge vector.
        std::vector<NodeIDSizeT> perNodeOffsets(grg->m_upPositions.size());
        for (NodeID i = 0; i < grg->numNodes(); i++) {
            for (const NodeID j : grg->getDownEdges(i)) {
                const size_t jStart = grg->m_upPositions[j + 1];
                const size_t jPos = jStart + perNodeOffsets[j];
                grg->m_upEdges.at(jPos) = i;
                perNodeOffsets[j]++;
            }
        }
    }
    if (!loadDownEdges) {
        std::fill(grg->m_downPositions.begin(), grg->m_downPositions.end(), 0);
        grg->m_downEdges.clear();
        grg->m_downEdges.shrink_to_fit();
        grg->finalize();
    }

    readGrgCommon(header, grg, inStream);
    if ((bool)(header.flags & MUTATION_IDS_ARE_ORDERED)) {
        grg->m_mutsAreOrdered = true;
    }
    return grg;
}

} // namespace grgl
