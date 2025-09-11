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
#include "grgl/serialize.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/file_vector.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/node_data.h"
#include "util.h"

#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <numeric>
#include <ostream>
#include <sstream>

#include "vbyte.h"

#define assert_deserialization(condition, msg)                                                                         \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::stringstream errMsg;                                                                                  \
            errMsg << "Failed (de)serialization: " << (msg);                                                           \
            throw SerializationFailure(msg);                                                                           \
        }                                                                                                              \
    } while (0)

// HOW TO UPDATE GRG FILE FORMAT
//
// There are two types of data that can be added to the GRG file: REQUIRED or OPTIONAL. The
// REQUIRED data is always there if the GRG file has a version equal-to-or-later-than when
// the type of data was added to the format specification. The OPTIONAL data is only there
// sometimes, and uses another piece of information in the file to determine whether it is
// there or not.
//
// When adding data of either type to the file format, the GRG_FILE_MINOR_VERSION should be
// incremented. This indicates a backwards/forwards-compatible change. The GRG_FILE_MAJOR_VERSION
// should only be incremented when the change is breaking, for example when the layout of an
// existing data structure (e.g., the edge lists) changes. This should be rare.
//
// REQUIRED data can be handled during deserialization by checking versionIsAtLeast(), and if
// the version of the file is at least as high as when the format was modified then the data
// must be there.
//
// OPTIONAL data is usually handled by adding a new GRG_FLAG_* bitmask value, and making sure
// that GRGFileHeader::flags contains that bit set when the data is present, and does not
// otherwise. The flags are always guaranteed to be set to 0 for older versions of the file
// format that have not yet added the new flag. Flags can only be "reclaimed" (i.e., change
// their meaning) when the MAJOR file format version changes.

namespace grgl {

static inline bool
versionIsAtLeast(const GRGFileHeader& header, const uint16_t versionMajor, const uint16_t versionMinor) {
    return (header.versionMajor > versionMajor) ||
           (header.versionMajor == versionMajor && header.versionMinor >= versionMinor);
}

static inline bool
versionIsUpTo(const GRGFileHeader& header, const uint16_t versionMajor, const uint16_t versionMinor) {
    return (header.versionMajor < versionMajor) ||
           (header.versionMajor == versionMajor && header.versionMinor <= versionMinor);
}

#define HAS_SMALL_CSR_IDX versionIsUpTo(header, 5, 1)

static inline bool hasFlag(const GRGFileHeader& header, const uint64_t flagValue) {
    return static_cast<bool>(header.flags & flagValue);
}

template <typename T> static inline void writeScalar(const T& intValue, std::ostream& outStream) {
    outStream.write(reinterpret_cast<const char*>(&intValue), sizeof(intValue));
}

static inline size_t writeVarInt32(uint32_t value, std::ostream& outStream) {
    static constexpr size_t maxSize32 = 5;
    static uint8_t buffer[maxSize32];
    const std::streamsize written = static_cast<std::streamsize>(vbyte_write_uint32(&buffer[0], value));
    release_assert(written <= maxSize32);
    outStream.write(reinterpret_cast<const char*>(&buffer[0]), written);
    return written;
}

static inline size_t writeVarInt64(uint64_t value, std::ostream& outStream) {
    static constexpr size_t maxSize64 = 10;
    static uint8_t buffer[maxSize64];
    const std::streamsize written = static_cast<std::streamsize>(vbyte_write_uint64(&buffer[0], value));
    release_assert(written <= maxSize64);
    outStream.write(reinterpret_cast<const char*>(&buffer[0]), written);
    return written;
}

static inline void writeString(const std::string& value, std::ostream& outStream) {
    static constexpr size_t maxSize32 = 5;
    static uint8_t buffer[maxSize32];
    const std::streamsize written = static_cast<std::streamsize>(vbyte_write_uint32(&buffer[0], value.size()));
    release_assert(written <= maxSize32);
    outStream.write(reinterpret_cast<const char*>(&buffer[0]), written);
    outStream.write(value.c_str(), static_cast<std::streamsize>(value.size()));
}

static void writeNodeID(const NodeID nodeId, std::ostream& outStream) { writeVarInt32(nodeId, outStream); }

// Renumber and simplify
// DFS visit each node, and assign it a new number. Since we only serialize down edges, we only
// need to renumber our down edges (which we can do).
//
// It is important that this traversal NOT modify the input GRG. Ideally, the passed-in GRG would be
// const (enforced by the compiler), but the use of FileVector prevents this. Should consider templating
// the serialization so we can get constness, because:
// 1. In-memory GRGs need to be shared by all threads for performance reasons. Therefore must be unmodified.
// 2. On-disk GRGs need to be separate copies of the same "GRG", as they will modify a bunch of book-keeping
//    about where they are reading from.
class RenumberAndWriteVisitor : public GRGVisitor {
public:
    explicit RenumberAndWriteVisitor(std::ostream& outStream, NodeIDSizeT numNodes, bool allowSimplify)
        : m_outStream(outStream),
          m_edgeCSR(numNodes),
          m_allowSimplify(allowSimplify) {}

    // Returns the _new_ ID for children, not the old NodeID.
    //
    // Note: none of the (TS->GRG) simplification comes from a parent node having only one child, it
    // all comes from a child only having one parent. E.g.:
    //      P1  ---P---  P2
    //        \ |  |  | /
    //          A  C  B
    // Only "C" can be simplified in ths above, where C's children get moved up to P. Is this really a
    // simplification that we want?
    bool getChildren(const grgl::GRGPtr& grg, const grgl::NodeID nodeId, NodeIDList& result, NodeIDSizeT& parentCoals) {
        bool hasChildren = false;
        const auto& children = grg->getDownEdges(nodeId);
        for (const auto childId : children) {
            if (!m_keepBeneath[childId]) {
                continue;
            }
            const auto& newChildId = m_nodeIdMap[childId];
            if (INVALID_NODE_ID == newChildId) {
                // We need to propagate coalescent counts up to parents, when we delete nodes. However, we only want
                // to do this for immediate children that have been deleted (because we are doing this recursively...)
                if (parentCoals != COAL_COUNT_NOT_SET) {
                    const NodeIDSizeT childCoals = grg->getNumIndividualCoals(childId);
                    if (childCoals != COAL_COUNT_NOT_SET) {
                        parentCoals += childCoals;
                    }
                }
                hasChildren |= getChildren(grg, childId, result, parentCoals);
            } else {
                hasChildren = true;
                result.push_back(newChildId);
            }
        }
        return hasChildren;
    }

    // Does the given node have mutations (after filtering) or other important properties?
    bool shouldKeep(const grgl::GRGPtr& grg, const grgl::NodeID nodeId) const {
        // We keep all missingness nodes _that are visited during traversal_.
        if (grg->hasMissingData() && m_isMissingnessNode.at(nodeId)) {
            return true;
        }
        if (grg->nodeHasMutations(nodeId)) {
            if (m_keepMutations.empty()) {
                return true;
            }
            for (auto mutId : grg->getMutationsForNode(nodeId)) {
                if (m_keepMutations[mutId]) {
                    return true;
                }
            }
        }
        return false;
    }

    // Set which samples we want to keep; if never called then we keep all samples.
    NodeIDSizeT setKeepSamples(const grgl::GRGPtr& grg, grgl::NodeIDList sampleIDList) {
        NodeIDSizeT numSamples = 0;
        NodeID prevSampleId = INVALID_NODE_ID;
        std::sort(sampleIDList.begin(), sampleIDList.end());
        m_nodeIdMap.resize(sampleIDList.back() + 1, INVALID_NODE_ID);
        release_assert(m_revIdMap.empty());
        for (const auto sampleId : sampleIDList) {
            if (!grg->isSample(sampleId)) {
                throw ApiMisuseFailure("Not a valid sampleId");
            }
            if (sampleId != prevSampleId) {
                const NodeID newSampleId = numSamples++;
                m_nodeIdMap.at(sampleId) = newSampleId;
                m_revIdMap.push_back(sampleId);
                release_assert(m_revIdMap.size() == newSampleId + 1);
            }
            prevSampleId = sampleId;
        }
        m_nodeCounter = numSamples;
        m_newSampleCount = numSamples;
        return numSamples;
    }

    // Set which mutations we want to keep; if never called then we keep all mutations.
    // Returns the list of seeds that need to be traversed (downwards) in order to ensure that
    // all of the children of the kept mutations are kept.
    NodeIDList setKeepMutations(const grgl::GRGPtr& grg, const grgl::NodeIDList& mutationIDList) {
        NodeIDList result;
        NodeIDSet keptMutNodes;
        if (!mutationIDList.empty()) {
            m_keepMutations.resize(grg->numMutations(), false);
            m_isMissingnessNode.resize(grg->numNodes(), false);
            for (const auto id : mutationIDList) {
                m_keepMutations.at(id) = true;
            }
            for (const auto& item : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
                if (m_keepMutations[std::get<1>(item)]) {
                    const NodeID mutNode = std::get<0>(item);
                    if (mutNode != INVALID_NODE_ID) {
                        keptMutNodes.emplace(mutNode);
                    }
                    const NodeID missingnessNode = std::get<2>(item);
                    if (missingnessNode != INVALID_NODE_ID) {
                        keptMutNodes.emplace(missingnessNode);
                        m_isMissingnessNode.at(missingnessNode) = true;
                    }
                }
            }
        } else {
            m_keepMutations.resize(0);
            m_keepMutations.shrink_to_fit();
            m_isMissingnessNode.resize(grg->numNodes(), false);
            for (const auto& item : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
                const NodeID missingnessNode = std::get<2>(item);
                m_isMissingnessNode.at(missingnessNode) = true;
            }
        }
        for (const auto& nodeId : keptMutNodes) {
            result.emplace_back(nodeId);
        }
        return std::move(result);
    }

    bool keepMutation(const MutationId mutId) { return m_keepMutations.empty() || m_keepMutations[mutId]; }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        constexpr size_t FLUSH_THRESHOLD = 1024 * 1024; // 1MB

        if (m_keepBeneath.empty()) {
            m_keepBeneath.resize(grg->numNodes(), false);
            // If nodeIdMap is not empty, we already mapped the samples. Otherwise we need to do it.
            const bool noSamplesMapped = m_nodeIdMap.empty();
            m_nodeIdMap.resize(grg->numNodes(), INVALID_NODE_ID);
            if (noSamplesMapped) {
                const NodeIDSizeT numSamples = grg->numSamples();
                m_newSampleCount = numSamples;
                m_nodeCounter = numSamples;
                std::iota(m_nodeIdMap.begin(), m_nodeIdMap.begin() + numSamples, 0);
                m_revIdMap.resize(grg->numSamples());
                for (NodeID sampleId = 0; sampleId < numSamples; sampleId++) {
                    m_newNodeData.setPopId(sampleId, grg->getPopulationId(nodeId));
                }
            }
            if (m_isMissingnessNode.empty() && grg->hasMissingData()) {
                m_isMissingnessNode.resize(grg->numNodes(), false);
                for (const auto& item : grg->getNodesAndMutations<GRG::NodeMutMiss>()) {
                    const NodeID missingnessNode = std::get<2>(item);
                    if (missingnessNode != INVALID_NODE_ID) {
                        m_isMissingnessNode.at(missingnessNode) = true;
                    }
                }
            }
        }
        if (dfsPass != DfsPass::DFS_PASS_THERE) {
            // A flag we use to track which nodes we have visited, any unvisited nodes are filtered out during
            // GRG serialization.
            m_keepBeneath[nodeId] = true;

            if (grg->isSample(nodeId)) {
                const NodeID newNodeId = m_nodeIdMap[nodeId];
                m_revIdMap[newNodeId] = nodeId;
                m_newNodeData.setPopId(newNodeId, grg->getPopulationId(nodeId));
            } else {
                // We should never visit a node more than once
                assert(m_nodeIdMap[nodeId] == INVALID_NODE_ID);

                NodeIDList children;
                NodeIDSizeT parentCoals = grg->getNumIndividualCoals(nodeId);
                getChildren(grg, nodeId, children, parentCoals);
                const auto numParents = grg->numUpEdges(nodeId);
                // Nodes meeting these criteria only make the graph larger, so simplify them out.
                const bool extraneousNode =
                    (children.size() <= 1 || numParents <= 1 || (children.size() == 2 && numParents == 2));
                // Nodes with no children but with mutations arise when removing samples, and we keep the mutation
                // but ditch the node.
                const bool keepMutNode = !children.empty() && shouldKeep(grg, nodeId);
                if (!m_allowSimplify || !extraneousNode || keepMutNode) {
                    // Save the position of the start of our edges.
                    const NodeID newNodeId = m_nodeCounter++;
                    m_nodeIdMap[nodeId] = newNodeId;
                    m_revIdMap.push_back(nodeId);
                    assert(m_revIdMap.size() == newNodeId + 1);

                    std::sort(children.begin(), children.end());
                    m_edgeCounter += children.size();
                    m_bytesWritten += m_edgeCSR.setData(newNodeId, children);

                    if (m_edgeCSR.inmemoryBucketSize() > FLUSH_THRESHOLD) {
                        m_edgeCSR.flushBuckets(m_outStream);
                        assert_deserialization(m_outStream.good(), "Writing GRG failed");
                    }
                    release_assert(newNodeId >= m_newSampleCount);
                    m_newNodeData.allocNumCoals((newNodeId - m_newSampleCount) + 1);
                    m_newNodeData.setNumCoals(m_newSampleCount, newNodeId, parentCoals);
                } else {
                    m_nodeIdMap[nodeId] = INVALID_NODE_ID;
                }
            }
        }
        return true;
    }

    void finalize() {
        release_assert(m_revIdMap.size() == m_nodeCounter);
        m_edgeCSR.finalizeNodes(m_nodeCounter);
        m_edgeCSR.flushBuckets(m_outStream);
        assert_deserialization(m_outStream.good(), "Writing GRG failed");
    }

    NodeIDSizeT getNumNodes() const { return m_nodeCounter; }

    EdgeSizeT getNumEdges() const { return m_edgeCounter; }

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

private:
    std::ostream& m_outStream;
    // Maps old NodeID to new NodeID
    std::vector<NodeIDSizeT> m_nodeIdMap;
    NodeIDSizeT m_newSampleCount{};
    NodeIDSizeT m_nodeCounter{};
    EdgeSizeT m_edgeCounter{};
    bool m_allowSimplify;
};

std::pair<NodeIDSizeT, EdgeSizeT>
simplifyAndSerialize(const GRGPtr& grg, std::ostream& outStream, const GRGOutputFilter& filter, bool allowSimplify) {
    assert_deserialization(outStream.good(), "Bad output stream");
    const size_t numMutations = grg->getMutations().size();
    release_assert(grg->getNodesAndMutations().size() == numMutations);
    const auto bpRange = grg->getBPRange();
    uint64_t flags = (!grg->isPhased() ? GRG_FLAG_UNPHASED : 0) | (grg->hasMissingData() ? GRG_FLAG_HAS_MISSING : 0);
    GRGFileHeader header = {
        GRG_FILE_MAGIC,
        GRG_FILE_MAJOR_VERSION,
        GRG_FILE_MINOR_VERSION,
        grg->getPloidy(),
        static_cast<uint16_t>(grg->getPopulations().size()),
        static_cast<uint64_t>(grg->numSamples()),
        static_cast<uint64_t>(numMutations),
        static_cast<uint64_t>(grg->numNodes()),
        static_cast<uint64_t>(grg->numEdges()),
        filter.isSpecified() ? filter.bpRange.first : grg->getSpecifiedBPRange().first,
        filter.isSpecified() ? filter.bpRange.second : grg->getSpecifiedBPRange().second,
        flags,
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
    // Save a spot for number of bytes written for edges (buckets)
    const size_t edgeBytePosition = outStream.tellp();
    writeScalar<uint64_t>(0, outStream);
    RenumberAndWriteVisitor visitor(outStream, grg->numNodes(), allowSimplify);
    if (filter.isSpecified()) {
        if (filter.direction == TraversalDirection::DIRECTION_UP) {
            header.sampleCount = visitor.setKeepSamples(grg, filter.seedList);
            grg->visitTopo(visitor, TraversalDirection::DIRECTION_UP, filter.seedList);
        } else {
            assert(filter.direction == TraversalDirection::DIRECTION_DOWN);
            const NodeIDList keptSeedNodes = visitor.setKeepMutations(grg, filter.seedList);
            grg->visitDfs(visitor, TraversalDirection::DIRECTION_DOWN, keptSeedNodes);
        }
    } else {
        fastCompleteDFS(grg, visitor);
    }
    visitor.finalize();
    auto savedPos = outStream.tellp();
    outStream.seekp(edgeBytePosition);
    writeScalar<uint64_t>(visitor.m_bytesWritten, outStream);
    outStream.seekp(savedPos);
    visitor.m_edgeCSR.flushIndexes(outStream);
    assert_deserialization(outStream.good(), "Writing GRG failed");

    header.nodeCount = visitor.getNumNodes();
    release_assert(header.edgeCount >= visitor.getNumEdges()); // Simplification should never increase edges
    header.edgeCount = visitor.getNumEdges();
    savedPos = outStream.tellp();
    outStream.seekp(0);
    outStream.write(reinterpret_cast<const char*>(&header), sizeof(header));
    assert_deserialization(outStream.good(), "Writing GRG failed");
    outStream.seekp(savedPos);

    ///////////// Mutations /////////////
    const auto mutToNode = grg->getMutationsToNodeOrdered<GRG::MutNodeMiss>();
    // Nodes in (new) MutationID order.
    {
        EagerFileVector<NodeID> mutNodes;
        mutNodes.reserve(mutToNode.size());
        header.mutationCount = 0;
        for (const auto& item : mutToNode) {
            const MutationId& mutId = std::get<0>(item);
            if (visitor.keepMutation(mutId)) {
                const NodeID nodeId = visitor.getNewID(std::get<1>(item));
                mutNodes.push_back(nodeId);
                header.mutationCount++;
            }
        }
        mutNodes.flush(outStream);
        assert_deserialization(outStream.good(), "Writing GRG failed");
    }
    // Mutations in (new) MutationID order.
    {
        std::vector<MutationId> deferredAlleles;
        for (const auto& item : mutToNode) {
            const MutationId& mutId = std::get<0>(item);
            if (visitor.keepMutation(mutId)) {
                const auto& mutation = grg->getMutationById(mutId);
                if (!mutation.isShort()) {
                    deferredAlleles.push_back(mutId);
                }
                writeScalar<Mutation>(mutation, outStream);
            }
        }
        assert_deserialization(outStream.good(), "Writing GRG failed");
        for (const auto mutId : deferredAlleles) {
            writeString(grg->getMutationById(mutId).alleleStorageString(), outStream);
        }
        assert_deserialization(outStream.good(), "Writing GRG failed");
    }

    // Populations
    if (!grg->getPopulations().empty()) {
        for (const auto& popDescription : grg->getPopulations()) {
            writeString(popDescription, outStream);
        }
        visitor.m_newNodeData.writePopIds(outStream);
    }
    // Individual coalescences.
    visitor.m_newNodeData.writeCoalCounts(outStream, header.sampleCount, header.nodeCount);
    assert_deserialization(outStream.good(), "Writing GRG failed");

    // Individual IDs are optional. If we have them, they are all stored together in an unsorted,
    // unencoded, CSR table.
    if (grg->hasIndividualIds()) {
        grg->m_individualIds.finalizeNodes();
        writeScalar<uint64_t>(grg->m_individualIds.numValueBytes(), outStream);
        grg->m_individualIds.flushBuckets(outStream, /*clearData=*/false);
        grg->m_individualIds.flushIndexes(outStream, /*clearData=*/false);
        header.flags |= GRG_FLAG_HAS_INDIV_IDS;
    }

    if (grg->hasMissingData()) {
        EagerFileVector<NodeID> missingnessNodes;
        missingnessNodes.reserve(mutToNode.size());
        for (const auto& item : mutToNode) {
            const MutationId& mutId = std::get<0>(item);
            if (visitor.keepMutation(mutId)) {
                const NodeID missNodeId = visitor.getNewID(std::get<2>(item));
                missingnessNodes.push_back(missNodeId);
            }
        }
        missingnessNodes.flush(outStream);
        assert_deserialization(outStream.good(), "Writing GRG failed");
    }

    // Rewrite the header, since some fields can change.
    outStream.seekp(0);
    outStream.write(reinterpret_cast<const char*>(&header), sizeof(header));
    assert_deserialization(outStream.good(), "Writing GRG failed");
    return {(NodeIDSizeT)header.nodeCount, (EdgeSizeT)header.edgeCount};
}

std::pair<NodeIDSizeT, EdgeSizeT> writeGrg(const GRGPtr& grg, std::ostream& out, bool allowSimplify) {
    GRGOutputFilter emptyFilter;
    return simplifyAndSerialize(grg, out, emptyFilter, allowSimplify);
}

template <typename T> static inline T readScalar(std::istream& inStream) {
    T simpleValue = 0;
    inStream.read(reinterpret_cast<char*>(&simpleValue), sizeof(simpleValue));
    return std::move(simpleValue);
}

// FIXME this is janky, mostly because vbyte_read_uint32 doesn't provide a file-based read.
static inline std::string readString(std::istream& inStream) {
    std::string result;
    static uint8_t buffer[5];
    inStream.read(reinterpret_cast<char*>(&buffer[0]), sizeof(buffer));
    uint32_t length = 0;
    size_t readBytes = vbyte_read_uint32(&buffer[0], &length);
    inStream.seekg(-std::streamoff(5 - readBytes), std::ios_base::cur);
    result.resize(length);
    // This is slightly sketchy. I should check the C++ standard to see if this
    // behavior is defined.
    inStream.read(const_cast<char*>(result.c_str()), static_cast<std::streamsize>(length));
    return result;
}

void readGrgCommon(const GRGFileHeader& header, const GRGPtr& grg, IFSPointer& inStream) {
    EagerFileVector<NodeID> mutNodes(inStream, inStream->tellg(), header.mutationCount);
    EagerFileVector<Mutation> mutations(inStream, inStream->tellg(), header.mutationCount);
    assert_deserialization(mutNodes.size() == mutations.size(), "Malformed GRG file (inconsistent mutations)");
    std::vector<MutationId> deferredAlleles;
    for (MutationId mutId = 0; mutId < header.mutationCount; mutId++) {
        if (!mutations[mutId].isShort()) {
            deferredAlleles.push_back(mutId);
        }
    }
    assert_deserialization(inStream->good(), "Malformed GRG file");
    for (const auto mutId : deferredAlleles) {
        mutations.ref(mutId).setAlleleStorageString(readString(*inStream));
    }
    grg->m_mutations = std::move(mutations);

    // Populations
    if (header.populationCount > 0) {
        for (size_t i = 0; i < header.populationCount; i++) {
            grg->addPopulation(std::move(readString(*inStream)));
            assert_deserialization(inStream->good(), "Malformed GRG file");
        }
        grg->m_nodeData.readPopIds(inStream, header.sampleCount);
        assert_deserialization(inStream->good(), "Malformed GRG file");
    }

    grg->m_nodeData.readCoalCounts(inStream, header.sampleCount, header.nodeCount);

    // Individual IDs are optional. If we have them, they are all stored together in an unsorted,
    // unencoded, CSR table.
    if (hasFlag(header, GRG_FLAG_HAS_INDIV_IDS)) {
        const size_t indivIdBytes = readScalar<uint64_t>(*inStream);
        EagerFileVector<uint8_t> stringValues(inStream, inStream->tellg(), indivIdBytes);
        EagerFileVector<EdgeSizeT> stringPointers;
        if (HAS_SMALL_CSR_IDX) {
            EagerFileVector<uint32_t> smallPointers(inStream, inStream->tellg(), grg->numIndividuals() + 1);
            stringPointers = convertTo64BitFV(smallPointers);
        } else {
            EagerFileVector<EdgeSizeT> bigPointers(inStream, inStream->tellg(), grg->numIndividuals() + 1);
            stringPointers = std::move(bigPointers);
        }
        CSRStringTable individualIds(std::move(stringPointers), std::move(stringValues), indivIdBytes);
        assert_deserialization(individualIds.numNodes() == grg->numIndividuals(), "Malformed GRG file: individual IDs");
        grg->m_individualIds = std::move(individualIds);
    }

    if (hasFlag(header, GRG_FLAG_HAS_MISSING)) {
        EagerFileVector<NodeID> missNodes(inStream, inStream->tellg(), header.mutationCount);
        assert_deserialization(missNodes.size() == mutNodes.size(), "Invalid missingness node size");
        for (MutationId mutId = 0; mutId < mutNodes.size(); mutId++) {
            grg->m_mutIdsByNodeIdAndMiss.emplace_back(mutNodes[mutId], mutId, missNodes[mutId]);
        }
    } else {
        for (MutationId mutId = 0; mutId < mutNodes.size(); mutId++) {
            grg->m_mutIdsByNodeId.emplace_back(mutNodes[mutId], mutId);
        }
    }
    grg->sortMutIdsByNodeID();

    assert_deserialization(inStream->good(), "Malformed GRG file");
    grg->m_mutsAreOrdered = true;
}

MutableGRGPtr readMutableGrg(IFSPointer& inStream) {
    assert_deserialization(inStream->good(), "Bad input stream");
    // Read header.
    GRGFileHeader header = {};
    inStream->read(reinterpret_cast<char*>(&header), sizeof(header));
    assert_deserialization(GRG_FILE_MAGIC == header.magic, "Invalid file header");
    assert_deserialization(GRG_FILE_MAJOR_VERSION == header.versionMajor, "Incompatible file major version");
    release_assert(versionIsAtLeast(header, GRG_FILE_MAJOR_VERSION, 0)); // Just testing functionality.
    assert_deserialization(header.nodeCount <= MAX_GRG_NODES, "Malformed GRG file");
    assert_deserialization(header.ploidy != 0, "Malformed GRG file: ploidy was 0");

    // Construct GRG and allocate all the nodes.
    MutableGRGPtr grg = std::make_shared<MutableGRG>(
        header.sampleCount, header.ploidy, !hasFlag(header, GRG_FLAG_UNPHASED), header.nodeCount);
    grg->setSpecifiedBPRange({header.rangeStart, header.rangeEnd});
    grg->makeNode(header.nodeCount - header.sampleCount);
    release_assert(grg->numNodes() == header.nodeCount);
    {
        const size_t edgeBytes = readScalar<uint64_t>(*inStream);
        const std::streamoff beforeEdges = inStream->tellg();
        LazyFileVector<uint8_t> edges(inStream, inStream->tellg(), edgeBytes);
        inStream->seekg(beforeEdges + (std::streamoff)edgeBytes);
        assert_deserialization(inStream->good(), "Malformed GRG file");
        EagerFileVector<EdgeSizeT> nodes;
        if (HAS_SMALL_CSR_IDX) {
            EagerFileVector<uint32_t> smallNodes(inStream, inStream->tellg(), header.nodeCount + 1);
            nodes = convertTo64BitFV(smallNodes);
        } else {
            EagerFileVector<EdgeSizeT> bigNodes(inStream, inStream->tellg(), header.nodeCount + 1);
            nodes = std::move(bigNodes);
        }
        LazyCSREdges32 edgeCSR(std::move(nodes), std::move(edges), header.edgeCount);
        assert_deserialization(inStream->good(), "Malformed GRG file");
        const std::streamoff afterNodes = inStream->tellg();

        for (NodeID nodeId = 0; nodeId < header.nodeCount; nodeId++) {
            NodeIDList edges;
            edgeCSR.getData(nodeId, edges);
            for (const NodeID childId : edges) {
                grg->connect(nodeId, childId);
            }
            grg->compact(nodeId);
        }
        assert_deserialization(inStream->good(), "Malformed GRG file");
        grg->compact();
        // The LazyFileVector above is seeking all over the place, so seek to the end of the
        // edges + nodes before reading the rest of the file.
        inStream->seekg(afterNodes);
    }

    readGrgCommon(header, grg, inStream);
    grg->m_nodesAreOrdered = true;
    return grg;
}

GRGPtr readImmutableGrg(IFSPointer& inStream, bool loadUpEdges) {
    assert_deserialization(inStream->good(), "Bad input stream");
    // Read header.
    GRGFileHeader header = {};
    inStream->read(reinterpret_cast<char*>(&header), sizeof(header));
    assert_deserialization(GRG_FILE_MAGIC == header.magic, "Invalid file header");
    assert_deserialization(GRG_FILE_MAJOR_VERSION == header.versionMajor, "Incompatible file major version");
    release_assert(versionIsAtLeast(header, GRG_FILE_MAJOR_VERSION, 0)); // Just testing functionality.
    assert_deserialization(header.ploidy != 0, "Malformed GRG file: ploidy was 0");

    // Construct GRG and allocate all the nodes.
    CSRGRGPtr grg = std::make_shared<CSRGRG>(
        header.sampleCount, header.nodeCount, header.ploidy, loadUpEdges, !hasFlag(header, GRG_FLAG_UNPHASED));
    grg->setSpecifiedBPRange({header.rangeStart, header.rangeEnd});

    const size_t edgeBytes = readScalar<uint64_t>(*inStream);
    EagerFileVector<uint8_t> edges(inStream, inStream->tellg(), edgeBytes);
    EagerFileVector<EdgeSizeT> nodes;
    if (HAS_SMALL_CSR_IDX) {
        EagerFileVector<uint32_t> smallNodes(inStream, inStream->tellg(), header.nodeCount + 1);
        nodes = convertTo64BitFV(smallNodes);
    } else {
        EagerFileVector<EdgeSizeT> bigNodes(inStream, inStream->tellg(), header.nodeCount + 1);
        nodes = std::move(bigNodes);
    }
    EagerCSREdges32 edgeCSR(std::move(nodes), std::move(edges), header.edgeCount);
    assert_deserialization(inStream->good(), "Malformed GRG file");

    if (loadUpEdges) {
        // Inverting the edges without using a ton of RAM is tricky, hence all the steps below.
        // Even now it's not very compact; to make it more compact we need to either (a) write
        // upEdges[] to disk before coverted or (b) use the same vector for upEdges and
        // upEdgeCSR, which means swapping the order.
        std::vector<NodeIDSizeT> upDegree(header.nodeCount, 0);
        for (NodeID srcId = 0; srcId < header.nodeCount; srcId++) {
            NodeIDList targets;
            edgeCSR.getData(srcId, targets);
            for (const NodeID tgtId : targets) {
                upDegree[tgtId] += 1;
            }
        }
        std::vector<EdgeSizeT> upIndexes(header.nodeCount, 0);
        EdgeSizeT prevIndex = 0;
        for (NodeID i = header.nodeCount; i > 0; i--) {
            const NodeID nodeId = i - 1;
            upIndexes[nodeId] = prevIndex + upDegree[nodeId];
            prevIndex = upIndexes[nodeId];
        }
        release_assert(upIndexes[0] == header.edgeCount);

        NodeIDList upEdges(header.edgeCount);
        for (NodeID srcId = 0; srcId < header.nodeCount; srcId++) {
            NodeIDList targets;
            edgeCSR.getData(srcId, targets);
            for (NodeID tgtId : targets) {
                EdgeSizeT& index = upIndexes[tgtId];
                upEdges[--index] = srcId;
            }
        }
        EagerCSREdges32 upEdgeCSR(header.nodeCount);
        for (NodeID nodeId = 0; nodeId < header.nodeCount; nodeId++) {
            const size_t index = upIndexes[nodeId];
            assert(upEdges.size() - index == upDegree[nodeId]);
            std::sort(&upEdges[index], &upEdges[upEdges.size()]);
            upEdgeCSR.setData(nodeId, upEdges.data(), {index, upEdges.size()});
            upEdges.resize(index);
        }
        upEdges.shrink_to_fit();
        upEdgeCSR.finalizeNodes();
        grg->m_upEdges = std::move(upEdgeCSR);
    }
    grg->m_downEdges = std::move(edgeCSR);

    readGrgCommon(header, grg, inStream);
    return grg;
}

} // namespace grgl
