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
#ifndef GRG_NODE_DATA_H
#define GRG_NODE_DATA_H

#include "file_vector.h"
#include "grgl/common.h"
#include "grgl/grgnode.h"
#include <cassert>
#include <limits>

namespace grgl {

using PopulationID = uint16_t;
constexpr PopulationID POPULATION_UNSPECIFIED = std::numeric_limits<PopulationID>::max();

constexpr NodeIDSizeT COAL_COUNT_NOT_SET = std::numeric_limits<NodeIDSizeT>::max();

class NodeDataContainer {
public:
    PopulationID getPopId(NodeID nodeId) {
        if (nodeId < m_nodeData_popId.size()) {
            return m_nodeData_popId[nodeId];
        }
        return POPULATION_UNSPECIFIED;
    }

    void setPopId(NodeID nodeId, PopulationID popId) {
        if (m_nodeData_popId.size() <= nodeId) {
            m_nodeData_popId.resize(nodeId + 1, POPULATION_UNSPECIFIED);
        }
        m_nodeData_popId.ref(nodeId) = popId;
    }

    NodeIDSizeT getNumCoals(const NodeIDSizeT numSamples, NodeID nodeId) {
        if (nodeId < numSamples) {
            return 0;
        }
        const NodeID index = nodeId - numSamples;
        if (index < m_nodeData_numCoals.size()) {
            const auto coalCount = m_nodeData_numCoals.read_atomic(index);
            assert(coalCount == COAL_COUNT_NOT_SET || coalCount < numSamples);
            return coalCount;
        }
        return COAL_COUNT_NOT_SET;
    }

    void allocNumCoals(const NodeIDSizeT nonSampleNodes) {
        m_nodeData_numCoals.resize(nonSampleNodes, COAL_COUNT_NOT_SET);
    }

    void setNumCoals(const NodeIDSizeT numSamples, NodeID nodeId, NodeIDSizeT coalCount) {
        if (nodeId >= numSamples) {
            const NodeID index = nodeId - numSamples;
            release_assert(index < m_nodeData_numCoals.size());
            assert(coalCount == COAL_COUNT_NOT_SET || coalCount < numSamples);
            m_nodeData_numCoals.store_atomic(index, coalCount);
        }
    }

    size_t writePopIds(std::ostream& outFile) { return m_nodeData_popId.flush(outFile, /*keepInRam=*/true); }

    void readPopIds(IFSPointer& inFile, const NodeIDSizeT numSamples) {
        m_nodeData_popId = std::move(EagerFileVector<PopulationID>(inFile, inFile->tellg(), numSamples));
    }

    size_t writeCoalCounts(std::ostream& outFile, const NodeIDSizeT numSamples, const NodeIDSizeT numNodes) {
        release_assert(numNodes - numSamples == m_nodeData_numCoals.size());
        return m_nodeData_numCoals.flush(outFile, /*keepInRam=*/true);
    }

    void readCoalCounts(IFSPointer& inFile, const NodeIDSizeT numSamples, const NodeIDSizeT numNodes) {
        const NodeIDSizeT numCoals = numNodes - numSamples;
        m_nodeData_numCoals = std::move(EagerFileVector<NodeIDSizeT>(inFile, inFile->tellg(), numCoals));
    }

private:
    // Node data. This will be replaced by a general purpose store where users can define their
    // own datatypes.
    // Number of individuals that coalesce at the (non-sample) node.
    EagerFileVector<NodeIDSizeT> m_nodeData_numCoals;
    // The ID of the population (sample node only)
    EagerFileVector<PopulationID> m_nodeData_popId;
};

} // namespace grgl

#endif /* GRG_NODE_DATA_H */
