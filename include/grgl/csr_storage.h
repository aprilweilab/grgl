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
#ifndef GRG_CSR_STORAGE_H
#define GRG_CSR_STORAGE_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "grgl/common.h"
#include "grgl/file_vector.h"
#include "grgl/grgnode.h"

#include "vbyte.h"

namespace grgl {

/**
 * Dispatcher class that can be template-specialized to the various combinations of T,Encoded,Sorted.
 * Used by CSRStorageImm.
 */
template <typename T, bool Encoded, bool Sorted> class VByteDispatcher {
public:
    inline size_t maxElemSize() const;
    inline size_t compressSingle(T in, uint8_t* out) const;
    inline size_t decompressSingle(const uint8_t* in, T* out) const;
    inline size_t compressData(const T* in, uint8_t* out, size_t length) const;
    inline size_t decompressData(const uint8_t* in, T* out, size_t length) const;
    inline T selectData(const uint8_t* in, size_t size, size_t index) const;
    inline size_t searchData(const uint8_t* in, size_t length, T value) const;
};

template <> inline size_t VByteDispatcher<uint8_t, false, false>::maxElemSize() const { return sizeof(uint8_t); }

template <> inline size_t VByteDispatcher<uint32_t, false, false>::maxElemSize() const { return sizeof(uint32_t); }

template <> inline size_t VByteDispatcher<uint32_t, true, true>::maxElemSize() const { return sizeof(uint32_t) + 1; }

template <> inline size_t VByteDispatcher<uint32_t, true, false>::maxElemSize() const { return sizeof(uint32_t) + 1; }

template <> inline size_t VByteDispatcher<uint64_t, false, false>::maxElemSize() const { return sizeof(uint64_t); }

template <> inline size_t VByteDispatcher<uint64_t, true, true>::maxElemSize() const { return sizeof(uint64_t) + 2; }

template <> inline size_t VByteDispatcher<uint64_t, true, false>::maxElemSize() const { return sizeof(uint64_t) + 2; }

template <> inline size_t VByteDispatcher<uint8_t, false, false>::compressSingle(uint8_t in, uint8_t* out) const {
    memcpy(out, &in, sizeof(in));
    return sizeof(in);
}

template <> inline size_t VByteDispatcher<uint32_t, false, false>::compressSingle(uint32_t in, uint8_t* out) const {
    memcpy(out, &in, sizeof(in));
    return sizeof(in);
}

template <> inline size_t VByteDispatcher<uint32_t, true, true>::compressSingle(uint32_t in, uint8_t* out) const {
    return vbyte_compress_unsorted32(&in, out, 1);
}

template <> inline size_t VByteDispatcher<uint32_t, true, false>::compressSingle(uint32_t in, uint8_t* out) const {
    return vbyte_compress_unsorted32(&in, out, 1);
}

template <> inline size_t VByteDispatcher<uint64_t, false, false>::compressSingle(uint64_t in, uint8_t* out) const {
    memcpy(out, &in, sizeof(in));
    return sizeof(in);
}

template <> inline size_t VByteDispatcher<uint64_t, true, true>::compressSingle(uint64_t in, uint8_t* out) const {
    return vbyte_compress_unsorted64(&in, out, 1);
}

template <> inline size_t VByteDispatcher<uint64_t, true, false>::compressSingle(uint64_t in, uint8_t* out) const {
    return vbyte_compress_unsorted64(&in, out, 1);
}

template <>
inline size_t VByteDispatcher<uint8_t, false, false>::decompressSingle(const uint8_t* in, uint8_t* out) const {
    *out = *in;
    return sizeof(uint8_t);
}

template <>
inline size_t VByteDispatcher<uint32_t, false, false>::decompressSingle(const uint8_t* in, uint32_t* out) const {
    *out = *(reinterpret_cast<const uint32_t*>(in));
    return sizeof(uint32_t);
}

template <>
inline size_t VByteDispatcher<uint32_t, true, true>::decompressSingle(const uint8_t* in, uint32_t* out) const {
    return vbyte_uncompress_unsorted32(in, out, 1);
}

template <>
inline size_t VByteDispatcher<uint32_t, true, false>::decompressSingle(const uint8_t* in, uint32_t* out) const {
    return vbyte_uncompress_unsorted32(in, out, 1);
}

template <>
inline size_t VByteDispatcher<uint64_t, false, false>::decompressSingle(const uint8_t* in, uint64_t* out) const {
    *out = *(reinterpret_cast<const uint64_t*>(in));
    return sizeof(uint64_t);
}

template <>
inline size_t VByteDispatcher<uint64_t, true, true>::decompressSingle(const uint8_t* in, uint64_t* out) const {
    return vbyte_uncompress_unsorted64(in, out, 1);
}

template <>
inline size_t VByteDispatcher<uint64_t, true, false>::decompressSingle(const uint8_t* in, uint64_t* out) const {
    return vbyte_uncompress_unsorted64(in, out, 1);
}

template <>
inline size_t
VByteDispatcher<uint8_t, false, false>::compressData(const uint8_t* in, uint8_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint8_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint32_t, false, false>::compressData(const uint32_t* in, uint8_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint32_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, true>::compressData(const uint32_t* in, uint8_t* out, size_t length) const {
    return vbyte_compress_sorted32(in, out, 0, length);
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, false>::compressData(const uint32_t* in, uint8_t* out, size_t length) const {
    return vbyte_compress_unsorted32(in, out, length);
}

template <>
inline size_t
VByteDispatcher<uint64_t, false, false>::compressData(const uint64_t* in, uint8_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint64_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, true>::compressData(const uint64_t* in, uint8_t* out, size_t length) const {
    return vbyte_compress_sorted64(in, out, 0, length);
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, false>::compressData(const uint64_t* in, uint8_t* out, size_t length) const {
    return vbyte_compress_unsorted64(in, out, length);
}

template <>
inline size_t
VByteDispatcher<uint8_t, false, false>::decompressData(const uint8_t* in, uint8_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint8_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint32_t, false, false>::decompressData(const uint8_t* in, uint32_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint32_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, true>::decompressData(const uint8_t* in, uint32_t* out, size_t length) const {
    return vbyte_uncompress_sorted32(in, out, 0, length);
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, false>::decompressData(const uint8_t* in, uint32_t* out, size_t length) const {
    return vbyte_uncompress_unsorted32(in, out, length);
}

template <>
inline size_t
VByteDispatcher<uint64_t, false, false>::decompressData(const uint8_t* in, uint64_t* out, size_t length) const {
    const size_t byteCount = length * sizeof(uint64_t);
    memcpy(out, in, byteCount);
    return byteCount;
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, true>::decompressData(const uint8_t* in, uint64_t* out, size_t length) const {
    return vbyte_uncompress_sorted64(in, out, 0, length);
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, false>::decompressData(const uint8_t* in, uint64_t* out, size_t length) const {
    return vbyte_uncompress_unsorted64(in, out, length);
}

template <>
inline uint8_t VByteDispatcher<uint8_t, false, false>::selectData(const uint8_t* in, size_t size, size_t index) const {
    release_assert(index < size);
    return in[index];
}

template <>
inline uint32_t
VByteDispatcher<uint32_t, false, false>::selectData(const uint8_t* in, size_t size, size_t index) const {
    release_assert(index < size);
    return (reinterpret_cast<const uint32_t*>(in))[index];
}

template <>
inline uint32_t VByteDispatcher<uint32_t, true, true>::selectData(const uint8_t* in, size_t size, size_t index) const {
    return vbyte_select_sorted32(in, size, 0, index);
}

template <>
inline uint32_t VByteDispatcher<uint32_t, true, false>::selectData(const uint8_t* in, size_t size, size_t index) const {
    return vbyte_select_unsorted32(in, size, index);
}

template <>
inline uint64_t
VByteDispatcher<uint64_t, false, false>::selectData(const uint8_t* in, size_t size, size_t index) const {
    release_assert(index < size);
    return (reinterpret_cast<const uint64_t*>(in))[index];
}

template <>
inline uint64_t VByteDispatcher<uint64_t, true, true>::selectData(const uint8_t* in, size_t size, size_t index) const {
    return vbyte_select_sorted64(in, size, 0, index);
}

template <>
inline uint64_t VByteDispatcher<uint64_t, true, false>::selectData(const uint8_t* in, size_t size, size_t index) const {
    return vbyte_select_unsorted64(in, size, index);
}

template <>
inline size_t
VByteDispatcher<uint8_t, false, false>::searchData(const uint8_t* in, size_t length, uint8_t value) const {
    for (size_t i = 0; i < length; i++) {
        if (in[i] == value) {
            return i;
        }
    }
    return std::numeric_limits<size_t>::max();
}

template <>
inline size_t
VByteDispatcher<uint32_t, false, false>::searchData(const uint8_t* in, size_t length, uint32_t value) const {
    const uint32_t* base = reinterpret_cast<const uint32_t*>(in);
    for (size_t i = 0; i < length; i++) {
        if (base[i] == value) {
            return i;
        }
    }
    return std::numeric_limits<size_t>::max();
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, true>::searchData(const uint8_t* in, size_t length, uint32_t value) const {
    uint32_t actual = 0;
    const size_t rv = vbyte_search_lower_bound_sorted32(in, length, value, 0, &actual);
    if (actual != value || rv >= length) {
        return std::numeric_limits<size_t>::max();
    }
    return rv;
}

template <>
inline size_t
VByteDispatcher<uint32_t, true, false>::searchData(const uint8_t* in, size_t length, uint32_t value) const {
    const size_t rv = vbyte_search_unsorted32(in, length, value);
    if (rv >= length) {
        return std::numeric_limits<size_t>::max();
    }
    return rv;
}

template <>
inline size_t
VByteDispatcher<uint64_t, false, false>::searchData(const uint8_t* in, size_t length, uint64_t value) const {
    const uint64_t* base = reinterpret_cast<const uint64_t*>(in);
    for (size_t i = 0; i < length; i++) {
        if (base[i] == value) {
            return i;
        }
    }
    return std::numeric_limits<size_t>::max();
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, true>::searchData(const uint8_t* in, size_t length, uint64_t value) const {
    uint64_t actual = 0;
    const size_t rv = vbyte_search_lower_bound_sorted64(in, length, value, 0, &actual);
    if (actual != value || rv >= length) {
        return std::numeric_limits<size_t>::max();
    }
    return rv;
}

template <>
inline size_t
VByteDispatcher<uint64_t, true, false>::searchData(const uint8_t* in, size_t length, uint64_t value) const {
    const size_t rv = vbyte_search_unsorted64(in, length, value);
    if (rv >= length) {
        return std::numeric_limits<size_t>::max();
    }
    return rv;
}

/**
 * Immutable compact sparse row (CSR) storage for use with graph data.
 * Template parameters determine the underlying storage type for the edge vector (the larger of the two
 * vectors, typically by at least an order of magnitude), the integer type (uint32_t or uint64_t), whether
 * the data is encoded using libvbyte, and whether the edge data is sorted (per node). I say "edge vector",
 * but it is really any data that maps from a node to a list of integers.
 */
template <template <class> class EdgeVect, typename IType, bool Encoded, bool Sorted> class CSRStorageImm {
public:
    explicit CSRStorageImm()
        : m_numNodes(0),
          m_appendIndex(0),
          m_numValues(0),
          m_reverse(false) {}

    explicit CSRStorageImm(NodeIDSizeT numNodes, bool reverse = false)
        : m_numNodes(numNodes),
          m_appendIndex(0),
          m_numValues(0),
          m_reverse(reverse) {
        // We use m_appendIndex, so this can be resize() instead of reserve(). That way you
        // can query while constructing, because there will always be a next item (numNodes + 1)
        // that has an index <= yours.
        m_nodeIndexes.resize(numNodes + 1);
    }

    explicit CSRStorageImm(EagerFileVector<NodeIDSizeT> nodeIndexes,
                           EdgeVect<uint8_t> nodeBuckets,
                           size_t numValues,
                           bool reverse = false)
        : m_nodeIndexes(std::move(nodeIndexes)),
          m_nodeBuckets(std::move(nodeBuckets)),
          m_numNodes(m_nodeIndexes.size() - 1),
          m_appendIndex(m_numNodes),
          m_numValues(numValues),
          m_reverse(reverse) {}

    size_t numNodes() const { return m_numNodes; }
    size_t numValues() const { return m_numValues; }
    size_t numValueBytes() { return m_nodeBuckets.size(); }

    size_t numValuesAt(NodeID nodeId) {
        assert(nodeId < m_numNodes);
        const size_t csrIndex = m_reverse ? ((m_numNodes - nodeId) - 1) : nodeId;
        if (m_nodeIndexes[csrIndex] >= m_nodeIndexes[csrIndex + 1]) {
            return 0;
        }
        size_t maxLength = m_dispatch.maxElemSize();
        const uint8_t* indata = m_nodeBuckets.data(m_nodeIndexes[csrIndex], maxLength);
        IType count = 0;
        const IType ctBytes = m_dispatch.decompressSingle(indata, &count);
        assert(ctBytes > 0);
        return count;
    }

    /**
     * Get the list of data for the given node.
     *
     * @param[in] node The node to retrieve data for.
     * @param[out] outList The IType's will be appended to this list.
     * @return The number of "OutType" items added to outList.
     */
    size_t getData(NodeID node, std::vector<IType>& outList) {
        const size_t csrIndex = m_reverse ? ((m_numNodes - node) - 1) : node;
        if (csrIndex >= m_nodeIndexes.size()) {
            throw ApiMisuseFailure("Index out of bounds (CSR data)");
        }
        IType count = 0;
        const size_t bucket = m_nodeIndexes[csrIndex];
        if (bucket >= m_nodeIndexes[csrIndex + 1]) {
            return 0;
        }
        const size_t maxElemSize = m_dispatch.maxElemSize();
        size_t maxLength = maxElemSize;
        const uint8_t* indata = m_nodeBuckets.data(bucket, maxLength);
        const IType ctBytes = m_dispatch.decompressSingle(indata, &count);
        const size_t bytesNeeded = count * sizeof(IType);
        const size_t oldSize = outList.size();
        outList.resize(oldSize + count);
        maxLength = maxElemSize * count;
        const size_t bytesProcessed = m_dispatch.decompressData(
            m_nodeBuckets.data(bucket + ctBytes, maxLength), (IType*)(outList.data() + oldSize), count);
        release_assert(bytesProcessed <= maxLength);
        return count;
    }

    /**
     * Get a single datum for the given node.
     *
     * @param[in] node The node being modified.
     * @param[in] index The index.
     */
    IType getDatum(NodeID node, size_t index) {
        const size_t csrIndex = m_reverse ? ((m_numNodes - node) - 1) : node;
        if (csrIndex >= m_nodeIndexes.size()) {
            throw ApiMisuseFailure("Index out of bounds (CSR data)");
        }
        IType count = 0;
        const size_t bucket = m_nodeIndexes[csrIndex];
        if (bucket >= m_nodeIndexes[csrIndex + 1]) {
            return 0;
        }
        const size_t maxElemSize = m_dispatch.maxElemSize();
        size_t maxLength = maxElemSize;
        const uint8_t* indata = m_nodeBuckets.data(bucket, maxLength);
        const IType ctBytes = m_dispatch.decompressSingle(indata, &count);
        const size_t bytesNeeded = count * sizeof(IType);
        maxLength = maxElemSize * count;
        return m_dispatch.selectData(m_nodeBuckets.data(bucket + ctBytes, maxLength), count, index);
    }

    /**
     * Set data for the given node, using the list of IType given.
     *
     * @param[in] node The node being modified.
     * @param[in] fullData The list of IType to set.
     */
    size_t setData(NodeID node, const IType* fullData, const std::pair<size_t, size_t> range) {
        const size_t inputSize = range.second - range.first;
        if (inputSize == 0) {
            return 0;
        }
        const size_t csrIndex = m_reverse ? ((m_numNodes - node) - 1) : node;
        if (csrIndex < m_appendIndex || m_appendIndex == m_numNodes) {
            throw ApiMisuseFailure("Data is already set, and immutable");
        }
        const size_t curSize = m_nodeBuckets.size();
        while (m_appendIndex < csrIndex + 1) {
            m_nodeIndexes.ref(m_appendIndex++) = curSize;
        }
        const size_t maxElemSize = m_dispatch.maxElemSize();
        const size_t inputBytes = inputSize * sizeof(IType);
        // The +1 is for the count
        size_t storageByteEst = (inputSize + 1) * maxElemSize;
        const size_t bucket = curSize;
        m_nodeBuckets.resize(curSize + storageByteEst);
        uint8_t* output = (uint8_t*)m_nodeBuckets.dataRef(bucket, storageByteEst);
        const IType ctBytes = m_dispatch.compressSingle(inputSize, output);
        const size_t written =
            m_dispatch.compressData((const IType*)(fullData + range.first), output + ctBytes, inputSize);
        m_nodeBuckets.resize(curSize + written + ctBytes);
        m_numValues += inputSize;
        return written + ctBytes;
    }

    /**
     * Set data for the given node, using the list of IType given.
     *
     * @param[in] node The node being modified.
     * @param[in] dataList The list of IType to set.
     */
    size_t setData(NodeID node, const std::vector<IType>& dataList) {
        return setData(node, dataList.data(), {0, dataList.size()});
    }

    /**
     * Set data for the next node in order, using the list of IType given.
     *
     * @param[in] dataList The vector of IType to set.
     */
    size_t appendData(const std::vector<IType>& dataList) {
        return setData(m_appendIndex, dataList, {0, dataList.size()});
    }

    /**
     * Set data for the next node in order, using the list of IType given.
     *
     * @param[in] dataList The array of IType to set.
     * @param[in] length The length of the input array, in number of elements (not, e.g., bytes).
     */
    size_t appendData(const IType* dataList, size_t length) { return setData(m_appendIndex, dataList, {0, length}); }

    size_t search(NodeID node, IType value) {
        const size_t csrIndex = m_reverse ? ((m_numNodes - node) - 1) : node;
        if (csrIndex >= m_nodeIndexes.size()) {
            throw ApiMisuseFailure("Requested node not initialized in CSR data");
        }
        IType count = 0;
        const size_t bucket = m_nodeIndexes[csrIndex];
        if (bucket >= m_nodeIndexes[csrIndex + 1]) {
            throw ApiMisuseFailure("Index out of bounds for CSR data");
        }
        const size_t maxElemSize = m_dispatch.maxElemSize();
        size_t maxLength = maxElemSize;
        const uint8_t* indata = m_nodeBuckets.data(bucket, maxLength);
        const IType ctBytes = vbyte_uncompress_unsorted32(indata, &count, 1);
        maxLength = maxElemSize * count;
        return m_dispatch.searchData(m_nodeBuckets.data(bucket + ctBytes, maxLength), count, value);
    }

    /**
     * Finalize the last spot in the index so the calculations work properly. You can shrink the number
     * nodes in the CSR at this point if desired.
     * Can be safely called multiple times, it will only finalize the structure the first time.
     */
    void finalizeNodes(const NodeIDSizeT numNodes = 0) {
        if (numNodes != 0) {
            release_assert(numNodes <= m_numNodes);
            m_numNodes = numNodes;
        }
        if (m_appendIndex <= m_numNodes) {
            m_nodeIndexes.resize(m_numNodes + 1);
            while (m_appendIndex < m_numNodes + 1) {
                m_nodeIndexes.ref(m_appendIndex++) = m_nodeBuckets.size();
            }
        } else {
            release_assert(m_appendIndex == m_numNodes + 1);
        }
    }

    size_t flushBuckets(std::ostream& out, bool clearData = true) { return m_nodeBuckets.flush(out, !clearData); }

    size_t inmemoryBucketSize() const { return m_nodeBuckets.size(); }

    size_t flushIndexes(std::ostream& out, bool clearData = true) { return m_nodeIndexes.flush(out, !clearData); }

protected:
    VByteDispatcher<IType, Encoded, Sorted> m_dispatch;

    // For N nodes, this vector has N elements. Each element is an index into the buckets
    EagerFileVector<NodeIDSizeT> m_nodeIndexes;
    EdgeVect<uint8_t> m_nodeBuckets;
    size_t m_numNodes;
    size_t m_appendIndex;
    size_t m_numValues;
    bool m_reverse;
};

} // namespace grgl

#endif /* GRG_CSR_STORAGE_H */
