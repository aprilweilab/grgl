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
#ifndef GRGL_FILE_VECTOR_H
#define GRGL_FILE_VECTOR_H

#include "grgl/common.h"
#include "vbyte.h"
#include <fstream>
#include <ios>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <sys/types.h>
#include <vector>

namespace grgl {

using IFSPointer = std::shared_ptr<std::istream>;

/**
 * A vector (array/buffer) of data associated with a file and position within
 * that file. Just an interface.
 */
template <typename T> class FileVector {
public:
    FileVector(IFSPointer file, std::streamoff start)
        : m_file(std::move(file)),
          m_start(start) {}

    virtual size_t size() const = 0;
    virtual bool empty() const = 0;
    virtual void resize(size_t newSize, T defaultValue = {}) = 0;
    virtual void reserve(size_t numElements) = 0;
    virtual T& ref(size_t index) = 0;
    virtual T& atRef(size_t index) = 0;
    virtual const T& operator[](size_t index) = 0;
    virtual const T& at(size_t index) = 0;
    virtual const T* data(size_t startIndex, size_t& length) = 0;
    virtual T* dataRef(size_t startIndex, size_t& length) = 0;
    virtual void push_back(T item) = 0;

protected:
    IFSPointer m_file;
    std::streamoff m_start;
};

/**
 * An eagerly loaded file vector. Gets copied into RAM immediately, as a
 * std::vector under the hood. Does not hold onto a file handle, since it only needs
 * a one-time read operation.
 */
template <typename T> class EagerFileVector : public FileVector<T> {
public:
    EagerFileVector()
        : FileVector<T>({}, 0) {}

    EagerFileVector(const IFSPointer& file, std::streamoff start, size_t length)
        : FileVector<T>({}, start) {
        release_assert(file->good());
        file->seekg(start);
        m_data.resize(length);
        // We expect/require 8-byte alignment on the vector, for atomic operations below.
        release_assert((((size_t)(void*)m_data.data()) & 0x7U) == 0x0);
        file->read((char*)m_data.data(), sizeof(T) * m_data.size());
        const size_t expectedBytes = sizeof(T) * m_data.size();
        if (file->gcount() != expectedBytes) {
            std::stringstream ssErr;
            ssErr << "Expected to read " << expectedBytes << " bytes, but only read " << file->gcount();
            throw BadInputFileFailure(ssErr.str().c_str());
        }
    }

    virtual ~EagerFileVector() = default;
    EagerFileVector(EagerFileVector<T>&& rhs) noexcept = default;
    EagerFileVector(const EagerFileVector<T>& rhs) = delete;
    EagerFileVector<T>& operator=(EagerFileVector<T>&& rhs) noexcept = default;
    EagerFileVector<T>& operator=(const EagerFileVector<T>& rhs) = delete;

    bool empty() const override { return m_data.empty() && m_flushedItems == 0; }
    size_t size() const override { return m_data.size() + m_flushedItems; }

    void resize(size_t newSize, T defaultValue = {}) override {
        api_exc_check(newSize >= m_flushedItems, "Cannot shrink EagerFileVector beyond what has been flushed");
        m_data.resize(newSize - m_flushedItems, defaultValue);
    }

    void reserve(size_t numElements) override {
        api_exc_check(numElements >= m_flushedItems, "Cannot shrink EagerFileVector beyond what has been flushed");
        m_data.reserve(numElements - m_flushedItems);
    }

    /**
     * Store a value atomically to the vector.
     * On the architectures we support, just the proper size and alignment of the write will
     * make it atomic. If you are using this in a lock-free algorithm, you likely want to use
     * memory barriers as well (this method does not provide ordering, just atomicity).
     */
    void store_atomic(const size_t index, T value) {
        static_assert(sizeof(void*) >= 4, "Architecture problem");
        static_assert(sizeof(T) == 4 || sizeof(T) == sizeof(void*), "Data type must be atomically read/writeable");
        *((volatile T*)&m_data[index - m_flushedItems]) = value;
    }

    /**
     * Read a value atomically from the vector.
     * On the architectures we support, just the proper size and alignment of the read will
     * make it atomic. If you are using this in a lock-free algorithm, you likely want to use
     * memory barriers as well (this method does not provide ordering, just atomicity).
     */
    T read_atomic(const size_t index) {
        static_assert(sizeof(void*) >= 4, "Architecture problem");
        static_assert(sizeof(T) == 4 || sizeof(T) == sizeof(void*), "Data type must be atomically read/writeable");
        return *((volatile T*)&m_data[index - m_flushedItems]);
    }

    T& ref(size_t index) override {
        assert(index >= m_flushedItems && index < size());
        return m_data[index - m_flushedItems];
    }

    const T& cref(size_t index) const {
        assert(index >= m_flushedItems && index < size());
        return m_data[index - m_flushedItems];
    }

    T& atRef(size_t index) override {
        api_exc_check(index >= m_flushedItems, "Cannot access flushed items in EagerFileVector");
        index -= m_flushedItems;
        return m_data.at(index);
    }

    const T& operator[](size_t index) override {
        assert(index >= m_flushedItems && index < size());
        return m_data[index - m_flushedItems];
    }

    const T& at(size_t index) override {
        api_exc_check(index >= m_flushedItems, "Cannot access flushed items in EagerFileVector");
        index -= m_flushedItems;
        return m_data.at(index);
    }

    const T* data(size_t startIndex, size_t& length) override {
        api_exc_check(startIndex >= m_flushedItems, "Cannot access flushed items in EagerFileVector");
        startIndex -= m_flushedItems;
        if (length + startIndex > m_data.size()) {
            release_assert(startIndex < m_data.size());
            length = m_data.size() - startIndex;
        }
        return (m_data.data() + startIndex);
    }

    T* dataRef(size_t startIndex, size_t& length) override {
        api_exc_check(startIndex >= m_flushedItems, "Cannot access flushed items in EagerFileVector");
        startIndex -= m_flushedItems;
        if (length + startIndex > m_data.size()) {
            release_assert(startIndex < m_data.size());
            length = m_data.size() - startIndex;
        }
        return (m_data.data() + startIndex);
    }

    void push_back(T item) override { m_data.push_back(std::move(item)); }

    /**
     * Flush all the data in the vector so far to the given file descriptor.
     *
     * @return The number of items of type T that were flushed.
     */
    size_t flush(std::ostream& out, const bool keepInRam = false) {
        const size_t written = sizeof(T) * m_data.size();
        out.write(reinterpret_cast<const char*>(m_data.data()), written);
        if (!keepInRam) {
            m_flushedItems += m_data.size();
            m_data.clear();
        }
        return written / sizeof(T);
    }

    const std::vector<T>& vector() const { return m_data; }

protected:
    std::vector<T> m_data;
    size_t m_flushedItems{};
};

/**
 * A lazy version of FileVector that only buffers a small amount of the data from disk
 * as it is used.
 */
template <typename T> class LazyFileVector : public FileVector<T> {
public:
    static constexpr size_t DEFAULT_BUFFER_AMOUNT = 1024 * 128; // 128KB

    LazyFileVector(IFSPointer file, std::streamoff start, size_t size, size_t bufferAmount = DEFAULT_BUFFER_AMOUNT)
        : FileVector<T>(file, start),
          m_fileOffset(0),
          m_size(size),
          m_bufferAmount(bufferAmount) {
        release_assert(this->m_file->good());
        const size_t readAmount = std::min(size, m_bufferAmount);
        this->m_file->seekg(start);
        m_buffer.resize(readAmount);
        const size_t expectedBytes = sizeof(T) * m_buffer.size();
        this->m_file->read((char*)m_buffer.data(), expectedBytes);
        if (file->gcount() != expectedBytes) {
            std::stringstream ssErr;
            ssErr << "Expected to read " << expectedBytes << " bytes, but only read " << file->gcount();
            throw BadInputFileFailure(ssErr.str().c_str());
        }
    }

    size_t size() const override { return m_size; }

    bool empty() const override { return m_size == 0; }

    void resize(size_t newSize, T defaultValue = {}) override {
        api_exc_check(false, "Resizing not supported on LazyFileVector");
    }

    void reserve(size_t numElements) override { api_exc_check(false, "Resizing not supported on LazyFileVector"); }

    T& ref(size_t index) override {
        size_t singleItem = 1;
        const size_t offsetInBuffer = updateBufferIfNeeded(index, singleItem);
        return m_buffer[offsetInBuffer];
    }

    T& atRef(size_t index) override {
        size_t singleItem = 1;
        const size_t offsetInBuffer = updateBufferIfNeeded(index, singleItem);
        release_assert(singleItem == 1);
        return m_buffer[offsetInBuffer];
    }

    const T& operator[](size_t index) override {
        size_t singleItem = 1;
        const size_t offsetInBuffer = updateBufferIfNeeded(index, singleItem);
        return m_buffer[offsetInBuffer];
    }

    const T& at(size_t index) override {
        size_t singleItem = 1;
        const size_t offsetInBuffer = updateBufferIfNeeded(index, singleItem);
        release_assert(singleItem == 1);
        return m_buffer[offsetInBuffer];
    }

    const T* data(size_t startIndex, size_t& length) override {
        const size_t offsetInBuffer = updateBufferIfNeeded(startIndex, length);
        return &m_buffer[offsetInBuffer];
    }

    T* dataRef(size_t startIndex, size_t& length) override {
        const size_t offsetInBuffer = updateBufferIfNeeded(startIndex, length);
        return &m_buffer[offsetInBuffer];
    }
    void push_back(T item) override { api_exc_check(false, "Appending not supported on LazyFileVector"); }

protected:
    inline size_t updateBufferIfNeeded(size_t index, size_t& length) {
        release_assert(index < m_size);
        if (m_fileOffset <= index) {
            const size_t bufferedAmount = m_buffer.size();
            const size_t offsetFromBuffer = index - m_fileOffset;
            // The most common case (hopefully!) the data is entirely within the buffer.
            if (offsetFromBuffer + length <= bufferedAmount) {
                return offsetFromBuffer;
            }
        }
        // The data we need overlaps with the buffer. We don't try to optimize this case
        // (but rather avoid it, on average). Just seek and read!
        const size_t maxRemaining = m_size - index;
        const size_t readAmount = std::min(std::max(length, m_bufferAmount), maxRemaining);
        this->m_buffer.resize(readAmount);
        this->m_file->seekg(this->m_start + index);
        m_fileOffset = index;
        this->m_file->read((char*)m_buffer.data(), sizeof(T) * readAmount);
        length = std::min(this->m_file->gcount() / sizeof(T), length);
        return 0;
    }

    // Buffer containing data that has been read off disk.
    std::vector<T> m_buffer;
    // The offset (in T units) into the vector on file that we have read/seeked to so far. I.e.,
    // if this is 0 then the actual stream offset in the entire file (not just this vector) is
    // m_start. If this is C, then we are C*sizeof(T) bytes into the vector, and m_start+(C*sizeof(T))
    // bytes into the file.
    size_t m_fileOffset;
    // Number of sizeof(T) items that the vector contains.
    size_t m_size;
    // The number of items we should read from disk at a time, into m_buffer. If client code ever
    // requests more items than this (via data()) then we will fill the buffer with those items.
    const size_t m_bufferAmount;
};

} // namespace grgl

#endif
