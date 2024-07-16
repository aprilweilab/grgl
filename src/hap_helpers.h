#ifndef GRGL_HAP_HELPERS_H
#define GRGL_HAP_HELPERS_H

#include <algorithm>
#include <cstdint>
#include <list>
#include <vector>

#include "grgl/mutation.h"
#include "similarity/mmh3.h"
#include "util.h"

namespace grgl {

using HapVectorT = uint32_t;
using HaplotypeVector = std::vector<HapVectorT>;
using MutationList = std::vector<Mutation>;

void dumpHash(const HaplotypeVector& hash);

inline bool bitwiseIsZero(const HaplotypeVector& gtHash) {
    return std::all_of(gtHash.cbegin(), gtHash.cend(), [](HapVectorT value) { return value == 0; });
}

inline HaplotypeVector bitwiseIntersect(const std::list<const HaplotypeVector*>& gtHashList) {
    release_assert(!gtHashList.empty());
    HaplotypeVector result = *gtHashList.front();
    for (const auto& gtHash : gtHashList) {
        const size_t minSize = std::min(result.size(), gtHash->size());
        for (size_t bucket = 0; bucket < minSize; bucket++) {
            result[bucket] &= (*gtHash)[bucket];
        }
        release_assert(minSize > 0);
        result.resize(minSize);
    }
    return std::move(result);
}

inline HaplotypeVector bitwiseUnion(std::list<HaplotypeVector>& gtHashList) {
    release_assert(gtHashList.size() > 1);
    HaplotypeVector result = gtHashList.front();
    for (const auto& gtHash : gtHashList) {
        release_assert(result.size() == gtHash.size());
        for (size_t bucket = 0; bucket < gtHash.size(); bucket++) {
            result[bucket] |= gtHash[bucket];
        }
    }
    return std::move(result);
}

inline uint32_t countBits(uint32_t value) {
    uint32_t numBits = value - ((value >> 1U) & 0x55555555U);
    numBits = (numBits & 0x33333333U) + ((numBits >> 2U) & 0x33333333U);
    return ((numBits + (numBits >> 4U) & 0xF0F0F0FU) * 0x1010101U) >> 24U;
}

inline uint32_t countBits(const HaplotypeVector& vect) {
    uint32_t count = 0;
    for (auto value : vect) {
        uint32_t numBits = value - ((value >> 1U) & 0x55555555U);
        numBits = (numBits & 0x33333333U) + ((numBits >> 2U) & 0x33333333U);
        count += ((numBits + (numBits >> 4U) & 0xF0F0F0FU) * 0x1010101U) >> 24U;
    }
    return count;
}

inline HaplotypeVector bitwiseIntersect(
        const HaplotypeVector& vect1,
        const HaplotypeVector& vect2,
        size_t& bitsSet) {
    bitsSet = 0;
    HaplotypeVector result = vect1;
    const size_t minSize = std::min(result.size(), vect2.size());
    for (size_t bucket = 0; bucket < minSize; bucket++) {
        result[bucket] &= vect2[bucket];
        bitsSet += countBits(result[bucket]);
    }
    release_assert(minSize > 0);
    result.resize(minSize);
    return std::move(result);
}

// Subtract second from first, modifying first. Return number of bits set in the result.
inline size_t bitwiseSubtract(HaplotypeVector& first, const HaplotypeVector& second) {
    size_t after = 0;
    const size_t minSize = std::min(first.size(), second.size());
    HaplotypeVector result(minSize);
    for (size_t bucket = 0; bucket < minSize; bucket++) {
        first[bucket] = first[bucket] & (first[bucket] ^ second[bucket]);
        after += countBits(first[bucket]);
    }
    return after;
}

inline size_t bitwiseHamming(const HaplotypeVector& hash1,
                             const HaplotypeVector& hash2,
                             size_t amount) {
    size_t dist = 0;
    for (size_t i = 0; i < amount; i++) {
        static_assert(sizeof(HapVectorT) == 4,
                      "Optimization for counting set bits assumes 32-bit ints");
        uint32_t xorValue = hash1[i] ^ hash2[i];
        // Fun times: https://graphics.stanford.edu/~seander/bithacks.htm
        uint32_t numBits = xorValue - ((xorValue >> 1U) & 0x55555555U);
        numBits = (numBits & 0x33333333U) + ((numBits >> 2U) & 0x33333333U);
        numBits = ((numBits + (numBits >> 4U) & 0xF0F0F0FU) * 0x1010101U) >> 24U;
        dist += numBits;
    }
    return dist;
}


inline size_t bitwiseHamming(const HaplotypeVector& hash1,
                             const HaplotypeVector& hash2) {
    size_t dist = 0;
    const size_t s1 = hash1.size();
    const size_t s2 = hash2.size();
    const size_t maxSize = std::max(s1, s2);
    for (size_t i = 0; i < maxSize; i++) {
        static_assert(sizeof(HapVectorT) == 4,
                      "Optimization for counting set bits assumes 32-bit ints");
        uint32_t xorValue;
        if (i >= s1) {
            xorValue = hash2[i];
        } else if (i >= s2) {
            xorValue = hash1[i];
        } else {
            xorValue = hash1[i] ^ hash2[i];
        }
        // Fun times: https://graphics.stanford.edu/~seander/bithacks.htm
        uint32_t numBits = xorValue - ((xorValue >> 1U) & 0x55555555U);
        numBits = (numBits & 0x33333333U) + ((numBits >> 2U) & 0x33333333U);
        numBits = ((numBits + (numBits >> 4U) & 0xF0F0F0FU) * 0x1010101U) >> 24U;
        dist += numBits;
    }
    return dist;
}

inline void setBit(HaplotypeVector& vect, const size_t bitIndex) {
    const size_t element = bitIndex / sizeof(HapVectorT);
    const HapVectorT mask = 0x1U << (bitIndex % sizeof(HapVectorT));
    vect[element] |= mask;
}

inline void setBitExpanding(HaplotypeVector& vect, const size_t bitIndex) {
    const size_t element = bitIndex / sizeof(HapVectorT);
    const HapVectorT mask = 0x1U << (bitIndex % sizeof(HapVectorT));
    if (element >= vect.size()) {
        vect.resize(element+1);
    }
    vect[element] |= mask;
}

inline size_t bloomFilterCapacity(size_t numBits) {
    return (numBits + (sizeof(HapVectorT)-1)) / sizeof(HapVectorT);
}

template <typename T>
void bloomFilterAddItem(HaplotypeVector& vect, const T& item) {
    const size_t K = vect.size() * sizeof(HapVectorT);
    const auto hashIdx = murmur3_32((const uint8_t*)&item, sizeof(item), 1);
    const size_t bitIndex = hashIdx % K;
    const size_t element = bitIndex / sizeof(HapVectorT);
    const HapVectorT mask = 0x1U << (bitIndex % sizeof(HapVectorT));
    vect[element] |= mask;
}

template <typename T>
void bloomFilterAddItem(HapVectorT* vect, const T& item, size_t vectSize) {
    const size_t K = vectSize * sizeof(HapVectorT);
    const auto hashIdx = murmur3_32((const uint8_t*)&item, sizeof(item), 1);
    const size_t bitIndex = hashIdx % K;
    const size_t element = bitIndex / sizeof(HapVectorT);
    const HapVectorT mask = 0x1U << (bitIndex % sizeof(HapVectorT));
    vect[element] |= mask;
}

inline size_t listwiseHamming(const std::vector<MutationId>& list1, const std::vector<MutationId>& list2) {
    size_t dist = 0;
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < list1.size() || i2 < list2.size()) {
        if (i1 < list1.size() && i2 < list2.size()) {
            const MutationId& id1 = list1[i1];
            const MutationId& id2 = list2[i2];
            if (id1 < id2) {
                i1++;
                dist++;
            } else if (id2 < id1) {
                i2++;
                dist++;
            } else {
                i1++; i2++;
            }
        } else if (i1 < list1.size()) {
            dist++;
            i1++;
        } else {
            dist++;
            i2++;
        }
    }
    return dist;
}

inline std::vector<MutationId>
listwiseIntersect(const std::vector<MutationId>& list1, const std::vector<MutationId>& list2) {
    std::vector<MutationId> result;
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < list1.size() || i2 < list2.size()) {
        if (i1 < list1.size() && i2 < list2.size()) {
            const MutationId& id1 = list1[i1];
            const MutationId& id2 = list2[i2];
            if (id1 < id2) {
                i1++;
            } else if (id2 < id1) {
                i2++;
            } else {
                result.emplace_back(id1);
                i1++; i2++;
            }
        } else if (i1 < list1.size()) {
            i1++;
        } else {
            i2++;
        }
    }
    return std::move(result);
}

inline std::vector<MutationId>
listwiseSubtract(const std::vector<MutationId>& list1, const std::vector<MutationId>& list2) {
    std::vector<MutationId> result;
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < list1.size() || i2 < list2.size()) {
        if (i1 < list1.size() && i2 < list2.size()) {
            const MutationId& id1 = list1[i1];
            const MutationId& id2 = list2[i2];
            if (id1 < id2) {
                result.emplace_back(id1);
                i1++;
            } else if (id2 < id1) {
                i2++;
            } else {
                i1++; i2++;
            }
        } else if (i1 < list1.size()) {
            result.emplace_back(list1[i1]);
            i1++;
        } else {
            i2++;
        }
    }
    return std::move(result);
}

}

#endif /* GRGL_HAP_HELPERS_H */