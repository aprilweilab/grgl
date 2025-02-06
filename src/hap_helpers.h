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
#ifndef GRGL_HAP_HELPERS_H
#define GRGL_HAP_HELPERS_H

#include <algorithm>
#include <cstdint>
#include <list>
#include <vector>

#include "grgl/mutation.h"
#include "util.h"

namespace grgl {

using HapVectorT = uint32_t;
using HaplotypeVector = std::vector<HapVectorT>;
using MutationList = std::vector<Mutation>;

void dumpHash(const HaplotypeVector& hash);

inline bool bitwiseIsZero(const HaplotypeVector& gtHash) {
    return std::all_of(gtHash.cbegin(), gtHash.cend(), [](HapVectorT value) { return value == 0; });
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

inline HaplotypeVector bitwiseIntersect(const HaplotypeVector& vect1, const HaplotypeVector& vect2, size_t& bitsSet) {
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

inline size_t bitwiseHamming(const HaplotypeVector& hash1, const HaplotypeVector& hash2) {
    size_t dist = 0;
    release_assert(hash1.size() == hash2.size());
    for (size_t i = 0; i < hash1.size(); i++) {
        static_assert(sizeof(HapVectorT) == 4, "Optimization for counting set bits assumes 32-bit ints");
        const uint32_t xorValue = hash1[i] ^ hash2[i];
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

inline size_t bloomFilterCapacity(size_t numBits) { return (numBits + (sizeof(HapVectorT) - 1)) / sizeof(HapVectorT); }

} // namespace grgl

#endif /* GRGL_HAP_HELPERS_H */
