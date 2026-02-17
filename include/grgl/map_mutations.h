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
#ifndef MAP_MUTATIONS_H
#define MAP_MUTATIONS_H

#include "grgl/common.h"
#include "grgl/mut_iterator.h"

#include <algorithm>
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace grgl {

class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

class BitVecCoverageSet;
class SparseCoverageSet;

// BatchMembership is a lightweight view over externally owned membership bits for a mutation batch.
// Membership bits are stored separately as uint32_t words.
class BatchMembership {
public:
    BatchMembership() = default;
    BatchMembership(uint32_t* ptr, size_t batchSizeBits)
        : m_words32(ptr),
          m_numBits(batchSizeBits) {}

    size_t numWords32() const { return (m_numBits + 31) / 32; }

    void setOnes(size_t mutBatchBits) {
        release_assert(mutBatchBits <= m_numBits);
        release_assert(m_words32 != nullptr);

        if (mutBatchBits == 0) {
            clear();
            return;
        }

        const size_t wordCount = (mutBatchBits + 31) / 32;

        std::fill(m_words32, m_words32 + wordCount, 0xFFFFFFFFU);

        const size_t rem = mutBatchBits % 32;
        if (rem != 0) {
            m_words32[wordCount - 1] &= ((1U << rem) - 1U);
        }

        std::fill(m_words32 + wordCount, m_words32 + numWords32(), 0U);
    }

    void clear() {
        release_assert(m_words32 != nullptr);
        std::fill(m_words32, m_words32 + numWords32(), 0U);
    }

    bool any() const {
        return std::any_of(m_words32, m_words32 + numWords32(), [](uint32_t word) { return word != 0U; });
    }

    void andAssign(const BatchMembership& other) {
        release_assert(numWords32() == other.numWords32());
        for (size_t i = 0; i < numWords32(); ++i) {
            m_words32[i] &= other.m_words32[i];
        }
    }

    void setBit(size_t bitIndex) {
        release_assert(bitIndex < m_numBits);
        const size_t word = bitIndex / 32;
        const size_t off = bitIndex % 32;
        m_words32[word] |= (1U << off);
    }

    void copyFrom(const BatchMembership& other) {
        release_assert(numWords32() == other.numWords32());
        std::copy(other.m_words32, other.m_words32 + numWords32(), m_words32);
    }

    class Iterator {
    public:
        Iterator(const BatchMembership* left, const BatchMembership* right, bool atEnd)
            : m_left(left),
              m_right(right),
              m_wordIndex(atEnd ? left->numWords32() : 0),
              m_currentWordIndex(atEnd ? left->numWords32() : 0) {
            if (!atEnd) {
                ++(*this);
            }
        }

        size_t operator*() const {
            const int bit = __builtin_ctz(m_currentWord);
            return (m_currentWordIndex * 32) + static_cast<size_t>(bit);
        }

        Iterator& operator++() {
            // consume current bit if present
            if (m_currentWord != 0U) {
                m_currentWord &= (m_currentWord - 1U);
                if (m_currentWord != 0U) {
                    return *this;
                }
            }

            // find next nonzero word (or XOR-diff word)
            while (m_wordIndex < m_left->numWords32()) {
                uint32_t word = m_left->m_words32[m_wordIndex];
                if (m_right != nullptr) {
                    word ^= m_right->m_words32[m_wordIndex];
                }
                m_currentWordIndex = m_wordIndex;
                ++m_wordIndex;
                m_currentWord = word;
                if (m_currentWord != 0U) {
                    return *this;
                }
            }

            m_currentWordIndex = m_left->numWords32();
            m_currentWord = 0U;
            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return m_currentWordIndex != other.m_currentWordIndex || m_currentWord != other.m_currentWord;
        }

    private:
        const BatchMembership* m_left;
        const BatchMembership* m_right; // nullptr => bits(), non-null => XOR
        size_t m_wordIndex;
        size_t m_currentWordIndex;
        uint32_t m_currentWord{};
    };

    Iterator begin() const { return {this, nullptr, false}; }
    Iterator end() const { return {this, nullptr, true}; }

    Iterator beginDiff(const BatchMembership& other) const { return {this, &other, false}; }
    Iterator endDiff(const BatchMembership& other) const { return {this, &other, true}; }

private:
    uint32_t* m_words32{};
    size_t m_numBits = 0;
};

struct MutationMappingStats {
    size_t totalMutations{};
    size_t emptyMutations{};
    size_t mutationsWithOneSample{};
    size_t mutationsWithNoCandidates{};
    size_t reusedNodes{};
    size_t reusedNodeCoverage{};
    size_t reusedExactly{};
    size_t singletonSampleEdges{}; // Direct edge from Mutation to sample
    size_t newTreeNodes{};
    size_t samplesProcessed{};
    size_t numCandidates{};
    size_t reuseSizeBiggerThanHistMax{};
    size_t numWithSingletons{};
    size_t maxSingletons{};
    size_t reusedMutNodes{};
    std::vector<size_t> reuseSizeHist;

    void print(std::ostream& outStream) const {
        outStream << "mutations: " << this->totalMutations << std::endl;
        outStream << "candidates: " << this->numCandidates << std::endl;
        outStream << "emptyMutations: " << this->emptyMutations << std::endl;
        outStream << "mutationsWithOneSample: " << this->mutationsWithOneSample << std::endl;
        outStream << "singletonSampleEdges: " << this->singletonSampleEdges << std::endl;
        outStream << "samplesProcessed: " << this->samplesProcessed << std::endl;
        outStream << "reusedNodes: " << this->reusedNodes << std::endl;
        outStream << "reusedExactly: " << this->reusedExactly << std::endl;
        outStream << "reusedNodeCoverage: " << this->reusedNodeCoverage << std::endl;
        outStream << "reusedMutNodes: " << this->reusedMutNodes << std::endl;
        outStream << "reuseSizeBiggerThanHistMax: " << this->reuseSizeBiggerThanHistMax << std::endl;
        outStream << "numWithSingletons: " << this->numWithSingletons << std::endl;
        outStream << "maxSingletons: " << this->maxSingletons << std::endl;
        outStream << "avgSingletons: " << (double)this->singletonSampleEdges / (double)this->numWithSingletons
                  << std::endl;
    }
};

class HaplotypeIndex;

MutationMappingStats mapMutations(const MutableGRGPtr& grg,
                                  MutationIterator& mutations,
                                  bool verbose = false,
                                  size_t mutationBatchSize = 64);

MutationMappingStats mapMutations(const MutableGRGPtr& grg,
                                  const std::vector<Mutation>& mutations,
                                  const std::vector<NodeIDList>& samples,
                                  bool verbose = false,
                                  size_t mutationBatchSize = 64);

}; // namespace grgl

#endif /* MAP_MUTATIONS_H */
