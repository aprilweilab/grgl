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
#ifndef GRG_MUTATION_H
#define GRG_MUTATION_H

#include "grgl/common.h"

#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>

namespace grgl {

using EXTERNAL_ID = int32_t;
constexpr EXTERNAL_ID NO_ORIGINAL_ID = -1;

class Mutation;

using MutationId = uint32_t;
constexpr MutationId INVALID_MUTATION_ID = std::numeric_limits<MutationId>::max();
constexpr MutationId MAX_MUTATION_ID = INVALID_MUTATION_ID - 1;

using BpPosition = uint64_t;
constexpr BpPosition INVALID_POSITION = std::numeric_limits<BpPosition>::max();
constexpr BpPosition MAX_POSITION = INVALID_POSITION - 1;

/**
 * A mutation in the GRG.
 */
class Mutation {
public:
    static constexpr char const* ALLELE_A = "A";
    static constexpr char const* ALLELE_C = "C";
    static constexpr char const* ALLELE_T = "T";
    static constexpr char const* ALLELE_G = "G";
    static constexpr char const* ALLELE_MISSING = ".";
    static constexpr char const* ALLELE_0 = "0";
    static constexpr char const* ALLELE_1 = "1";

    Mutation()
        : m_alleleStorage({}),
          m_position(INVALID_POSITION),
          m_time(-1.0),
          m_originalId(NO_ORIGINAL_ID),
          m_refPosition(0) {}

    Mutation(BpPosition position, std::string allele)
        : m_alleleStorage({}),
          m_position(position),
          m_time(-1.0),
          m_originalId(NO_ORIGINAL_ID),
          m_refPosition(0) {
        if (allele.size() > std::numeric_limits<uint16_t>::max()) {
            throw ApiMisuseFailure("Allele length cannot exceed (2^16)-1");
        }
        setAlleleValues(std::move(allele), "");
    }

    Mutation(BpPosition position,
             std::string allele,
             const std::string& refAllele,
             float time = -1.0,
             EXTERNAL_ID originalId = NO_ORIGINAL_ID)
        : m_alleleStorage({}),
          m_position(position),
          m_time(time),
          m_originalId(originalId),
          m_refPosition(0) {
        if (allele.size() > std::numeric_limits<uint16_t>::max()) {
            throw ApiMisuseFailure("Allele length cannot exceed (2^16)-1");
        }
        setAlleleValues(std::move(allele), refAllele);
    }

    ~Mutation() { cleanup(); }

    void cleanup() {
        if (m_alleleStorage._str & PTR_STOLEN_BIT) {
            char* ptr = ((char*)(m_alleleStorage._str & ~PTR_STOLEN_BIT));
            delete[] ptr;
        }
        m_alleleStorage._str = 0;
    }

    Mutation(const Mutation& rhs)
        : m_alleleStorage({}),
          m_position(rhs.m_position),
          m_time(rhs.m_time),
          m_originalId(rhs.m_originalId),
          m_refPosition(0) {
        setAlleleValues(std::move(rhs.getAllele()), rhs.getRefAllele());
    }

    Mutation& operator=(const Mutation& rhs) {
        if (this != &rhs) {
            this->cleanup();
            this->m_position = rhs.m_position;
            this->m_time = rhs.m_time;
            this->m_originalId = rhs.m_originalId;
            this->m_refPosition = 0;
            setAlleleValues(std::move(rhs.getAllele()), rhs.getRefAllele());
        }
        return *this;
    }

    Mutation& operator=(Mutation&& rhs) noexcept {
        if (this != &rhs) {
            this->cleanup();
            this->m_position = rhs.m_position;
            this->m_time = rhs.m_time;
            this->m_originalId = rhs.m_originalId;
            this->m_refPosition = 0;
            setAlleleValues(std::move(rhs.getAllele()), rhs.getRefAllele());
        }
        return *this;
    }

    Mutation(Mutation&& rhs) noexcept
        : m_alleleStorage({}),
          m_position(rhs.m_position),
          m_time(rhs.m_time),
          m_originalId(rhs.m_originalId),
          m_refPosition(0) {
        setAlleleValues(std::move(rhs.getAllele()), rhs.getRefAllele());
    }

    bool isEmpty() const { return m_position == INVALID_POSITION; }

    BpPosition getPosition() const { return m_position; }

    float getTime() const { return m_time; }

    void setTime(float time) { m_time = time; }

    std::string getAllele() const {
        std::string storage = alleleStorageString();
        return storage.substr(0, m_refPosition);
    }

    std::string getRefAllele() const {
        std::string storage = alleleStorageString();
        return storage.substr(m_refPosition);
    }

    EXTERNAL_ID getOriginalId() const { return m_originalId; }

    bool operator==(const Mutation& rhs) const {
        if (this->m_position != rhs.m_position) {
            return false;
        }
        std::string storage = alleleStorageString();
        std::string storageRhs = rhs.alleleStorageString();
        return storage == storageRhs && this->m_originalId == rhs.m_originalId;
    }

    bool operator<(const Mutation& rhs) const {
        if (this->m_position == rhs.m_position) {
            std::string storage = alleleStorageString();
            std::string storageRhs = rhs.alleleStorageString();
            if (storage == storageRhs) {
                return this->m_originalId < rhs.m_originalId;
            }
            return storage < storageRhs;
        }
        return this->m_position < rhs.m_position;
    }

protected:
    void setAlleleValues(std::string altAllele, std::string refAllele) {
        const size_t totalBytes = altAllele.size() + refAllele.size();
        m_refPosition = altAllele.size();
        if (totalBytes > 7) {
            m_alleleStorage._str = (size_t) new char[totalBytes + 1];
            std::memcpy((void*)m_alleleStorage._str, altAllele.c_str(), altAllele.size());
            std::memcpy((void*)(m_alleleStorage._str + m_refPosition), refAllele.c_str(), refAllele.size());
            *(char*)(m_alleleStorage._str + totalBytes) = 0; // null term
            m_alleleStorage._str |= PTR_STOLEN_BIT;
        } else {
            for (size_t i = 0; i < altAllele.size(); i++) {
                m_alleleStorage._b[i + 1] = altAllele[i];
            }
            for (size_t i = altAllele.size(); i < totalBytes; i++) {
                m_alleleStorage._b[i + 1] = refAllele[i - altAllele.size()];
            }
        }
    }

    std::string alleleStorageString() const {
        std::string result;
        if (m_alleleStorage._str & PTR_STOLEN_BIT) {
            result = std::string((char*)(m_alleleStorage._str & ~PTR_STOLEN_BIT));
        } else {
            for (size_t i = 1; i < sizeof(m_alleleStorage); i++) {
                if (0 == m_alleleStorage._b[i]) {
                    break;
                }
                result.push_back(m_alleleStorage._b[i]);
            }
        }
        return std::move(result);
    }

    // It is safe to steal the least-significant bit from a heap-allocated pointer
    // since those will always be at least 4-byte aligned. Technically we could steal 2 bits...
    static constexpr size_t PTR_STOLEN_BIT = 0x1;
    // Long alleles use _str, short alleles (<=7 characters) use _b. We never use _b[0], because
    // the bit we can steal from the pointer is 0x1 which would conflict with character values.
    union {
        size_t _str;
        char _b[8];
    } m_alleleStorage; // 8
    static_assert(sizeof(m_alleleStorage) == 8, "Allele storage size changed");
    BpPosition m_position;    // 8
    float m_time;             // 4
    EXTERNAL_ID m_originalId; // 4
    uint16_t m_refPosition;   // 2
    // Have 6 spare bits here due to alignment
};
static_assert(sizeof(Mutation) == 32, "Mutation size changed");

/**
 * std::less-style comparison for Mutations that only looks at position and allele (ignores
 * all other fields).
 * Typically used for correctness comparisons.
 */
struct MutationLtPosAllele {
    bool operator()(const Mutation& lhs, const Mutation& rhs) const {
        if (lhs.getPosition() == rhs.getPosition()) {
            return lhs.getAllele() < rhs.getAllele();
        }
        return lhs.getPosition() < rhs.getPosition();
    }
};

} // namespace grgl

/**
 * Implement std::hash that uses only position/allele
 */
template <> struct std::hash<grgl::Mutation> {
    std::size_t operator()(const grgl::Mutation& mut) const noexcept {
        std::size_t h1 = std::hash<grgl::BpPosition>{}(mut.getPosition());
        std::size_t h2 = std::hash<std::string>{}(mut.getAllele());
        return grgl::hash_combine(h1, h2);
    }
};

#endif /* GRG_MUTATION_H */
