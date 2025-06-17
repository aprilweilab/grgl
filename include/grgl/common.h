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
#ifndef GRGL_COMMON_H
#define GRGL_COMMON_H

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace grgl {

#define release_assert(condition)                                                                                      \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cerr << "release_assert(" #condition ") failed at " << __FILE__ << ":" << __LINE__ << std::endl;      \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)

/**
 * Exception thrown when the GRG API is used incorrectly (bad parameters, incorrect preconditions, etc.).
 */
class ApiMisuseFailure : public std::runtime_error {
public:
    explicit ApiMisuseFailure(char const* const message)
        : std::runtime_error(message) {}
};

#define api_exc_check(condition, msg)                                                                                  \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::stringstream ssErr;                                                                                   \
            ssErr << msg;                                                                                              \
            throw ApiMisuseFailure(ssErr.str().c_str());                                                               \
        }                                                                                                              \
    } while (0)

/**
 * Exception thrown when the input file does not meet conditions needed for GRG.
 */
class BadInputFileFailure : public std::runtime_error {
public:
    explicit BadInputFileFailure(char const* const message)
        : std::runtime_error(message) {}
};

class FloatRange {
public:
    static constexpr double RANGE_UNSPECIFIED = std::numeric_limits<double>::quiet_NaN();

    FloatRange()
        : m_start(RANGE_UNSPECIFIED),
          m_end(RANGE_UNSPECIFIED) {}

    FloatRange(const double start, const double end)
        : m_start(start),
          m_end(end) {}

    ~FloatRange() = default;
    FloatRange(const FloatRange&) = default;
    FloatRange(FloatRange&&) = default;
    FloatRange& operator=(const FloatRange&) = default;
    FloatRange& operator=(FloatRange&&) = default;

    double start() const { return m_start; }

    double end() const { return m_end; }

    inline bool isUnspecified() const { return std::isnan(m_start) || std::isnan(m_end); }

    bool isNormalized() const { return !isUnspecified() && m_end <= 1.0; }

    FloatRange normalized(size_t absStart, size_t absEnd) const {
        if (isUnspecified()) {
            return {};
        }
        // The +1 is because our range checking is all (inclusive, exclusive], but
        // we need our upper range to be 1.0, so we make sure that the max value is
        // <1.0
        const size_t span = (absEnd - absStart) + 1;
        if (!isNormalized() && span > 0) {
            return {(m_start - (double)absStart) / (double)span, (m_end - (double)absStart) / (double)span};
        }
        return *this;
    }

    FloatRange denormalized(size_t absStart, size_t absEnd) const {
        if (isUnspecified()) {
            return {};
        }
        if (isNormalized()) {
            const size_t span = (absEnd - absStart) + 1;
            return {(m_start * (double)span) + (double)absStart, (m_end * (double)span) + (double)absStart};
        }
        return *this;
    }

    bool contains(double position) const {
        assert(isUnspecified() || m_end > 1.0 || position <= 1.0);
        // An unspecified range includes everything.
        return isUnspecified() || (m_start <= position && m_end > position);
    }

private:
    double m_start;
    double m_end;
};

inline bool operator==(const FloatRange& a, const FloatRange& b) {
    return (a.isUnspecified() && b.isUnspecified()) || (a.start() == b.start() && a.end() == b.end());
}

// Ah, magic: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#header.boost.functional.hash_hpp
inline std::size_t hash_combine(std::size_t hash1, std::size_t hash2) {
    return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6U) + (hash1 >> 2U));
}

// Optional: set if the GRG has individual IDs, otherwise it is unset.
constexpr uint64_t GRG_FLAG_HAS_INDIV_IDS = 0x1;

constexpr uint64_t GRG_FILE_MAGIC = 0xE9366C64DDC8C5B0;
constexpr uint16_t GRG_FILE_MAJOR_VERSION = 5;
constexpr uint16_t GRG_FILE_MINOR_VERSION = 1;

#pragma pack(push, 1)
struct GRGFileHeader {
    uint64_t magic;
    uint16_t versionMajor;
    uint16_t versionMinor;
    uint16_t ploidy;
    uint16_t populationCount;
    uint64_t sampleCount;
    uint64_t mutationCount;
    uint64_t nodeCount;
    uint64_t edgeCount;
    uint64_t rangeStart; // in BP
    uint64_t rangeEnd;   // in BP
    uint64_t flags;
    uint64_t unused[7];
};
#pragma pack(pop)
static_assert(sizeof(GRGFileHeader) == 128, "GRG header size changed");

template <typename T> inline T roundUpToMultiple(const T input, const T multiple) {
    return ((input + (multiple - 1)) / multiple) * multiple;
}

} // namespace grgl

#endif /* GRGL_COMMON_H */
