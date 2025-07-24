#ifndef GRG_HAP_SEGMENT_H
#define GRG_HAP_SEGMENT_H

#include "hap_helpers.h"
#include "hap_index.h"
#include <limits>

// A "hap segment" is a small segment of a haplotype that is fixed width and enumerable.
// For a given dataset, and a hap segment length K (variants), there are theoretically
// 2^K possible hap segments, but in reality there are way fewer than that, even in
// dataset with heterogenous populations.

// Implementation:
//  Each segment is represented by a number between 0...P, where P is the number of unique
//  hap segments. The hap segment itself is captured by a bitvector of length K. Then we
//  can do Hamming distance on two samples by comparing their list of hap segments:
//   if the segment number is the same, distance is 0 over that segment, otherwise lookup
//   the bitvector and do the usual hamming distance calculation on it.

// The distance calculations are faster this way, and we can store exact haplotypes with
// the same amount of RAM we previously used for approximate haplotypes. Another advantage
// is that this representation lends itself to a heuristic algorithm for phasing as well:
//   * Step 1: phase data within each window of length K, by minimizing the total number
//             of unique hap segments.
//   * Step 2: the hap segments are now fixed, but the pair of haplotypes in a diploid can
//             swap hap segments between themselves, in a way that minimizes distances to
//             the nearest neighbors of all haplotypes.

template <> struct std::hash<grgl::HaplotypeVector> {
    size_t operator()(const grgl::HaplotypeVector& idList) const noexcept {
        size_t hashValue = 42;
        for (const auto& id : idList) {
            hashValue = grgl::hash_combine(hashValue, std::hash<grgl::HapVectorT>{}(id));
        }
        return hashValue;
    }
};

namespace grgl {

using HapDist = uint8_t;
constexpr size_t INVALID_HAP_DIST = std::numeric_limits<HapDist>::max();
constexpr size_t MAX_FAST_HAP_LENGTH = INVALID_HAP_DIST - 1;

using HapIdx = uint32_t;
using HapIdxList = std::vector<HapIdx>;
constexpr HapIdx INVALID_HAP_IDX = std::numeric_limits<HapIdx>::max();
constexpr size_t MAX_FAST_HAPLOTYPES = INVALID_HAP_IDX - 1;
// The 0th hap index always represents the "zero" haplotype (that is, no mutations)
constexpr HapIdx ZERO_HAP_IDX = 0;

constexpr size_t HAP_INVALID_BIT_COUNT = std::numeric_limits<size_t>::max();

/**
 * Compressed representation of a single Window containing a bunch of "hap segments"
 * that are enumerated by uniqueness. These are forced to be short (<= 254 variants
 * within the window), and enumerated as 32-bit integers, which means they can compress
 * about 8x.
 */
struct HaplotypeWindow {
    size_t firstMutationIdx;
    size_t lastMutationIdx;
    std::vector<HaplotypeVector> orderedHaps;
    std::vector<size_t> orderedHapCounts;

    // Distance matrix between the first X haplotypes for the window. Where
    // X is stored in precomputed.
    // FIXME: this uses X^2 bytes, but it is symmetric so we only need (X^2)/2
    std::vector<HapDist> distMatrix;
    size_t precomputed;

    // All unique haplotype segments
    std::unordered_map<HaplotypeVector, HapIdx> hapMap;

    size_t addOrGetHaplotype(HaplotypeVector newHapVec, size_t bitCount = HAP_INVALID_BIT_COUNT) {
        size_t nextHapIdx = this->orderedHaps.size();
        const auto inserted = this->hapMap.emplace(newHapVec, nextHapIdx);
        if (inserted.second) {
            if (bitCount == HAP_INVALID_BIT_COUNT) {
                bitCount = countBits(newHapVec);
            }
            this->orderedHaps.push_back(std::move(newHapVec));
            this->orderedHapCounts.push_back(bitCount);
            assert(this->orderedHaps.size() == this->hapMap.size());
        } else {
            nextHapIdx = inserted.first->second;
        }
        return nextHapIdx;
    }

    // Force computation of the distance, even if it is already precomputed.
    inline size_t computeDistance(const HapIdx hapIdx1, const HapIdx hapIdx2) {
        if (hapIdx1 == ZERO_HAP_IDX) {
            return this->orderedHapCounts.at(hapIdx2);
        }
        if (hapIdx2 == ZERO_HAP_IDX) {
            return this->orderedHapCounts.at(hapIdx1);
        }
        const HaplotypeVector& h1_i = this->orderedHaps.at(hapIdx1);
        const HaplotypeVector& h2_i = this->orderedHaps.at(hapIdx2);
        return bitwiseHamming(h1_i, h2_i);
    };

    // Get the distance by either looking it up or computing it.
    inline size_t getDistance(const HapIdx hapIdx1, const HapIdx hapIdx2) {
        if (hapIdx1 == hapIdx2) {
            return 0;
        }
        if (hapIdx1 >= this->precomputed || hapIdx2 >= this->precomputed) {
            return this->computeDistance(hapIdx1, hapIdx2);
        }
        if (hapIdx1 > hapIdx2) {
            auto& h1h2Cell = this->distMatrix.at((hapIdx2 * this->precomputed) + hapIdx1);
            if (h1h2Cell == INVALID_HAP_DIST) {
                h1h2Cell = this->computeDistance(hapIdx1, hapIdx2);
            }
            return h1h2Cell;
        }
        if (hapIdx2 > hapIdx1) {
            auto& h1h2Cell = this->distMatrix.at((hapIdx1 * this->precomputed) + hapIdx2);
            if (h1h2Cell == INVALID_HAP_DIST) {
                h1h2Cell = this->computeDistance(hapIdx1, hapIdx2);
            }
            return h1h2Cell;
        }
        release_assert(false); // Unreachable.
    }
};

} // namespace grgl

#endif /* GRG_HAP_SEGMENT_H */
