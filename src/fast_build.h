/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2025 April Wei
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
#ifndef GRG_FAST_BUILD_H
#define GRG_FAST_BUILD_H

#include <memory>

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "grgl/mutation.h"
#include "util.h"

namespace grgl {

class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

/* Flags that change how the GRG build operation works.
 */
enum {
    GBF_EMPTY = 0x0U,
    GBF_NO_INDIVIDUAL_IDS = 0x1U,
    GBF_VERBOSE_OUTPUT = 0x2U,
    GBF_NO_TREE_MAP = 0x4U,
    GBF_TREES_FASTER1 = 0x8U,  // Build more trees than optimal, for speed
    GBF_TREES_FASTER2 = 0x10U, // Build even more trees than optimal, for speed
};
using GrgBuildFlags = uint64_t;

// See TODO in fast_build.cpp
constexpr size_t FIXED_HAP_LENGTH = 254;

/**
 * Add the non-genotype information (if any) to the GRG, such as the base-pair range, individuals,
 * and populations.
 */
void addExtraInfoToGRG(MutableGRGPtr& grg,
                       grgl::MutationIterator& mutationIterator,
                       GrgBuildFlags buildFlags,
                       const std::map<std::string, std::string>& indivIdToPop);

/**
 * Given an input file, and a genome range, construct a GRG for that range that simultaneously
 * maps mutations.
 *
 * @param[in] sampleFile The input file (IGD, BGEN, or VCF).
 * @param[in] genomeRange The range to apply this to. MUST BE IN VARIANT COUNTS, NOT BASEPAIR.
 * @param[in] hapLength The number of variants per "haplotype": this determines the granularity of
 *      the algorithm and how fast/memory efficient it is. Larger == more coarse (larger graph in
 *      the end) but faster and less memory.
 * @param[in] buildFlags
 * @param[in] itFlags
 * @param[in] noTreeBuildThreshold Any variants with frequency/count below this threshold do not
 *      contribute to the tree building. If you are mapping mutations directly, then they will be
 *      mapped by creating a Mutation node and one edge per sample.
 * @param[in] indivIdToPop Map from individual identifier to population identifier.
 * @param[in] rebuildProportion When the nearest-neighbor index reached this proportion of deleted
 *      items, rebuild it. Default: 0.10.
 */
MutableGRGPtr fastGRGFromSamples(const std::string& filePrefix,
                                 const std::string& sampleFile,
                                 FloatRange& genomeRange,
                                 GrgBuildFlags buildFlags,
                                 MutationIteratorFlags itFlags,
                                 size_t treeCount,
                                 double noTreeBelowThreshold,
                                 const std::map<std::string, std::string>& indivIdToPop,
                                 double rebuildProportion = 0.10);

} // namespace grgl

#endif /* GRG_FAST_BUILD_H */
