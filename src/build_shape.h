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
#ifndef GRG_BUILD_SHAPE_H
#define GRG_BUILD_SHAPE_H

#include <memory>

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "grgl/mutation.h"
#include "util.h"

namespace grgl {

class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

enum {
    GBF_EMPTY = 0x0U,
    GBF_NO_INDIVIDUAL_IDS = 0x1U,
    GBF_VERBOSE_OUTPUT = 0x2U,
};
using GrgBuildFlags = uint64_t;

/**
 * Given a VCF file, and a genome range, construct a GRG for that range - but do not map the
 * mutations. The resulting graph just has empty nodes and sample nodes, and is connected.
 */
MutableGRGPtr createEmptyGRGFromSamples(const std::string& sampleFile,
                                        FloatRange& genomeRange,
                                        size_t bitsPerMutation,
                                        GrgBuildFlags buildFlags,
                                        MutationIteratorFlags itFlags,
                                        double dropBelowThreshold,
                                        const std::map<std::string, std::string>& indivIdToPop,
                                        size_t tripletLevels = 0);

} // namespace grgl

#endif /* GRG_BUILD_SHAPE_H */
