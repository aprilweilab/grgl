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
#ifndef GRG_WINDOWING_H
#define GRG_WINDOWING_H

#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

#include "grgl/grgnode.h"
#include "grgl/mutation.h"

namespace grgl {

// We support up to ~65,000 windows, which for human chromosome 2 is about 800bp
// per window on average.
using WindowId = uint16_t;
static constexpr WindowId MAX_WINDOW_ID = std::numeric_limits<WindowId>::max();

struct Window {
    BpPosition begin; // Inclusive.
    BpPosition end;   // Exclusive.
};

using WindowList = std::vector<Window>;

/**
 * Given a range, split it up into equal-sized windows.
 *
 * @param[in] bpRange The base-pair range (usually from GRG::getBPRange() or GRG::getSpecifiedBPRange()).
 * @param[in] numWindows How many windows to use total.
 * @param[in] overlap How many windows overlap. numWindows must be divisible by overlap.
 */
WindowList windowByBP(std::pair<BpPosition, BpPosition> bpRange, size_t bpPerWindow, size_t overlap);

/**
 * Given a GRG, split it up into windows of w cM each.
 */
WindowList
windowByCM(std::pair<BpPosition, BpPosition> bpRange, const std::string& mapFile, double cmPerWindow, size_t overlap);

class GRG;
using GRGPtr = std::shared_ptr<GRG>;

/**
 * Given a GRG and Window, produce a list of MutationIDs and their corresponding NodeIDs for those mutations,
 * which fall within the Window.
 */
std::vector<std::pair<MutationId, NodeID>> mutationNodePairsForWindow(const GRGPtr& grg, const Window& window);

/**
 * Given a GRG and Window, produce a list of MutationIDs which fall within the Window.
 */
std::vector<MutationId> mutationsForWindow(const GRGPtr& grg, const Window& window);

} // namespace grgl

#endif /* GRG_WINDOWING_H */
