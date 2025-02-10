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
#include "grgl/windowing.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/mutation.h"
#include "util.h"
#include <limits>
#include <sstream>

namespace grgl {

WindowList windowByBP(std::pair<BpPosition, BpPosition> bpRange, size_t bpPerWindow, size_t overlap) {
    if (overlap == 0) {
        overlap = 1;
    }
    const size_t effectiveBpPerWindow = bpPerWindow / overlap;
    const size_t numWindows = ((bpRange.second - bpRange.first) + (effectiveBpPerWindow - 1)) / effectiveBpPerWindow;
    if (numWindows > MAX_WINDOW_ID) {
        std::stringstream errStr;
        errStr << "Number of windows must be less than or equal to " << MAX_WINDOW_ID << " (bpPerWindow too small)";
        throw ApiMisuseFailure(errStr.str().c_str());
    }
    WindowList windows(numWindows);
    for (size_t i = 0; i < windows.size(); i++) {
        if (i > 0) {
            windows[i].begin = windows[i - 1].begin + effectiveBpPerWindow;
        } else {
            windows[i].begin = bpRange.first;
        }
        windows[i].end = windows[i].begin + bpPerWindow;
    }
    assert(windows[numWindows - 1].end >= bpRange.second);
    return std::move(windows);
}

static inline std::pair<BpPosition, double> parseRow(const std::pair<std::string, std::string>& row,
                                                     const size_t lineNum) {
    std::pair<BpPosition, double> result;

#define ERR_LINE(message)                                                                                              \
    do {                                                                                                               \
        std::stringstream ssErr;                                                                                       \
        ssErr << message << " at line " << lineNum;                                                                    \
        throw BadInputFileFailure(ssErr.str().c_str());                                                                \
    } while (0)

    if (!parseExactUint64(row.first, result.first)) {
        ERR_LINE("Invalid position value ");
    }
    if (!parseExactDouble(row.second, result.second)) {
        ERR_LINE("Invalid cumulative cM value ");
    }
    return result;
}

WindowList
windowByCM(std::pair<BpPosition, BpPosition> bpRange, const std::string& mapFile, double cmPerWindow, size_t overlap) {
    if (overlap == 0) {
        overlap = 1;
    }
    using MapType = std::vector<std::pair<std::string, std::string>>;
    const MapType stringPairs = loadMapFromTSV<MapType>(mapFile, "Position(bp)", "Map(cM)");

    const double cmTillNextWindow = cmPerWindow / (double)overlap;
    BpPosition prevPos = 0;
    double prevCm = 0.0;
    double startCm = -1.0;
    double stopCm = -1.0;
    size_t lineNum = 1;
    std::list<Window> pendingWindows;
    double nextWindowCm = std::numeric_limits<double>::max();
    WindowList windows;
    for (const auto& posCmPair : stringPairs) {
        const auto bpAndCm = parseRow(posCmPair, lineNum);
        if (bpAndCm.first >= bpRange.first && startCm == -1.0) {
            // How far into the current span does our start position go.
            const double x = (double)(bpRange.first - prevPos) / (double)(bpAndCm.first - prevPos);
            startCm = (x * (bpAndCm.second - prevCm)) + prevCm;
            nextWindowCm = startCm + cmTillNextWindow;
            pendingWindows.push_back({bpRange.first, INVALID_POSITION});
        }
        if (bpAndCm.first >= bpRange.second && stopCm == -1.0) {
            // How far into the current span does our end position go.
            const double x = (double)(bpRange.second - prevPos) / (double)(bpAndCm.first - prevPos);
            stopCm = (x * (bpAndCm.second - prevCm)) + prevCm;
            break;
        }
        while (bpAndCm.second >= nextWindowCm) {
            const double x = (double)(nextWindowCm - prevCm) / (double)(bpAndCm.second - prevCm);
            const BpPosition bpPos = (BpPosition)(x * (double)(bpAndCm.first - prevPos)) + prevPos;
            if (pendingWindows.size() == overlap) {
                Window w = pendingWindows.front();
                pendingWindows.pop_front();
                w.end = bpPos;
                windows.push_back(w);
            }
            pendingWindows.push_back({bpPos, INVALID_POSITION});
            assert(pendingWindows.size() <= overlap);
            nextWindowCm += cmTillNextWindow;
        }
        prevPos = bpAndCm.first;
        prevCm = bpAndCm.second;
        lineNum++;
    }
    if (stopCm == -1.0) {
        stopCm = prevCm;
    }
    while (!pendingWindows.empty()) {
        Window w = pendingWindows.front();
        pendingWindows.pop_front();
        w.end = bpRange.second;
        windows.push_back(w);
    }

    return std::move(windows);
}

std::vector<std::pair<MutationId, NodeID>> mutationNodePairsForWindow(const GRGPtr& grg, const Window& window) {
    std::vector<std::pair<MutationId, NodeID>> result;
    for (const auto& nodeAndMutId : grg->getNodeMutationPairs()) {
        const auto& mutation = grg->getMutationById(nodeAndMutId.second);
        if (mutation.getPosition() >= window.begin && mutation.getPosition() < window.end) {
            result.emplace_back(nodeAndMutId.second, nodeAndMutId.first);
        }
    }
    return std::move(result);
}

std::vector<MutationId> mutationsForWindow(const GRGPtr& grg, const Window& window) {
    std::vector<MutationId> result;
    // TODO: this can be sped up if the mutations are known to be ordered by ID, which is the case for
    // most GRGs (unless we are in the middle of building a GRG, which we can just reject that scenario).
    for (MutationId mutId = 0; mutId < grg->numMutations(); mutId++) {
        const auto& mutation = grg->getMutationById(mutId);
        if (mutation.getPosition() >= window.begin && mutation.getPosition() < window.end) {
            result.emplace_back(mutId);
        }
    }
    return std::move(result);
}

} // namespace grgl
