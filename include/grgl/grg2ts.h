/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2026 April Wei
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
#ifndef GRG2TS_H
#define GRG2TS_H

#include <memory>
#include <tskit.h>

namespace grgl {

class GRG;
using GRGPtr = std::shared_ptr<GRG>;

/**
 * Convert a GRG into a tskit tree-sequence object.
 *
 * @param[in] grg The GRG object.
 * @param[in] treeRange Optional range of tree indices. If provided, only trees within that range
 *      [first, last) are converted (first is inclusive, last is exclusive).
 */
void convertGRGToTreeSeq(GRGPtr& grg, tsk_treeseq_t* outTS, std::pair<size_t, size_t> treeRange = {});

} // namespace grgl

#endif /* GRG2TS_H */
