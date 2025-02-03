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
#ifndef GRGL_CALCULATIONS_H
#define GRGL_CALCULATIONS_H

#include "grgl/grg.h"
#include <iosfwd>

static constexpr char const* USE_RANDOM_PHENOTYPE = "<<random>>";

void emitAlleleFrequency(grgl::GRGPtr& grg,
                         std::ostream& outStream,
                         std::pair<uint32_t, uint32_t> bpRange,
                         const grgl::NodeIDList& onlySamples);

void emitZygosityInfo(grgl::GRGPtr& grg,
                      std::ostream& outStream,
                      std::pair<uint32_t, uint32_t> bpRange,
                      const grgl::NodeIDList& onlySamples);

void emitBeta(const grgl::GRGPtr& grg, const std::string& phenotype, std::ostream& outStream, bool betaOnly = false);

#endif /* GRGL_CALCULATIONS_H */
