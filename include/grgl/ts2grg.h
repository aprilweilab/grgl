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
#ifndef TS2GRG_H
#define TS2GRG_H

#include <functional>
#include <list>
#include <memory>
#include <stdexcept>
#include <tskit.h>

namespace grgl {

/**
 * Exception thrown when a call to Tskit fails unexpectedly.
 */
class TskitApiFailure : public std::runtime_error {
public:
    explicit TskitApiFailure(char const* const message)
        : std::runtime_error(message) {}
};

class Mutation;
class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

/**
 * Convert a tskit tree-sequence object into a GRG.
 */
MutableGRGPtr convertTreeSeqToGRG(const tsk_treeseq_t* treeSeq,
                                  bool binaryMutations = false,
                                  bool useNodeTimes = false,
                                  bool maintainTopology = false,
                                  bool computeCoals = false);

} // namespace grgl

#endif /* TS2GRG_H */
