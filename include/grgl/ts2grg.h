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
 *
 * @param[in] treeSeq The tskit TreeSequence object.
 * @param[in] binaryMutations Set to true to emit Mutations as 0 (ref) and 1 (any alt).
 * @param[in] useNodeTimes Set to true to set the Mutation's time based on its tskit node time,
 *      instead of its tskit mutation time.
 * @param[in] maintainTopology Set to false to save some time and space when converting to GRG, at
 *      the cost of not capturing all topology changes due to recombinations "unseen" by Mutations
 *      above them. WARNING: This causes inaccurate sample-to-mutation mapping when there are
 *      back mutations in the trees (nested mutations at the same site).
 * @param[in] computeCoals Set to true to compute the number of individuals that coalesce at each
 *      node. Required if you want to compute the _exact_ variance for each mutation when performing
 *      diploid genotype matrix operations (e.g., GWAS) later. However, you can often just use the
 *      binomial approximation to variance and leave this as false. Computing the coalescences is
 *      expensive on large TreeSequences.
 * @param[in] treeRange Optional range of tree indices. If provided, only trees within that range
 *      [first, last) are converted (first is inclusive, last is exclusive).
 */
MutableGRGPtr convertTreeSeqToGRG(const tsk_treeseq_t* treeSeq,
                                  bool binaryMutations = false,
                                  bool useNodeTimes = false,
                                  bool maintainTopology = true,
                                  bool computeCoals = false,
                                  std::pair<size_t, size_t> treeRange = {});

} // namespace grgl

#endif /* TS2GRG_H */
