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
#ifndef MAP_MUTATIONS_H
#define MAP_MUTATIONS_H

#include "grgl/common.h"
#include "grgl/mut_iterator.h"

#include <iosfwd>

namespace grgl {

class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

struct MutationMappingStats {
    size_t totalMutations{};
    size_t emptyMutations{};
    size_t mutationsWithOneSample{};
    size_t mutationsWithNoCandidates{};
    size_t reusedNodes{};
    size_t reusedNodeCoverage{};
    size_t reusedExactly{};
    size_t singletonSampleEdges{}; // Direct edge from Mutation to sample
    size_t newTreeNodes{};
    size_t samplesProcessed{};
    size_t numCandidates{};
    size_t reuseSizeBiggerThanHistMax{};
    size_t numWithSingletons{};
    size_t maxSingletons{};
    size_t reusedMutNodes{};
    std::vector<size_t> reuseSizeHist;

    void print(std::ostream& outStream) const {
        outStream << "mutations: " << this->totalMutations << std::endl;
        outStream << "candidates: " << this->numCandidates << std::endl;
        outStream << "emptyMutations: " << this->emptyMutations << std::endl;
        outStream << "mutationsWithOneSample: " << this->mutationsWithOneSample << std::endl;
        outStream << "singletonSampleEdges: " << this->singletonSampleEdges << std::endl;
        outStream << "samplesProcessed: " << this->samplesProcessed << std::endl;
        outStream << "reusedNodes: " << this->reusedNodes << std::endl;
        outStream << "reusedExactly: " << this->reusedExactly << std::endl;
        outStream << "reusedNodeCoverage: " << this->reusedNodeCoverage << std::endl;
        outStream << "reusedMutNodes: " << this->reusedMutNodes << std::endl;
        outStream << "reuseSizeBiggerThanHistMax: " << this->reuseSizeBiggerThanHistMax << std::endl;
        outStream << "numWithSingletons: " << this->numWithSingletons << std::endl;
        outStream << "maxSingletons: " << this->maxSingletons << std::endl;
        outStream << "avgSingletons: " << (double)this->singletonSampleEdges / (double)this->numWithSingletons
                  << std::endl;
    }
};

class HaplotypeIndex;

MutationMappingStats mapMutations(const MutableGRGPtr& grg, MutationIterator& mutations);

}; // namespace grgl

#endif /* MAP_MUTATIONS_H */
