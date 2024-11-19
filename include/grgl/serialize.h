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
 * should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef GRG_SERIALIZE_H
#define GRG_SERIALIZE_H

#include "grgl/grgnode.h"
#include "grgl/visitor.h"
#include <iosfwd>
#include <memory>
#include <stdexcept>

namespace grgl {

/**
 * Exception thrown when a GRG fails to serialize/deserialize.
 */
class SerializationFailure : public std::runtime_error {
public:
    explicit SerializationFailure(char const* const message)
        : std::runtime_error(message) {}
};

class GRG;
class MutableGRG;
using GRGPtr = std::shared_ptr<GRG>;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

/**
 * Simple class to encapsulate the filtering of a GRG to shrink it.
 */
class GRGOutputFilter {
public:
    GRGOutputFilter()
        : direction(TraversalDirection::DIRECTION_UP),
          seedList() {}

    GRGOutputFilter(TraversalDirection dir, const NodeIDList& seeds)
        : direction(dir),
          seedList(seeds) {}

    bool isSpecified() const { return !seedList.empty(); }

    TraversalDirection direction;
    NodeIDList seedList;
    std::pair<BpPosition, BpPosition> bpRange;
};

/**
 * Serialize the GRG to the given outstream.
 *
 * @param[in] grg The GRG to be serialized.
 * @param[in] out The (binary) output stream.
 * @param[in] useVarInt (optional) `true` by default. Set to `false` to avoid using variable-sized
 *      integer encoding in the serialization. The only reason to do this is if you are using another
 *      library/program to read GRG files and it does not support variable integer encoding.
 */
std::pair<NodeIDSizeT, size_t>
writeGrg(const GRGPtr& grg, std::ostream& out, bool useVarInt = true, bool allowSimplify = true);

std::pair<NodeIDSizeT, size_t> simplifyAndSerialize(
    const GRGPtr& grg, std::ostream& outStream, const GRGOutputFilter& filter, bool useVarInt, bool allowSimplify);

/**
 * Deserialize the GRG from the given input stream.
 *
 * @param[in] inStream The (binary) input stream.
 */
MutableGRGPtr readMutableGrg(std::istream& inStream);
GRGPtr readImmutableGrg(std::istream& inStream, bool loadUpEdges = true, bool loadDownEdges = true);

}; // namespace grgl

#endif /* GRG_SERIALIZE_H */
