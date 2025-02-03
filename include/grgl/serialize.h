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
using IFSPointer = std::shared_ptr<std::istream>;

/**
 * Serialize the GRG to the given outstream.
 *
 * @param[in] grg The GRG to be serialized.
 * @param[in] out The (binary) output stream.
 */
std::pair<NodeIDSizeT, size_t> writeGrg(const GRGPtr& grg, std::ostream& out, bool allowSimplify = true);

class GRGOutputFilter;

std::pair<NodeIDSizeT, size_t> simplifyAndSerialize(const GRGPtr& grg,
                                                    std::ostream& outStream,
                                                    const GRGOutputFilter& filter,
                                                    bool allowSimplify = true);

/**
 * Deserialize the GRG from the given input stream.
 *
 * @param[in] inStream The (binary) input stream.
 */
MutableGRGPtr readMutableGrg(IFSPointer& inStream);
GRGPtr readImmutableGrg(IFSPointer& inStream, bool loadUpEdges = true);

}; // namespace grgl

#endif /* GRG_SERIALIZE_H */
