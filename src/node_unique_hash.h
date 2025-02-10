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
#ifndef GRG_NODE_UNIQUE_HASH_H
#define GRG_NODE_UNIQUE_HASH_H

#include "grgl/grgnode.h"

#include <string>
#include <unordered_map>

#include "picohash.h"

namespace grgl {

using HashDigest = std::string;
using DigestToNode = std::unordered_map<HashDigest, NodeID>;

template <typename OrderedContainer>
inline HashDigest hashNodeSet(const OrderedContainer& nodeIdSet, size_t length = PICOHASH_MD5_DIGEST_LENGTH) {
    HashDigest result;
    result.resize(PICOHASH_MD5_DIGEST_LENGTH);
    picohash_ctx_t hashContext;
    picohash_init_md5(&hashContext);
    for (const auto nodeId : nodeIdSet) {
        picohash_update(&hashContext, &nodeId, sizeof(nodeId));
    }
    picohash_final(&hashContext, &result.front());
    if (length < PICOHASH_MD5_DIGEST_LENGTH) {
        return std::move(result.substr(length));
    }
    return std::move(result);
}

} // namespace grgl

#endif
