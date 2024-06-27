#ifndef GRG_NODE_UNIQUE_HASH_H
#define GRG_NODE_UNIQUE_HASH_H

#include "grgl/grgnode.h"

#include <unordered_map>
#include <string>

#include "picohash.h"

namespace grgl {

using HashDigest = std::string;
using DigestToNode = std::unordered_map<HashDigest, NodeID>;

template <typename OrderedContainer>
inline HashDigest hashNodeSet(
        const OrderedContainer& nodeIdSet,
        size_t length = PICOHASH_MD5_DIGEST_LENGTH) {
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

}

#endif