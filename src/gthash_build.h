#ifndef GRG_GTHASH_BUILD_H
#define GRG_GTHASH_BUILD_H

#include <memory>

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "gthash_index.h"
#include "util.h"

#define POSITION_KMER_SIZE (10)

namespace grgl {

class MutableGRG;
using MutableGRGPtr = std::shared_ptr<MutableGRG>;

/**
 * Given a VCF file, and a genome range, construct a GRG for that range - but do not map the
 * mutations. The resulting graph just has empty nodes and sample nodes, and is connected.
 */
MutableGRGPtr createEmptyGRGFromSamples(const std::string& sampleFile,
                                        FloatRange& genomeRange,
                                        size_t bitsPerMutation,
                                        bool useBinaryMuts,
                                        bool emitMissingData,
                                        bool flipRefMajor,
                                        double dropBelowThreshold);

}

#endif /* GRG_GTHASH_BUILD_H */
