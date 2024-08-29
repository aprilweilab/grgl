#ifndef GRG_BUILD_SHAPE_H
#define GRG_BUILD_SHAPE_H

#include <memory>

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "util.h"

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
                                        double dropBelowThreshold,
                                        const std::map<std::string, std::string>& indivIdToPop,
                                        size_t tripletLevels = 0);

}

#endif /* GRG_BUILD_SHAPE_H */