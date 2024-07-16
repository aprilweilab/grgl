#ifndef GRGL_CALCULATIONS_H
#define GRGL_CALCULATIONS_H

#include <iosfwd>
#include "grgl/grg.h"

void emitAlleleFrequency(grgl::GRGPtr& grg,
                        std::ostream& outStream,
                        std::pair<uint32_t, uint32_t> bpRange,
                        const grgl::NodeIDList& onlySamples);

#endif /* GRGL_CALCULATIONS_H */