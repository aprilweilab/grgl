#ifndef GRGL_CALCULATIONS_H
#define GRGL_CALCULATIONS_H

#include <iosfwd>
#include "grgl/grg.h"

static constexpr char const* USE_RANDOM_PHENOTYPE = "<<random>>";

void emitAlleleFrequency(grgl::GRGPtr& grg,
                         std::ostream& outStream,
                         std::pair<uint32_t, uint32_t> bpRange,
                         const grgl::NodeIDList& onlySamples);

void emitZygosityInfo(grgl::GRGPtr& grg,
                      std::ostream& outStream,
                      std::pair<uint32_t, uint32_t> bpRange,
                      const grgl::NodeIDList& onlySamples);

void emitBeta(const grgl::GRGPtr& grg, 
              const std::string& phenotype, 
              std::ostream& outStream,
              bool betaOnly = false);

#endif /* GRGL_CALCULATIONS_H */