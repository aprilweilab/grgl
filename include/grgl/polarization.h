#ifndef GRG_POLARIZATION_H
#define GRG_POLARIZATION_H

#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"

#include <string>
#include <utility>
#include <vector>

namespace grgl {

struct PolarizationStats {
    size_t totalSeen{0};
    size_t emitted{0};
    size_t alreadyPolarized{0};
    size_t swapped{0};
    size_t droppedUnknown{0};
    size_t inconsistent{0};

    void reset() {
        totalSeen = 0;
        emitted = 0;
        alreadyPolarized = 0;
        swapped = 0;
        droppedUnknown = 0;
        inconsistent = 0;
    }
};

/**
 * Polarize a mutation and immediately remap it into the GRG with the provided carrier lists.
 */
bool polarizeMutation(const MutableGRGPtr& grg,
                      MutationId mutId,
                      const std::string& ancestralAllele,
                      PolarizationStats& stats);

/**
 * Polarize a batch of mutations, remapping any that flip in one mapMutations call.
 * Returns per-mutation success (true if emitted/kept).
 */
std::vector<bool> polarizeMutations(const MutableGRGPtr& grg,
                                    const std::vector<std::pair<MutationId, std::string>>& batch,
                                    PolarizationStats& stats);

} // namespace grgl

#endif // GRG_POLARIZATION_H
