#ifndef GRGL_TRANSFORM_H
#define GRGL_TRANSFORM_H

#include "grg.h"

namespace grgl {

/**
 * Make a single pass over the mutable GRG and reduce the size by building hierarchy where there
 * should be some.
 *
 * @param[in] mutGRG The MutableGRG that will be modified.
 * @return The number of edges that were removed.
 */
size_t reduceGRG(const MutableGRGPtr& mutGRG);

/**
 * Reduce a GRG until one of the conditions is met.
 *
 * @param[in] mutGRG The MutableGRG that will be modified.
 * @param[in] iterations Maximum number of graph iterations.
 * @param[in] minDropped Minimum number of edges dropped by a single iteration.
 * @param[in] fractionDropped Total fraction of edges that have been dropped.
 * @param[in] verbose Set to true to get stdout information per iteration. Default: false.
 * @return The number of iterations that were performed.
 */
size_t reduceGRGUntil(
    const MutableGRGPtr& mutGRG, size_t iterations, size_t minDropped, double fractionDropped, bool verbose = false);

/**
 * Calculate any missing coalescence information in a GRG. You can use this to compute coalescences
 * for a large graph, though it is likely to scale poorly (RAM-wise). Generally, this is used
 * to "fix up" the coalescences after modifying a GRG -- any modifications that break coalescence
 * should be marked with COAL_COUNT_NOT_SET and this will fix them.
 *
 * @param[in] grg The graph to traverse.
 * @return The number of nodes that had their coalescence counts updated.
 */
size_t calculateMissingCoals(const GRGPtr& grg);

} // namespace grgl

#endif /* GRGL_TRANSFORM_H */