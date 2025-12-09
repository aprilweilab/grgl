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
size_t reduceGRG(const grgl::MutableGRGPtr& mutGRG);

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
size_t reduceGRGUntil(const grgl::MutableGRGPtr& mutGRG,
                      size_t iterations,
                      size_t minDropped,
                      double fractionDropped,
                      bool verbose = false);

} // namespace grgl

#endif /* GRGL_TRANSFORM_H */