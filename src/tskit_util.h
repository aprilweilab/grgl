#ifndef GRGL_TSKIT_UTIL_H
#define GRGL_TSKIT_UTIL_H

#include <tskit.h>

#define TSKIT_OK_OR_EXIT(ok, msg) do { \
    int tskit_ok_val = (ok); \
    if (tskit_ok_val != 0) { \
        std::cerr << msg << ": " << tsk_strerror(tskit_ok_val) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
} while(0)

namespace grgl {

/**
 * The given node is not connected to the given marginal tree.
 */
inline bool tsIsDisconnected(const tsk_tree_t *tree, const tsk_id_t treeNode) {
    return (TSK_NULL == tree->parent[treeNode]) && (0 == tree->num_children[treeNode]);
}

}

#endif /* GRGL_TSKIT_UTIL_H */