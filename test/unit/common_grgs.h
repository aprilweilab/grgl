#ifndef COMMON_GRGS_H
#define COMMON_GRGS_H

#include "grgl/grg.h"
#include <memory>

using namespace grgl;

// Simple full binary tree of depth 3
//               |---pn6---|
//          |--pn4--|   |--pn5--|
//         s0      s1  s2      s3
inline MutableGRGPtr depth3BinTree() {
    const size_t nSamples = 4;
    MutableGRGPtr grg = std::make_shared<MutableGRG>(nSamples, 2);
    NodeID pn4 = grg->makeNode();
    NodeID pn5 = grg->makeNode();
    NodeID pn6 = grg->makeNode();
    grg->connect((NodeID)4, (NodeID)0);
    grg->connect((NodeID)4, (NodeID)1);
    grg->connect((NodeID)5, (NodeID)2);
    grg->connect((NodeID)5, (NodeID)3);
    grg->connect((NodeID)6, (NodeID)4);
    grg->connect((NodeID)6, (NodeID)5);
    return grg;
}

#endif /* COMMON_GRGS_H */