#ifndef COMMON_GRGS_H
#define COMMON_GRGS_H

#include "grgl/grg.h"
#include <memory>

using namespace grgl;

// Simple full binary tree of depth 3
//               |---pn6---|
//          |--pn4--|   |--pn5--|
//         s0      s1  s2      s3
inline MutableGRGPtr depth3BinTree(const bool keepNodeOrder = true) {
    const size_t nSamples = 4;
    MutableGRGPtr grg = std::make_shared<MutableGRG>(nSamples, 2);
    NodeID pn4 = grg->makeNode(/*count=*/1, /*forceOrdered=*/keepNodeOrder);
    NodeID pn5 = grg->makeNode(/*count=*/1, /*forceOrdered=*/keepNodeOrder);
    NodeID pn6 = grg->makeNode(/*count=*/1, /*forceOrdered=*/keepNodeOrder);
    grg->connect((NodeID)4, (NodeID)0);
    grg->connect((NodeID)4, (NodeID)1);
    grg->connect((NodeID)5, (NodeID)2);
    grg->connect((NodeID)5, (NodeID)3);
    grg->connect((NodeID)6, (NodeID)4);
    grg->connect((NodeID)6, (NodeID)5);
    return grg;
}

// 8-sample GRG.
// pn8 -> (s0, s1)
// pn9 -> (pn8, s2)
// pn10 -> (s3, s4, s7)
// pn11 -> (pn8, pn10, s6)
// pn12 -> (pn10, s5)
// pn13 -> (pn12, pn9)
inline MutableGRGPtr sample8Grg(const bool keepNodeOrder = true) {
    const size_t nSamples = 8;
    MutableGRGPtr grg = std::make_shared<MutableGRG>(nSamples, 2);
    // We force this graph to maintain the topological flag, which means we have to ensure
    // we never create an edge from n1->n2 where n1<=n2.
    for (size_t i = 8; i <= 13; i++) {
        grg->makeNode(/*count=*/1, /*forceOrdered=*/keepNodeOrder);
    }
    grg->connect((NodeID)8, (NodeID)0);
    grg->connect((NodeID)8, (NodeID)1);
    grg->connect((NodeID)9, (NodeID)8);
    grg->connect((NodeID)9, (NodeID)2);
    grg->connect((NodeID)10, (NodeID)3);
    grg->connect((NodeID)10, (NodeID)4);
    grg->connect((NodeID)10, (NodeID)7);
    grg->connect((NodeID)11, (NodeID)8);
    grg->connect((NodeID)11, (NodeID)10);
    grg->connect((NodeID)11, (NodeID)6);
    grg->connect((NodeID)12, (NodeID)10);
    grg->connect((NodeID)12, (NodeID)5);
    grg->connect((NodeID)13, (NodeID)12);
    grg->connect((NodeID)13, (NodeID)9);
    return grg;
}

#endif /* COMMON_GRGS_H */