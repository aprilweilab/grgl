#include <gtest/gtest.h>

#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"

using namespace grgl;

TEST(GRG, Mutations) {
    Mutation m1(1.1, "A", "G");
    ASSERT_TRUE(!m1.isEmpty());
    Mutation m2(5.1, "T", "G");
    ASSERT_TRUE(!m2.isEmpty());

    {
        MutableGRG g1(/*numSamples=*/2, /*ploidy=*/1);
        g1.addMutation(m1, 1);
        ASSERT_EQ(g1.getUnmappedMutations().size(), 0);
        ASSERT_EQ(g1.getMutationsForNode(0).size(), 0);
        ASSERT_EQ(g1.getMutationsForNode(1).size(), 1);
        ASSERT_EQ(g1.getNodeMutationPairs().size(), 1);

        g1.addMutation(m2, INVALID_NODE_ID);
        ASSERT_EQ(g1.getUnmappedMutations().size(), 1);
        ASSERT_EQ(g1.getMutationsForNode(0).size(), 0);
        ASSERT_EQ(g1.getMutationsForNode(1).size(), 1);
        ASSERT_EQ(g1.getNodeMutationPairs().size(), 2);
    }
}
