#include <gtest/gtest.h>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"

#include "common_grgs.h"
#include "grgl/visitor.h"

#include <fstream>

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

TEST(GRG, DotProductGood) {
    GRGPtr grg = depth3BinTree();
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    // We haven't serialized the GRG ever, so the nodes are unordered.
    ASSERT_FALSE(grg->nodesAreOrdered());

    Mutation m1(5, "A", "G");
    grg->addMutation(m1, 1); // Attach directly to sample 1
    Mutation m2(6, "A", "G");
    grg->addMutation(m2, 5); // Affects the two samples beneath 5
    Mutation m3(7, "A", "G");
    grg->addMutation(m3, 6); // Affects all four samples beneath 6

    // Top-down dot-product
    std::vector<double> mutValues(3, 1.0); // Input vector [1, 2, 1]
    mutValues.at(1) = 2.0;
    auto result1 = grg->dotProduct(mutValues, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result1.size(), grg->numSamples());
    std::vector<double> expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result1, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    std::vector<double> sampleValues(4, 1.0); // Input vector [1, 1, 1, 1]
    auto result2 = grg->dotProduct(sampleValues, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result2, expect);

    // Now we serialize and deserialize the graph, resulting in an _ordered_ GRG which means
    // the methods invoked under the hood will be slightly different.
    const char * const testFile = "test.grg.dotproductgood.grg";
    std::ofstream outStream(testFile);
    writeGrg(grg, outStream, false, false);
    outStream.close();
    std::ifstream inStream(testFile);
    grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    // Top-down dot-product
    auto result3 = grg->dotProduct(mutValues, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result3.size(), grg->numSamples());
    expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result3, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    auto result4 = grg->dotProduct(sampleValues, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result4.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result4, expect);
}

TEST(GRG, DotProductBad) {
    GRGPtr grg = depth3BinTree();
    Mutation m1(5, "A", "G");
    grg->addMutation(m1, 1); // Attach directly to sample 1
    Mutation m2(6, "A", "G");
    grg->addMutation(m2, 5); // Affects the two samples beneath 5
    Mutation m3(7, "A", "G");

	// Wrong input size
    std::vector<double> mutValues(4);
    ASSERT_THROW(
        grg->dotProduct(mutValues, TraversalDirection::DIRECTION_DOWN),
        ApiMisuseFailure
    );

	// Wrong input size
    std::vector<double> sampleValues(3);
    ASSERT_THROW(
        grg->dotProduct(sampleValues, TraversalDirection::DIRECTION_UP),
        ApiMisuseFailure
    );
}