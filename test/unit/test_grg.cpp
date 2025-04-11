#include <gtest/gtest.h>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"

#include "common_grgs.h"
#include "common_visitors.h"
#include "grgl/visitor.h"

#include <fstream>

namespace grgl {

TEST(GRG, Mutations) {
    Mutation m1(1, "A", "G");
    ASSERT_TRUE(!m1.isEmpty());
    Mutation m2(5, "T", "G");
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
    GRGPtr grg = depth3BinTree(/*keepNodeOrder=*/false);
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
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
    auto result1 = grg->matMul(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result1.size(), grg->numSamples());
    std::vector<double> expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result1, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    std::vector<double> sampleValues(4, 1.0); // Input vector [1, 1, 1, 1]
    auto result2 = grg->matMul(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result2, expect);

    // Now we serialize and deserialize the graph, resulting in an _ordered_ GRG which means
    // the methods invoked under the hood will be slightly different.
    const char * const testFile = "test.grg.dotproductgood.grg";
    std::ofstream outStream(testFile);
    writeGrg(grg, outStream, false);
    outStream.close();
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    // Top-down dot-product
    auto result3 = grg->matMul(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result3.size(), grg->numSamples());
    expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result3, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    auto result4 = grg->matMul(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result4.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result4, expect);
}

TEST(GRG, MatrixMultGood) {
    GRGPtr grg = depth3BinTree(/*keepNodeOrder=*/false);
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    ASSERT_FALSE(grg->nodesAreOrdered());

    Mutation m1(5, "A", "G");
    grg->addMutation(m1, 1); // Attach directly to sample 1
    Mutation m2(6, "A", "G");
    grg->addMutation(m2, 5); // Affects the two samples beneath 5
    Mutation m3(7, "A", "G");
    grg->addMutation(m3, 6); // Affects all four samples beneath 6

    // Top-down multiplication
    std::vector<double> mutValues(6, 1.0);
    mutValues.at(1) = 2.0;
    mutValues.at(3) = 0.0;
    auto result1 = grg->matMul(mutValues, 2, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result1.size(), 2*grg->numSamples());
    std::vector<double> expect = {1.0, 2.0, 3.0, 3.0, 1.0, 1.0, 2.0, 2.0};
    ASSERT_EQ(result1, expect);

    // Bottom-up multiplication. This is just allele freq counts, and then the allele
    // frequencies (second row)
    std::vector<double> sampleValues(8, 1.0);
    for (size_t i = 4; i < sampleValues.size(); i++) {
        sampleValues[i] /= (double)grg->numSamples();
    }
    auto result2 = grg->matMul(sampleValues, 2, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), 2*grg->numMutations());
    expect = {1.0, 2.0, 4.0, 0.25, 0.5, 1.0};
    for (size_t i = 0; i < result2.size(); i++) {
        ASSERT_NEAR(result2.at(i), expect.at(i), 0.001);
    }
}

TEST(GRG, MatMulBad) {
    GRGPtr grg = depth3BinTree();
    Mutation m1(5, "A", "G");
    grg->addMutation(m1, 1); // Attach directly to sample 1
    Mutation m2(6, "A", "G");
    grg->addMutation(m2, 5); // Affects the two samples beneath 5
    Mutation m3(7, "A", "G");

	// Wrong input size
    std::vector<double> mutValues(4);
    ASSERT_THROW(
        grg->matMul(mutValues, 1, TraversalDirection::DIRECTION_DOWN),
        ApiMisuseFailure
    );

	// Wrong input size
    std::vector<double> sampleValues(3);
    ASSERT_THROW(
        grg->matMul(sampleValues, 1, TraversalDirection::DIRECTION_UP),
        ApiMisuseFailure
    );

    // Input size is not divisible by # of rows
    std::vector<double> sampleValues2(9);
    ASSERT_THROW(
        grg->matMul(sampleValues2, 2, TraversalDirection::DIRECTION_UP),
        ApiMisuseFailure
    );

	// Wrong input size
    std::vector<double> sampleValues3(9);
    ASSERT_THROW(
        grg->matMul(sampleValues3, 9, TraversalDirection::DIRECTION_UP),
        ApiMisuseFailure
    );
}

TEST(GRG, NoUpEdges) {
    GRGPtr grg = depth3BinTree();
    const NodeIDList rootsBefore = grg->getRootNodes();

	// Serialize and deserialize the GRG
    const char * const testFile = "test.no_up_edges.grg";
    std::ofstream outStream(testFile);
    writeGrg(grg, outStream, /*allowSimplify=*/false);
    outStream.close();
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    grg = readImmutableGrg(inStream, /*loadUpEdges=*/false);

    const NodeIDList rootsAfter = grg->getRootNodes();
    ASSERT_EQ(rootsBefore, rootsAfter);
}

TEST(GRG, TestTopoVisit) {
    GRGPtr grg = depth3BinTree();
    class TestVisitor : public GRGVisitor {
    public:
        bool visit(const grgl::GRGPtr& grg,
                   const grgl::NodeID nodeId,
                   const grgl::TraversalDirection direction,
                   const grgl::DfsPass dfsPass) override {
            m_visited++;
            if (nodeId >= 4) {
                return false;
            }
            return true;
        }
        NodeIDSizeT m_visited{};
    };
    TestVisitor visitor;
    grg->visitTopo(visitor, TraversalDirection::DIRECTION_UP, grg->getSampleNodes());
    ASSERT_EQ(visitor.m_visited, 6);

    TestVisitor visitorDense;
    grg->visitTopoNodeOrderedDense(visitorDense, TraversalDirection::DIRECTION_UP, grg->getSampleNodes());
    ASSERT_EQ(visitorDense.m_visited, 6);

    TestVisitor visitorSparse;
    grg->visitTopoNodeOrderedSparse(visitorSparse, TraversalDirection::DIRECTION_UP, grg->getSampleNodes());
    ASSERT_EQ(visitorSparse.m_visited, 6);
}

TEST(GRG, TestFrontier) {
    GRGPtr grg = sample8Grg();

    const NodeIDList downSeeds = {13, 11};
    FrontierVisitor downVisitor(downSeeds);
    grg->visitTopoNodeOrderedDense(downVisitor, TraversalDirection::DIRECTION_DOWN, downSeeds);
    const NodeIDList expectedDown = {10, 8};
    ASSERT_EQ(downVisitor.m_frontier, expectedDown);

    const NodeIDList upSeeds = {1, 3};
    FrontierVisitor upVisitor(upSeeds);
    grg->visitTopoNodeOrderedDense(upVisitor, TraversalDirection::DIRECTION_UP, upSeeds);
    const NodeIDList expectedUp = {11, 13};
    ASSERT_EQ(upVisitor.m_frontier, expectedUp);
}

}