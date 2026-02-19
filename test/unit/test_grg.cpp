#include <gtest/gtest.h>

#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"
#include "grgl/transform.h"

#include "common_grgs.h"
#include "common_visitors.h"
#include "grgl/visitor.h"
#include "grg_helpers.h"

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
        ASSERT_EQ(g1.getNodesAndMutations().size(), 1);

        g1.addMutation(m2, INVALID_NODE_ID);
        ASSERT_EQ(g1.getUnmappedMutations().size(), 1);
        ASSERT_EQ(g1.getMutationsForNode(0).size(), 0);
        ASSERT_EQ(g1.getMutationsForNode(1).size(), 1);
        ASSERT_EQ(g1.getNodesAndMutations().size(), 2);
    }
}

TEST(GRG, DotProductGood) {
    GRGPtr grg = depth3BinTree();
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
    GRGPtr grg = depth3BinTree();
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

TEST(GRGHelper, CoalsForParent) {
    GRGPtr grg = sample8Grg();
    UnmanagedNodeData<NodeIDList> nodeToIndivs;
    nodeToIndivs.add(grg, 10, NodeIDListUPtr(new NodeIDList({0, 2})));
    nodeToIndivs.add(grg, 11, NodeIDListUPtr(new NodeIDList({0, 3})));
    nodeToIndivs.add(grg, 12, NodeIDListUPtr(new NodeIDList({3})));
    NodeIDList children = {10,11,12};
    NodeIDList uncoalesced;
    const NodeIDSizeT coals = getCoalsForParent(grg, nodeToIndivs, 14, children, uncoalesced);
    ASSERT_EQ(coals, 2);
    NodeIDList expected = {2};
    ASSERT_EQ(uncoalesced, expected);

    // Coalesce the last one
    children = {13, 9};
    nodeToIndivs.add(grg, 13, NodeIDListUPtr(new NodeIDList({2})));
    nodeToIndivs.add(grg, 9, NodeIDListUPtr(new NodeIDList({1})));
    const NodeIDSizeT coals2 = getCoalsForParent(grg, nodeToIndivs, 14, children, uncoalesced);
    ASSERT_EQ(coals2, 1);
    ASSERT_EQ(uncoalesced.size(), 1);

    // Start over, and test the no coalescences case
    nodeToIndivs.m_data.clear();
    nodeToIndivs.add(grg, 10, NodeIDListUPtr(new NodeIDList({0, 1})));
    nodeToIndivs.add(grg, 11, NodeIDListUPtr(new NodeIDList({2})));
    nodeToIndivs.add(grg, 12, NodeIDListUPtr(new NodeIDList({3})));
    children = {10,11,12};
    uncoalesced.clear();
    const NodeIDSizeT coals3 = getCoalsForParent(grg, nodeToIndivs, 14, children, uncoalesced);
    ASSERT_EQ(coals3, 0);
    expected = {0, 1, 2, 3};
    ASSERT_EQ(uncoalesced, expected);
}

TEST(GRG, MatmulXTX) {
    GRGPtr grg = depth3BinTree();
    ASSERT_EQ(grg->numSamples(), 4);
    ASSERT_EQ(grg->numNodes(), 7);
    // Lets add some mutations just so we can test matrix multiplication
    grg->addMutation(Mutation(1, "A", "G"), 6);
    grg->addMutation(Mutation(2, "T", "G"), 5);
    grg->addMutation(Mutation(3, "A", "C"), 4);
    ASSERT_EQ(grg->numMutations(), 3);
    ASSERT_EQ(grg->getPloidy(), 2);

    // TODO: we need a way to drop all the individual coalescences, so that this test
    // is more legit -- also just in general, not all GRGs need coalescences.
    //ASSERT_FALSE(grg->hasIndividualCoals());

    std::vector<int64_t> inputVect = {1, 1, 1, 1};
    std::vector<int64_t> outputVect(3);

    // When we don't have coalescence counts, this should fail!
    auto doMatmul = [&]() {
        grg->matrixMultiplication<int64_t, int64_t, false>(
        inputVect.data(), 4, 1, TraversalDirection::DIRECTION_UP, outputVect.data(), 3,
        false, false, nullptr, grgl::GRG::NIE_XTX, nullptr);
    };
    ASSERT_THROW(doMatmul(), grgl::ApiMisuseFailure);

    // Now calculate the coalescences and run again.
    calculateMissingCoals(grg);
    outputVect[0] = 0;
    outputVect[1] = 0;
    outputVect[2] = 0;
    doMatmul();
    std::vector<int64_t> expected = {8, 4, 4};
    ASSERT_EQ(expected, outputVect);
}


TEST(GRG, NegativeNodes) {
    GRGPtr grg = depth3BinTree();
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    ASSERT_TRUE(grg->nodesAreOrdered());

    MutableGRGPtr mutGRG = std::dynamic_pointer_cast<MutableGRG>(grg);
    // We can maintain the topological order if we:
    // 1: Add a new node and connect it as a parent
    const NodeID node1 = mutGRG->makeNode();
    mutGRG->connect(node1, grg->numNodes() - 2);
    ASSERT_TRUE(grg->nodesAreOrdered());
    // 2: Add a new "negative" node and connect it as a child
    // Nodes are no longer ordered by ID, but they are still trivially topologically ordered
    const NodeID node2 = mutGRG->makeNode(1, /*negative=*/true);
    mutGRG->connect(0, -node2);
    ASSERT_FALSE(grg->nodesAreOrdered());
    ASSERT_TRUE(grg->nodesAreTopo());

    // But we will break it if we add a negative node as a parent!
    mutGRG->connect(-node2, 3);
    ASSERT_FALSE(grg->nodesAreTopo());

    // Just check some helper methods
    const NodeID smallNode = 1;
    const NodeID bigNode = MAX_GRG_NODES - 1;
    const NodeID smallNeg = -smallNode;
    const NodeID bigNeg = -bigNode;
    // Type casting should work
    ASSERT_TRUE((SignedNodeID)bigNeg < (SignedNodeID)smallNeg);
    ASSERT_TRUE(nodeIsNegative(bigNeg));
    ASSERT_TRUE(nodeIsNegative(smallNeg));
    ASSERT_FALSE(nodeIsNegative(bigNode));
    ASSERT_FALSE(nodeIsNegative(smallNode));
    ASSERT_EQ(nodeStripNegative(bigNeg), bigNode);
    ASSERT_EQ(nodeStripNegative(smallNeg), smallNode);
}

TEST(GRG, ChangeSamples) {
    MutableGRGPtr grg = std::dynamic_pointer_cast<MutableGRG>(depth3BinTree());
    ASSERT_EQ(grg->numSamples(), 4);
    ASSERT_EQ(grg->numNodes(), 7);
    // Lets add some mutations just so we can test matrix multiplication
    grg->addMutation(Mutation(1, "A", "G"), 6);
    grg->addMutation(Mutation(2, "T", "G"), 5);
    grg->addMutation(Mutation(3, "A", "C"), 4);
    ASSERT_EQ(grg->numMutations(), 3);

    // Create 4 new nodes, beneath the old samples, and turn them into samples.
    NodeID newest = grg->makeNode(4, /*negative=*/true);
    ASSERT_EQ(grg->numNodes(), 11);
    grg->connect(0, -newest);
    newest++;
    grg->connect(1, -newest);
    newest++;
    grg->connect(2, -newest);
    newest++;
    grg->connect(3, -newest);
    ASSERT_LT(newest, grg->numNodes());
    // Nodes are no longer ordered according to NodeID, but they are still trivially
    // topologically ordered.
    ASSERT_FALSE(grg->nodesAreOrdered());
    ASSERT_TRUE(grg->nodesAreTopo());

    // Get the old samples
    NodeIDList expected = {0, 1, 2, 3};
    ASSERT_EQ(expected, grg->getSampleNodes());

    // Test a matmul with the old samples
    std::vector<int64_t> inputVect = {11, 41, 37, 73};
    std::vector<int64_t> outputVect(3);
    grg->matrixMultiplication<int64_t, int64_t, false>(inputVect.data(), 4, 1, TraversalDirection::DIRECTION_UP, outputVect.data(), 3);
    const std::vector<int64_t> expectedMat = {162, 110, 52};
    ASSERT_EQ(expectedMat, outputVect);

    // Make new samples
    NodeIDList newSamples = {7, 8, 9, 10};
    ASSERT_EQ(newest, newSamples.back());
    grg->setSamples(newSamples);
    ASSERT_EQ(newSamples, grg->getSampleNodes());
    ASSERT_TRUE(grg->isSample(7));
    ASSERT_TRUE(grg->isSample(10));
    ASSERT_FALSE(grg->isSample(0));
    ASSERT_FALSE(grg->isSample(3));

    // Get the topological order
    NodeIDList topo = grg->getOrderedNodes(TraversalDirection::DIRECTION_UP);
    std::vector<bool> saw(grg->numNodes(), false);
    ASSERT_EQ(topo.size(), grg->numNodes());
    for (NodeID node : topo) {
        ASSERT_FALSE(saw.at(node));
        switch (node) {
            case 0:
                ASSERT_TRUE(saw.at(7));
                break;
            case 1:
                ASSERT_TRUE(saw.at(8));
                break;
            case 2:
                ASSERT_TRUE(saw.at(9));
                break;
            case 3:
                ASSERT_TRUE(saw.at(10));
                break;
            default:
                break;
        }
        saw.at(node) = true;
    }

    // Test a matmul with the new samples!
    outputVect[0] = 0;
    outputVect[1] = 0;
    outputVect[2] = 0;
    grg->matrixMultiplication<int64_t, int64_t, false>(inputVect.data(), 4, 1, TraversalDirection::DIRECTION_UP, outputVect.data(), 3);
    ASSERT_EQ(expectedMat, outputVect);

    // Now truly break the topological order and ensure that matmul STILL works (it will use
    // a visitor)
    grg->connect(grg->numNodes() - 1, grg->makeNode());
    ASSERT_FALSE(grg->nodesAreTopo());
    outputVect[0] = 0;
    outputVect[1] = 0;
    outputVect[2] = 0;
    grg->matrixMultiplication<int64_t, int64_t, false>(inputVect.data(), 4, 1, TraversalDirection::DIRECTION_UP, outputVect.data(), 3);
    ASSERT_EQ(expectedMat, outputVect);

    // Now write this GRG to disk, _WITHOUT_ simplification. Reload it, and we should still get
    // the same matmul result
    std::string testFile = "test.change_samples.grg";
    saveGRG(grg, testFile, /*allowSimplify=*/false);
    GRGPtr grg2 = loadImmutableGRG(testFile);
    ASSERT_TRUE(grg2->nodesAreOrdered()); // Serializing re-ordered everything.
    ASSERT_EQ(grg2->numSamples(), 4);
    ASSERT_EQ(grg2->numNodes(), 11);
    outputVect[0] = 0;
    outputVect[1] = 0;
    outputVect[2] = 0;
    grg2->matrixMultiplication<int64_t, int64_t, false>(inputVect.data(), 4, 1, TraversalDirection::DIRECTION_UP, outputVect.data(), 3);
    ASSERT_EQ(expectedMat, outputVect);

    // Samples should have been renumbered!
    expected = {0, 1, 2, 3};
    ASSERT_EQ(expected, grg2->getSampleNodes());
}

} // namespace grgl
