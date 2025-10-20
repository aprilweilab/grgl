// #define GRGL_CUDA_ENABLED
#ifdef GRGL_CUDA_ENABLED
#include <gtest/gtest.h>

#include "grgl/cuda_matmul.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"

#include "common_grgs.h"
#include "common_visitors.h"
#include "grgl/visitor.h"

#include <fstream>

TEST(CudaGRG, CudaCSRConstruction) {
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

    
    // Now we serialize and deserialize the graph, resulting in an _ordered_ GRG which means
    // the methods invoked under the hood will be slightly different.
    const char * const testFile = "test.grg.dotproductgood.grg";
    std::ofstream outStream(testFile);
    writeGrg(grg, outStream, false);
    outStream.close();
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    GPUCsrVisitor<uint32_t, double> visitor(grg->numNodes(), grg->numEdges());
    auto seeds = grg->getRootNodes();
    grg->visitDfs(
        visitor,
        TraversalDirection::DIRECTION_DOWN,
        seeds
    );
    CudaCSR<uint32_t, double> gpu_csr;
    visitor.constructCSR(grg.get(), gpu_csr);

    ASSERT_EQ(gpu_csr.nnz, grg->numEdges());
    ASSERT_EQ(gpu_csr.num_rows, grg->numNodes());
    ASSERT_GE(gpu_csr.max_height, 3);
    // ASSERT_EQ(gpu_csr.max_height, 3);
}


TEST(CudaGRG, MatMult) {
    GRGPtr grg = depth3BinTree(false);
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    ASSERT_FALSE(grg->nodesAreOrdered());

    Mutation m1(5, "A", "G");
    grg->addMutation(m1, 1); // Attach directly to sample 1
    Mutation m2(6, "A", "G");
    grg->addMutation(m2, 5); // Affects the two samples beneath 5
    Mutation m3(7, "A", "G");
    grg->addMutation(m3, 6); // Affects all four samples beneath 6

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
    std::vector<double> mutValues(3, 1.0); // Input vector [1, 2, 1]
    mutValues[1] = 2.0;
    auto result3 = grg->matMulGPU(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result3.size(), grg->numSamples());
    std::vector<double> expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result3, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    std::vector<double> sampleValues(4, 1.0); 
    auto result4 = grg->matMulGPU(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result4.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result4, expect);
}



TEST(CudaGRG, MatMultFromEnvVar) {
    // get testFile from env var GRGL_TEST_GRG_FILE if exists
    const char* env_p = std::getenv("GRGL_TEST_GRG_FILE");
    if (env_p == nullptr) {
        std::cout << "Environment variable GRGL_TEST_GRG_FILE not set, skipping" << std::endl;
        return;
    }
    const char * const testFile = (env_p != nullptr) ? env_p : "test.grg.dotproductgood.grg";
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    auto grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    // Top-down dot-product
    std::vector<double> mutValues(grg->numMutations(), 1.0);
    mutValues[0] = 2.0; // make sure input vector is not all ones
    mutValues[1] = 3.0;
    mutValues[grg->numMutations() - 1] = 4.0;
    auto result0 = grg->matMulGPU(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    auto result1 = grg->matMul(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);


    // Bottom-up dot-product 
    std::vector<double> sampleValues(grg->numSamples(), 1.0);
    sampleValues[0] = 2.0; // make sure input vector is not all ones
    sampleValues[1] = 3.0;
    sampleValues[grg->numSamples() - 1] = 4.0;
    auto result2 = grg->matMulGPU(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    auto result3 = grg->matMul(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "matMul function successfully executed on GRG from file: " << testFile << std::endl;
}

TEST(CudaGRG, MatrixMultMultiRow) {
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

    
    // Now we serialize and deserialize the graph, resulting in an _ordered_ GRG which means
    // the methods invoked under the hood will be slightly different.
    const char * const testFile = "test.grg.dotproductgood.grg";
    std::ofstream outStream(testFile);
    writeGrg(grg, outStream, false);
    outStream.close();
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    // Top-down multiplication
    std::vector<double> mutValues(6, 1.0);
    mutValues.at(1) = 2.0;
    mutValues.at(3) = 0.0;
    auto result1 = grg->matMulGPU(mutValues, 2, TraversalDirection::DIRECTION_DOWN);
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

TEST(CudaGRG, MatMultMultiRowFromEnvVar) {
    // get Row count from env var GRGL_TEST_MATMUL_ROWS if exists
    char* env_p = std::getenv("GRGL_TEST_MATMUL_ROWS");
    size_t rowCount = 1;
    if (env_p != nullptr) {
        rowCount = std::stoul(env_p);
    } else {
        return;
    }
    std::cout << "Using row count " << rowCount << " for MatMultMultiRowFromEnvVar test" << std::endl;

    // get testFile from env var GRGL_TEST_GRG_FILE if exists
    env_p = std::getenv("GRGL_TEST_GRG_FILE");
    if (env_p == nullptr) {
        std::cout << "Environment variable GRGL_TEST_GRG_FILE not set, skipping" << std::endl;
        return;
    }
    const char * const testFile = (env_p != nullptr) ? env_p : "test.grg.dotproductgood.grg";
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    auto grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    // Top-down dot-product
    std::vector<double> mutValues(grg->numMutations() * rowCount, 1.0);
    mutValues[0] = 2.0; // make sure input vector is not all ones
    mutValues[1] = 3.0;
    mutValues[grg->numMutations() * rowCount - 1] = 4.0;
    auto result0 = grg->matMulGPU(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN);
    auto result1 = grg->matMul(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    std::vector<double> sampleValues(grg->numSamples() * rowCount, 1.0);
    sampleValues[0] = 2.0; // make sure input vector is not all ones
    sampleValues[1] = 3.0;
    sampleValues[grg->numSamples() * rowCount - 1] = 4.0;
    auto result2 = grg->matMulGPU(sampleValues, rowCount, TraversalDirection::DIRECTION_UP);
    auto result3 = grg->matMul(sampleValues, rowCount, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "matMul function successfully executed on GRG from file: " << testFile << std::endl;
}

#endif // GRGL_CUDA_ENABLED