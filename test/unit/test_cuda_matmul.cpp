// #define GRGL_CUDA_ENABLED
#ifdef GRGL_CUDA_ENABLED

#include <gtest/gtest.h>

#include "grgl/cuda/gpu_grg.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"

#include "common_grgs.h"
#include "common_visitors.h"
#include "grgl/visitor.h"

#include <fstream>

TEST(GPUGRG, Construction) {
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
    CHECK_CUDA_LAST_ERROR();
    GPUGRG gpu_grg = convertGRGToGPUGRG(grg);

    ASSERT_EQ(gpu_grg.numNodes(), grg->numNodes());
    ASSERT_EQ(gpu_grg.numEdges(), grg->numEdges());
    ASSERT_EQ(gpu_grg.numSamples(), grg->numSamples());
    ASSERT_EQ(gpu_grg.numMutations(), grg->numMutations());
    ASSERT_GE(gpu_grg.maxHeight(), 3);

    CHECK_CUDA_LAST_ERROR();


}

TEST(GPUGRG, StoreAndLoad) {
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

    GPUGRG gpu_grg = convertGRGToGPUGRG(grg);
    CHECK_CUDA_LAST_ERROR();
    ASSERT_EQ(gpu_grg.numNodes(), grg->numNodes());
    ASSERT_EQ(gpu_grg.numEdges(), grg->numEdges());
    ASSERT_EQ(gpu_grg.numSamples(), grg->numSamples());
    ASSERT_EQ(gpu_grg.numMutations(), grg->numMutations());
    ASSERT_GE(gpu_grg.maxHeight(), 3);

    const char * const gpuGrgFile = "test.gpu_grg.storeload.gpugrg";
    storeGPUGRGToDisk(gpu_grg, gpuGrgFile);
    GPUGRG loaded_gpu_grg = loadGPUGRGFromDisk(gpuGrgFile);

    CHECK_CUDA_LAST_ERROR();
    ASSERT_EQ(loaded_gpu_grg.numNodes(), gpu_grg.numNodes());
    ASSERT_EQ(loaded_gpu_grg.numEdges(), gpu_grg.numEdges());
    ASSERT_EQ(loaded_gpu_grg.numSamples(), gpu_grg.numSamples());
    ASSERT_EQ(loaded_gpu_grg.numMutations(), gpu_grg.numMutations());
    ASSERT_EQ(loaded_gpu_grg.maxHeight(), gpu_grg.maxHeight());

    CHECK_CUDA_LAST_ERROR();
}

// This test would do mat mul
// downwards and upwards, single element and multiple element
TEST(GPUGRG, MatMult) {
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

    GPUGRG gpu_grg = convertGRGToGPUGRG(grg);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Convert CUDA kernel launch error: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("CUDA kernel launch failed");
    }

    ASSERT_EQ(gpu_grg.numNodes(), grg->numNodes());
    ASSERT_EQ(gpu_grg.numEdges(), grg->numEdges());
    ASSERT_EQ(gpu_grg.numSamples(), grg->numSamples());
    ASSERT_EQ(gpu_grg.numMutations(), grg->numMutations());
    ASSERT_GE(gpu_grg.maxHeight(), 3);

    std::vector<double> mutValues(3, 1.0); // Input vector [1, 2, 1]
    mutValues[1] = 2.0;
    auto result1 = gpu_grg.matMulBlocking(mutValues, 1, TraversalDirection::DIRECTION_DOWN, false);
    ASSERT_EQ(result1.size(), grg->numSamples());
    std::vector<double> expect = {1.0, 2.0, 3.0, 3.0};
    ASSERT_EQ(result1, expect);

    // Bottom-up dot-product (this is just allele freq counts)
    std::vector<double> sampleValues(4, 1.0); 
    auto result2 = gpu_grg.matMulBlocking(sampleValues, 1, TraversalDirection::DIRECTION_UP, false);
    ASSERT_EQ(result2.size(), grg->numMutations());
    expect = {1.0, 2.0, 4.0};
    ASSERT_EQ(result2, expect);

    // Top-down multiplication
    std::vector<double> mutValues2(6, 1.0);
    mutValues2.at(1) = 2.0;
    mutValues2.at(3) = 0.0;
    auto result3 = gpu_grg.matMulBlocking(mutValues2, 2, TraversalDirection::DIRECTION_DOWN, false);
    ASSERT_EQ(result3.size(), 2*grg->numSamples());
    expect = {1.0, 2.0, 3.0, 3.0, 1.0, 1.0, 2.0, 2.0};
    ASSERT_EQ(result3, expect);

    
    // Bottom-up multiplication. This is just allele freq counts, and then the allele
    // frequencies (second row)
    std::vector<double> sampleValues2(8, 1.0);
    for (size_t i = 4; i < sampleValues2.size(); i++) {
        sampleValues2[i] /= (double)grg->numSamples();
    }
    auto result4 = gpu_grg.matMulBlocking(sampleValues2, 2, TraversalDirection::DIRECTION_UP, false);
    ASSERT_EQ(result4.size(), 2*grg->numMutations());
    expect = {1.0, 2.0, 4.0, 0.25, 0.5, 1.0};
    for (size_t i = 0; i < result4.size(); i++) {
        ASSERT_NEAR(result4.at(i), expect.at(i), 0.001);
    }
}

TEST(GPUGRG, MatMulEnvVar) {
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

    GPUGRG gpu_grg = convertGRGToGPUGRG(grg);
    std::cout << "Loaded and converted GRG from file: " << testFile << std::endl;

    // Top-down dot-product
    std::vector<double> mutValues(grg->numMutations(), 1.0);
    mutValues[0] = 2.0; // make sure input vector is not all ones
    mutValues[1] = 3.0;
    mutValues[grg->numMutations() - 1] = 4.0;
    auto result0 = gpu_grg.matMulBlocking(mutValues, 1, TraversalDirection::DIRECTION_DOWN, false);
    auto result1 = grg->matMul(mutValues, 1, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    std::vector<double> sampleValues(grg->numSamples(), 1.0);
    sampleValues[0] = 2.0; // make sure input vector is not all ones
    sampleValues[1] = 3.0;
    sampleValues[grg->numSamples() - 1] = 4.0;
    auto result2 = gpu_grg.matMulBlocking(sampleValues, 1, TraversalDirection::DIRECTION_UP, false);
    auto result3 = grg->matMul(sampleValues, 1, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    
    std::cout << "matMul function successfully executed on GRG from file: " << testFile << std::endl;

    // EXTRA XTX tests
    // Top-down dot-product
    result0 = gpu_grg.matMulBlocking(mutValues, 1, TraversalDirection::DIRECTION_DOWN, false, GPUGRG::NIE_XTX);
    result1.clear();
    result1.resize(grg->numSamples(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        mutValues.data(),
        grg->numMutations(),
        1,
        TraversalDirection::DIRECTION_DOWN,
        result1.data(),
        grg->numSamples(),
        false,
        false,
        nullptr,
        GRG::NIE_XTX
    );
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    result2 = gpu_grg.matMulBlocking(sampleValues, 1, TraversalDirection::DIRECTION_UP, false, GPUGRG::NIE_XTX);
    result3.clear();
    result3.resize(grg->numMutations(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        sampleValues.data(),
        grg->numSamples(),
        1,
        TraversalDirection::DIRECTION_UP,
        result3.data(),
        grg->numMutations(),
        false,
        false,
        nullptr,
        GRG::NIE_XTX
    );
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "XTX initialization tests completed successfully" << std::endl;
}

TEST(GPUGRG, MatMulEnvVarMultiRow) {
    char* env_p = std::getenv("GRGL_TEST_MATMUL_ROWS");
    size_t rowCount = 1;
    if (env_p != nullptr) {
        rowCount = std::stoul(env_p);
    } else {
        return;
    }
    std::cout << "Using row count " << rowCount << " for MatMultMultiRowFromEnvVar test" << std::endl;


    env_p = std::getenv("GRGL_TEST_GRG_FILE");
    if (env_p == nullptr) {
        std::cout << "Environment variable GRGL_TEST_GRG_FILE not set, skipping" << std::endl;
        return;
    }

    const char * const testFile = (env_p != nullptr) ? env_p : "test.grg.dotproductgood.grg";
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    auto grg = readImmutableGrg(inStream);
    ASSERT_TRUE(grg->nodesAreOrdered());

    GPUGRG gpu_grg_a = convertGRGToGPUGRG(grg);
    std::cout << "Loaded and converted GRG from file: " << testFile << std::endl;
    const char * const gpuGrgFile = "test.gpu_grg.storeload.gpugrg";
    storeGPUGRGToDisk(gpu_grg_a, gpuGrgFile);
    GPUGRG gpu_grg = loadGPUGRGFromDisk(gpuGrgFile);
    std::cout << "Loaded GPUGRG from disk: " << gpuGrgFile << std::endl;

    // Top-down dot-product
    std::vector<double> mutValues(grg->numMutations() * rowCount, 2.0);
    mutValues[0] = 2.0; // make sure input vector is not all ones
    mutValues[1] = 3.0;
    mutValues[grg->numMutations() - 1] = 4.0;
    mutValues[grg->numMutations() * rowCount - 1] = 5.0;
    auto result0 = gpu_grg.matMulBlocking(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, false);
    auto result1 = grg->matMul(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN);
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    std::vector<double> sampleValues(grg->numSamples() * rowCount, 4.0);
    sampleValues[0] = 2.0; // make sure input vector is not all ones
    sampleValues[1] = 3.0;
    sampleValues[grg->numSamples() - 1] = 4.0;
    sampleValues[grg->numSamples() * rowCount - 1] = 5.0;
    auto result2 = gpu_grg.matMulBlocking(sampleValues, rowCount, TraversalDirection::DIRECTION_UP, false);
    auto result3 = grg->matMul(sampleValues, rowCount, TraversalDirection::DIRECTION_UP);
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "Basic Tests Completed" << std::endl;

    /* Start Init Tests */
    size_t numCols;
    size_t outSize;

    // Init by vector
    // Top-down dot-product
    numCols = grg->numMutations();
    outSize = rowCount * grg->numSamples();
    std::vector<double> initVector(rowCount, 1.75);
    initVector[0] = 1.0;
    initVector[3] = 2.25;
    initVector[rowCount-2] = 10.0;
    result0 = gpu_grg.matMulBlocking(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, false, GPUGRG::NIE_VECTOR, initVector.data());
    result1.clear();
    result1.resize(rowCount * grg->numSamples(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        mutValues.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_DOWN,
        result1.data(),
        outSize,
        false,
        false,
        initVector.data(),
        GRG::NIE_VECTOR
    );
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    numCols = grg->numSamples();
    outSize = rowCount * grg->numMutations();
    result2 = gpu_grg.matMulBlocking(sampleValues, rowCount, TraversalDirection::DIRECTION_UP, false, GPUGRG::NIE_VECTOR, initVector.data());
    result3.clear();
    result3.resize(rowCount * grg->numMutations(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        sampleValues.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_UP,
        result3.data(),
        outSize,
        false,
        false,
        initVector.data(),
        GRG::NIE_VECTOR
    );
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "Init by Vector Test Completed." << std::endl;

    // Init by matrix
    // Top-down dot-product
    numCols = grg->numMutations();
    outSize = rowCount * grg->numSamples();
    std::vector<double> initMatrix(grg->numNodes() * rowCount, 2.5);
    initMatrix[0] = 5.0;
    initMatrix[grg->numNodes() - 1] = 6.0;
    initMatrix[grg->numNodes() * rowCount - 2] = 12.0;
    result0 = gpu_grg.matMulBlocking(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, false, GPUGRG::NIE_MATRIX, initMatrix.data());
    result1.clear();
    result1.resize(rowCount * grg->numSamples(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        mutValues.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_DOWN,
        result1.data(),
        outSize,
        false,
        false,
        initMatrix.data(),
        GRG::NIE_MATRIX
    );
    ASSERT_EQ(result0.size(), result1.size());
    ASSERT_EQ(result0, result1);

    // Bottom-up dot-product
    numCols = grg->numSamples();
    outSize = rowCount * grg->numMutations();
    result2 = gpu_grg.matMulBlocking(sampleValues, rowCount, TraversalDirection::DIRECTION_UP, false, GPUGRG::NIE_MATRIX, initMatrix.data());
    result3.clear();
    result3.resize(rowCount * grg->numMutations(), 0.0);
    grg->matrixMultiplication<double, double, false>(
        sampleValues.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_UP,
        result3.data(),
        outSize,
        false,
        false,
        initMatrix.data(),
        GRG::NIE_MATRIX
    );
    ASSERT_EQ(result2.size(), result3.size());
    ASSERT_EQ(result2, result3);
    std::cout << "Init by Matrix Test Completed." << std::endl;

    /* Start ByIndividual Tests */
    // Top-down dot-product with byIndividual
    numCols = grg->numMutations();
    outSize = rowCount * grg->numSamples() / grg->getPloidy();
    auto result4 = gpu_grg.matMulBlocking(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, true);
    std::vector<double> result5(outSize);
    grg->matrixMultiplication<double, double, false>(
        mutValues.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_DOWN,
        result5.data(),
        outSize,
        false,
        true
    );
    ASSERT_EQ(result4.size(), result5.size());
    ASSERT_EQ(result4, result5);

    numCols = grg->numSamples() / grg->getPloidy();
    outSize = rowCount * grg->numMutations();
    std::vector<double> sampleValuesInd(numCols * rowCount, 3.25);
    sampleValuesInd[0] = 2.0; // make sure input vector is not all ones
    sampleValuesInd[1] = 3.0;
    sampleValuesInd[numCols - 3] = 4.0;
    sampleValuesInd[numCols * rowCount - 2] = 5.0;

    std::vector<double> result6 = gpu_grg.matMulBlocking(sampleValuesInd, rowCount, TraversalDirection::DIRECTION_UP, true);
    std::vector<double> result7(outSize);
    grg->matrixMultiplication<double, double, false>(
        sampleValuesInd.data(),
        numCols,
        rowCount,
        TraversalDirection::DIRECTION_UP,
        result7.data(),
        outSize,
        false,
        true
    );
    ASSERT_EQ(result6.size(), result7.size());
    ASSERT_EQ(result6, result7);
    std::cout << "ByIndividual Test Completed." << std::endl;

}

TEST(GPUGRG, BenchMatMul) {
    char* env_p;
    bool direct_gpu = false;

    env_p = std::getenv("GRGL_BENCH_GPUGRG_FILE");
    if (env_p != nullptr) {
        direct_gpu = true;
    } 
    env_p = std::getenv("GRGL_BENCH_GRG_FILE");
    if (env_p == nullptr) {
        std::cout << "Environment variable GRGL_BENCH_GRG_FILE not set, skipping" << std::endl;
        return;
    }
    
    size_t iterations = 10;
    env_p = std::getenv("GRGL_BENCH_ITERATIONS");
    if (env_p != nullptr) {
        iterations = std::stoul(env_p);
    } else {
        return;
    }
    std::cout << "Using " << iterations << " iterations for BenchMatMul test" << std::endl;

    env_p = std::getenv("GRGL_BENCH_MATMUL_ROWS");
    size_t rowCount;
    if (env_p != nullptr) {
        rowCount = std::stoul(env_p);
    } else {
        rowCount = 1;
    }
    std::cout << "Using row count " << rowCount << " for MatMultMultiRowFromEnvVar test" << std::endl;

    if (direct_gpu) {
        env_p = std::getenv("GRGL_BENCH_GPUGRG_FILE");
        GPUGRG gpu_grg = loadGPUGRGFromDisk(env_p);
        env_p = std::getenv("GRGL_BENCH_GRG_FILE");
        const char * const testFile = (env_p != nullptr) ? env_p : "test.grg.dotproductgood.grg";
        grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
        auto grg = readImmutableGrg(inStream);
        ASSERT_TRUE(grg->nodesAreOrdered());

        std::cout << "Loaded GPUGRG from disk: " << env_p << std::endl;
        std::vector<double> mutValues(grg->numMutations() * rowCount, 2.0);
        mutValues[0] = 2.0; // make sure input vector is not all ones
        mutValues[1] = 3.0;
        mutValues[grg->numMutations() - 1] = 4.0;
        mutValues[grg->numMutations() * rowCount - 1] = 5.0;
        std::cout << "Starting benchmark DOWN matMul with " << iterations << " iterations." << std::endl;
        auto result0 = gpu_grg.matMulPerf(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, iterations);
        auto result1 = grg->matMul(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN);
        // ASSERT_EQ(result0.size(), result1.size());
        /*
        for (size_t i = 0; i < result0.size(); i++) {
            if (result0[i] != result1[i]) {
                std::cout << "Difference at index " << i << ": GPU result = " << result0[i]
                        << ", CPU result = " << result1[i] << std::endl;
            }
        }
            */
        // ASSERT_EQ(result0, result1);
        

        // Bottom-up dot-product

        std::vector<double> sampleValues(grg->numSamples() * rowCount, 4.0);
        sampleValues[0] = 2.0; // make sure input vector is not all ones
        sampleValues[1] = 3.0;
        sampleValues[grg->numSamples() - 1] = 4.0;
        sampleValues[grg->numSamples() * rowCount - 1] = 5.0;
        std::cout << "Starting benchmark UP matMul with " << iterations << " iterations." << std::endl;
        auto result2 = gpu_grg.matMulPerf(sampleValues, rowCount, TraversalDirection::DIRECTION_UP, iterations);
        auto result3 = grg->matMul(sampleValues, rowCount, TraversalDirection::DIRECTION_UP);
    } else {
        env_p = std::getenv("GRGL_BENCH_GRG_FILE");
        if (env_p == nullptr) {
            std::cout << "Environment variable GRGL_BENCH_GRG_FILE not set, skipping" << std::endl;
            return;
        }
        const char * const testFile = (env_p != nullptr) ? env_p : "test.grg.dotproductgood.grg";
        grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
        auto grg = readImmutableGrg(inStream);
        ASSERT_TRUE(grg->nodesAreOrdered());
        GPUGRG gpu_grg = convertGRGToGPUGRG(grg);

        // Top-down dot-product
        
        std::vector<double> mutValues(grg->numMutations() * rowCount, 2.0);
        mutValues[0] = 2.0; // make sure input vector is not all ones
        mutValues[1] = 3.0;
        mutValues[grg->numMutations() - 1] = 4.0;
        mutValues[grg->numMutations() * rowCount - 1] = 5.0;
        std::cout << "Starting benchmark DOWN matMul with " << iterations << " iterations." << std::endl;
        auto result0 = gpu_grg.matMulPerf(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN, iterations);
        auto result1 = grg->matMul(mutValues, rowCount, TraversalDirection::DIRECTION_DOWN);
        // ASSERT_EQ(result0.size(), result1.size());
        /*
        for (size_t i = 0; i < result0.size(); i++) {
            if (result0[i] != result1[i]) {
                std::cout << "Difference at index " << i << ": GPU result = " << result0[i]
                        << ", CPU result = " << result1[i] << std::endl;
            }
        }
            */
        // ASSERT_EQ(result0, result1);
        

        // Bottom-up dot-product

        std::vector<double> sampleValues(grg->numSamples() * rowCount, 4.0);
        sampleValues[0] = 2.0; // make sure input vector is not all ones
        sampleValues[1] = 3.0;
        sampleValues[grg->numSamples() - 1] = 4.0;
        sampleValues[grg->numSamples() * rowCount - 1] = 5.0;
        std::cout << "Starting benchmark UP matMul with " << iterations << " iterations." << std::endl;
        auto result2 = gpu_grg.matMulPerf(sampleValues, rowCount, TraversalDirection::DIRECTION_UP, iterations);
        auto result3 = grg->matMul(sampleValues, rowCount, TraversalDirection::DIRECTION_UP);
    }
    // ASSERT_EQ(result2.size(), result3.size());
    // manually compare all elements and print any differences
    /*
    for (size_t i = 0; i < result2.size(); i++) {
        if (result2[i] != result3[i]) {
            std::cout << "Difference at index " << i << ": GPU result = " << result2[i]
                      << ", CPU result = " << result3[i] << std::endl;
        }
    }
    */

    // ASSERT_EQ(result2, result3);
    
    // std::cout << "matMul function successfully executed on GRG from file: " << testFile << std::endl;
}

#endif // GRGL_CUDA_ENABLED