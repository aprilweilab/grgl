#include <gtest/gtest.h>

#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/file_vector.h"

#include <chrono>
#include <limits>
#include <random>
#include <fstream>

using namespace grgl;

using CSRBasic32 = CSRStorageImm<EagerFileVector, uint32_t, false, false>;
using CSRBasic64 = CSRStorageImm<EagerFileVector, uint64_t, false, false>;

TEST(CSRStorageImm, Success) {
    CSRBasic32 storage1(100);
    for (size_t i = 0; i < 100; i++) {
        storage1.setData(i, {(uint32_t)(i*2), (uint32_t)(i*3)});
    }
    storage1.finalizeNodes();

    std::vector<uint32_t> results;
    storage1.getData(55, results);
    std::vector<uint32_t> expect = {110, 165};
    ASSERT_EQ(results, expect);

    storage1.getData(0, results);
    expect = {110, 165, 0, 0};
    ASSERT_EQ(results, expect);

    ASSERT_EQ(storage1.getDatum(99, 1), 99*3);
    ASSERT_EQ(storage1.getDatum(11, 0), 22);

    CSRBasic32 reverse(50, true);
    for (size_t i = 50; i > 0; i--) {
        const size_t j = i-1;
        reverse.setData(j, {(uint32_t)(j+2), (uint32_t)(j+9)});
    }
    ASSERT_EQ(reverse.getDatum(6, 0), 8);
}

TEST(CSRStorageImm, Success64) {
    CSRBasic64 storage1(100);
    for (size_t i = 0; i < 99; i++) {
        storage1.setData(i, {(uint64_t)(i*2), (uint64_t)(i*3)});
    }
    storage1.setData(99,
        {std::numeric_limits<uint64_t>::max()-1, std::numeric_limits<uint64_t>::max()});
    storage1.finalizeNodes();

    std::vector<uint64_t> results;
    storage1.getData(55, results);
    std::vector<uint64_t> expect = {110, 165};
    ASSERT_EQ(results, expect);

    storage1.getData(0, results);
    expect = {110, 165, 0, 0};
    ASSERT_EQ(results, expect);

    ASSERT_EQ(storage1.getDatum(99, 1), std::numeric_limits<uint64_t>::max());
    ASSERT_EQ(storage1.getDatum(11, 0), 22);
}

TEST(CSRStorageImm, Fail) {
    CSRBasic32 storage1(100);
    storage1.setData(0, {99, 101, 103});
    // It is ok to skip indexes, but the order/sizes must be maintained. Here we added 0, then 3, which means that 1,2 are
    // set to be empty. So then when we try to set 1 to {9} it fails.
    storage1.setData(3, {0});
    ASSERT_THROW(
        storage1.setData(1, {9}),
        ApiMisuseFailure);
    // Wrong dataset size
    ASSERT_THROW(
        storage1.setData(0, {0, 1}),
        ApiMisuseFailure);

    CSRBasic32 storage2(3, true);
    storage2.setData(1, {99, 101, 103});
    // Wrong data; 1 already got set to empty because we set 0 first.
    ASSERT_THROW(
        storage2.setData(2, {99, 101, 103}),
        ApiMisuseFailure);
    storage2.setData(0, {100000});
}

#define PERFORM_TIMING 0

#if PERFORM_TIMING
TEST(CSRStorageImm, Speed) {
    const size_t seed = 42;
    std::mt19937 generator(seed);

    std::uniform_int_distribution<uint32_t> sampler(0, 2000000);

    const size_t numNodes = 1000000;

    CSRStorageImm<EagerFileVector, uin32_t, false, false> uncompressed(numNodes);
    CSRStorageImm<EagerFileVector, uin32_t, true, true> compressed(numNodes);
    size_t total = 0;
    for (size_t i = 0; i < numNodes; i++) {
        std::vector<uint32_t> edges(500);
        for (size_t j = 0; j < 500; j++) {
            edges[j] = sampler(generator);
        }
        std::sort(edges.begin(), edges.end());
        uncompressed.setData(i, edges);
        total += compressed.setData(i, edges);
    }
    std::cout << "Uncompressed: " << (numNodes * 500) * 4 << " bytes\n";
    std::cout << "Compressed: " << total << " bytes\n";

    auto operationStartTime = std::chrono::high_resolution_clock::now();
#define START_TIMING_OPERATION() operationStartTime = std::chrono::high_resolution_clock::now();
#define EMIT_TIMING_MESSAGE(msg)                                                                                       \
    do {                                                                                                               \
        std::cerr << msg                                                                                               \
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - \
                                                                           operationStartTime)                         \
                         .count()                                                                                      \
                  << " ms" << std::endl;                                                                               \
    } while (0)

    START_TIMING_OPERATION();
    std::vector<uint32_t> compEdges;
    for (size_t i = 0; i < numNodes; i++) {
        compEdges.clear();
        compressed.getData(i, compEdges);
        ASSERT_EQ(compEdges.size(), 500);
    }
    EMIT_TIMING_MESSAGE("Compressed traversal took ");

    START_TIMING_OPERATION();
    std::vector<uint32_t> uncEdges;
    for (size_t i = 0; i < numNodes; i++) {
        uncEdges.clear();
        uncompressed.getData(i, uncEdges);
        ASSERT_EQ(uncEdges.size(), 500);
    }
    EMIT_TIMING_MESSAGE("Uncompressed traversal took ");

    ASSERT_EQ(uncEdges, compEdges);
}

#endif

TEST(CSRStorageImm, Writing) {
    const size_t seed = 42;
    std::mt19937 generator(seed);

    std::uniform_int_distribution<uint32_t> sampler(0, 2000000);

    const size_t numNodes = 100;
    const size_t numEdges = 500;

    CSRStorageImm<EagerFileVector, uint32_t, false, false> uncompressed(numNodes);
    CSRStorageImm<EagerFileVector, uint32_t, true, true> compressed(numNodes);
    size_t total = 0;
    for (size_t i = 0; i < numNodes; i++) {
        std::vector<uint32_t> edges(numEdges);
        for (size_t j = 0; j < numEdges; j++) {
            edges[j] = sampler(generator);
        }
        std::sort(edges.begin(), edges.end());
        uncompressed.setData(i, edges);
        total += compressed.setData(i, edges);
    }

    std::ofstream outFile("compressed.test.csr", std::ios::binary);
    const uint64_t valueBytes2 = compressed.numValueBytes();
    outFile.write((const char*)&valueBytes2, sizeof(uint64_t));
    ASSERT_EQ(valueBytes2, compressed.flushBuckets(outFile));
    const uint64_t indexBytes2 = compressed.numNodes();
    ASSERT_EQ(indexBytes2+1, compressed.flushIndexes(outFile));
    outFile.close();

    auto infile = std::make_shared<std::ifstream>("compressed.test.csr", std::ios::binary);
    ASSERT_TRUE(infile->good());
    uint64_t valueBytesRd = 0;
    infile->read((char*)&valueBytesRd, sizeof(uint64_t));
    EagerFileVector<uint8_t> edges(infile, sizeof(uint64_t), valueBytesRd);
    EagerFileVector<uint32_t> nodes(infile, sizeof(uint64_t)+valueBytesRd, numNodes+1);
    CSRStorageImm<EagerFileVector, uint32_t, true, true> compressedRead(
        std::move(nodes), std::move(edges), 0);
    for (size_t i = 0; i < numNodes; i++) {
        std::vector<uint32_t> fromRam;
        std::vector<uint32_t> fromDisk;
        uncompressed.getData(i, fromRam);
        compressedRead.getData(i, fromDisk);
        ASSERT_EQ(fromRam.size(), fromDisk.size());
        ASSERT_EQ(fromRam, fromDisk);
    }
}