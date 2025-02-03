#include <gtest/gtest.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <limits>
#include <streambuf>
#include <unistd.h>
#include <vector>

#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"
#include "common_grgs.h"
#include "grgl/visitor.h"
#include "vbyte.h"

constexpr size_t DEFAULT_BUFFER_SIZE = 1024;

/**
 * Utility class for reading/writing serialized data to memory instead
 * of disk. Just makes testing easier.
 *
 * This doesn't work with seeking, so is only useful for simple tests.
 */
class InMemBuffer : public std::streambuf
{
public:
    explicit InMemBuffer(size_t size)
            : m_buffer(new char[size]) {
        reset(size);
    }

    void reset(size_t size) {
        m_bufSize = size;
        setp(m_buffer, m_buffer + m_bufSize);
        setg(m_buffer, m_buffer, m_buffer + m_bufSize);
    }

    size_t bytesWritten() const {
        return static_cast<size_t>(pptr() - pbase());
    }

    ~InMemBuffer() override {
        delete [] m_buffer;
    }

    char *m_buffer;
    size_t m_bufSize;
};

using namespace grgl;

// Test our wrapper around libvbyte (and libvbyte as well)
TEST(Serialization, VByteDispatcher) {
    uint8_t outBuf[1024];
    {
        VByteDispatcher<uint64_t, true, false> disp;
        const uint64_t big = 0xfffffffffU;
        const uint64_t small = 1009;
        uint64_t result = 0;

        disp.compressSingle(big, &outBuf[0]);
        disp.decompressSingle(&outBuf[0], &result);
        ASSERT_EQ(big, result);
        ASSERT_GE(disp.maxElemSize(), sizeof(uint64_t));

		std::vector<uint64_t> results;
		std::vector<uint64_t> values = {small, big};
        disp.compressData(values.data(), &outBuf[0], values.size());
        results.resize(values.size());
        disp.decompressData(&outBuf[0], results.data(), results.size());
        ASSERT_EQ(values, results);
    }

    {
        VByteDispatcher<uint32_t, true, true> disp;
        const uint32_t big = 0xfffffffU;
        const uint32_t small = 12999;
        uint32_t result = 0;

        disp.compressSingle(big, &outBuf[0]);
        disp.decompressSingle(&outBuf[0], &result);
        ASSERT_EQ(big, result);
        ASSERT_GE(disp.maxElemSize(), sizeof(uint32_t));

		std::vector<uint32_t> results;
		std::vector<uint32_t> values = {small, big};
        disp.compressData(values.data(), &outBuf[0], values.size());
        results.resize(values.size());
        disp.decompressData(&outBuf[0], results.data(), results.size());
        ASSERT_EQ(values, results);
    }

    {
        VByteDispatcher<uint32_t, false, false> disp;
        const uint32_t big = std::numeric_limits<uint32_t>::max();
        const uint32_t small = 1;
        uint32_t result = 0;

        disp.compressSingle(big, &outBuf[0]);
        disp.decompressSingle(&outBuf[0], &result);
        ASSERT_EQ(big, result);
        ASSERT_EQ(disp.maxElemSize(), sizeof(uint32_t));

		std::vector<uint32_t> results;
		std::vector<uint32_t> values = {small, big};
        disp.compressData(values.data(), &outBuf[0], values.size());
        results.resize(values.size());
        disp.decompressData(&outBuf[0], results.data(), results.size());
        ASSERT_EQ(values, results);
    }
}

TEST(Serialization, SimpleRoundTrip) {
    const bool allowSimplify = false;
    Mutation testMut = {999, "ACCC", "T", 19.1};
    const char * const testFile = "test.serialize.roundtrip.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    grg->addPopulation("pop1");
    grg->addPopulation("pop2");
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    grg->addMutation(testMut, 6);
    writeGrg(grg, outStream, allowSimplify);
    outStream.close();

    // Deserialize it and check.
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_EQ(grg2->numEdges(), 6);
    ASSERT_EQ(grg2->numNodes(), 7);
    auto actualMuts = grg2->getMutationsForNode(6);
    ASSERT_EQ(actualMuts.size(), 1);
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getAllele(), testMut.getAllele());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getPosition(), testMut.getPosition());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getRefAllele(), testMut.getRefAllele());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getTime(), testMut.getTime());

    // Serialize/deserialize it again, disable variable-sized int encoding.
    const char * const testFile2 = "test.serialize.roundtrip2.grg";
    std::ofstream outStream2(testFile2);
    writeGrg(grg, outStream2, allowSimplify);
    outStream2.close();
    grgl::IFSPointer inStream2 = std::make_shared<std::ifstream>(testFile2);
    GRGPtr grg3 = readImmutableGrg(inStream2);
    ASSERT_TRUE(grg3->numEdges() == 6);
    ASSERT_TRUE(grg3->numNodes() == 7);

    // We cannot compare serialized buffers, because we make use of unordered
    // sets in the GRG, and do not sort them prior to serialization.
}

TEST(Serialization, WriteFailure) {
    InMemBuffer buffer(4); // Too small for GRG to fit!
    std::ostream outStream(&buffer);
    GRGPtr grg = depth3BinTree();
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    ASSERT_THROW(writeGrg(grg, outStream), grgl::SerializationFailure);
}

TEST(Serialization, ReadFailure) {
    InMemBuffer buffer(DEFAULT_BUFFER_SIZE);
    grgl::IFSPointer inStream = std::make_shared<std::istream>(&buffer);
    // Buffer has random memory contents; will fail magic number check
    ASSERT_THROW(readMutableGrg(inStream), grgl::SerializationFailure);
}

TEST(Serialization, LongAllele) {
    const bool allowSimplify = true;
    Mutation testMut = {1, "this is a really long allele", "and so is this!", 19.1};
    const char * const testFile = "test.serialize.longallele.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    grg->addMutation(testMut, 6);
    writeGrg(grg, outStream, allowSimplify);
    outStream.close();

    // Deserialize it and check.
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_LE(grg2->numEdges(), 6);
    ASSERT_LE(grg2->numNodes(), 7);
    grgl::Mutation mut = grg2->getMutationById(0);
    ASSERT_EQ(mut.getAllele(), testMut.getAllele());
    ASSERT_EQ(mut.getPosition(), testMut.getPosition());
    ASSERT_EQ(mut.getRefAllele(), testMut.getRefAllele());
    ASSERT_EQ(mut.getTime(), testMut.getTime());
}

TEST(Serialization, SubsetMuts) {
    const bool allowSimplify = false;
    Mutation mut1 = {999, "ACCC", "T", 19.1};
    Mutation mut2 = {10001, "C", "A", 20.0};
    const char * const testFile = "test.serialize.subsetmuts.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    grg->addPopulation("pop1");
    grg->addPopulation("pop2");
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    MutationId kept = grg->addMutation(mut1, 4);
    grg->addMutation(mut2, 6);

    GRGOutputFilter filter = {TraversalDirection::DIRECTION_DOWN, {kept}};
    auto result = simplifyAndSerialize(grg, outStream, filter);
    NodeIDSizeT resultNodes = result.first;
    size_t resultEdges = result.second;
    ASSERT_EQ(resultNodes, 5); // The 4 samples + the node w/ mut1
    ASSERT_EQ(resultEdges, 2); // mut1 only reaches 2 samples
    outStream.close();

    // Deserialize it and check.
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_EQ(grg2->numEdges(), 2);
    ASSERT_EQ(grg2->numNodes(), 5);
    ASSERT_EQ(grg2->numMutations(), 1);
    Mutation loadedMut = grg2->getMutationById(0);
    ASSERT_EQ(loadedMut.getPosition(), 999);
    for (const auto& nodeMuts : grg2->getNodeMutationPairs()) {
        auto children = grg2->getDownEdges(nodeMuts.first);
        ASSERT_EQ(children.size(), 2);
        NodeIDList expected = {0, 1};
        ASSERT_EQ(children, expected);
    }
    ASSERT_EQ(grg2->getPopulations().size(), 2);
}

TEST(Serialization, SubsetSamples) {
    const bool allowSimplify = false;
    Mutation mut1 = {999, "ACCC", "T", 19.1};
    Mutation mut2 = {10001, "C", "A", 20.0};
    const char * const testFile = "test.serialize.subsetsamples.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    grg->addPopulation("pop1");
    grg->addPopulation("pop2");
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    grg->addMutation(mut1, 4);
    grg->addMutation(mut2, 6);

	// Only keep samples 2 and 3, which means mut1 should become unreachable.
    GRGOutputFilter filter = {TraversalDirection::DIRECTION_UP, {2,3}};
    auto result = simplifyAndSerialize(grg, outStream, filter);
    NodeIDSizeT resultNodes = result.first;
    size_t resultEdges = result.second;
    ASSERT_EQ(resultNodes, 3); // The 2 samples kept + the two nodes above them get simplified into 1
    ASSERT_EQ(resultEdges, 2);
    outStream.close();

    // Deserialize it and check.
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_EQ(grg2->numEdges(), 2);
    ASSERT_EQ(grg2->numNodes(), 3);
    ASSERT_EQ(grg2->numMutations(), 2);
    Mutation loadedMut = grg2->getMutationById(0);
    ASSERT_EQ(loadedMut.getPosition(), 999);
    loadedMut = grg2->getMutationById(1);
    ASSERT_EQ(loadedMut.getPosition(), 10001);
    for (const auto& nodeMuts : grg2->getNodeMutationPairs()) {
        if (nodeMuts.second == 0) {
            ASSERT_EQ(nodeMuts.first, INVALID_NODE_ID);
        } else {
            auto children = grg2->getDownEdges(nodeMuts.first);
            ASSERT_EQ(children.size(), 2);
            NodeIDList expected = {0, 1}; // Samples got renumbered.
            ASSERT_EQ(children, expected);
        }
    }
}

TEST(Serialization, SubsetSamples2) {
    const bool allowSimplify = false;
    Mutation mut1 = {999, "ACCC", "T", 19.1};
    Mutation mut2 = {10001, "C", "A", 20.0};
    const char * const testFile = "test.serialize.subsetsamples.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    grg->addPopulation("pop1");
    grg->addPopulation("pop2");
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    grg->addMutation(mut1, 4);
    grg->addMutation(mut2, 6);

	// Only keep samples 0 and 1, which means both muts are reachable
    GRGOutputFilter filter = {TraversalDirection::DIRECTION_UP, {0,1}};
    auto result = simplifyAndSerialize(grg, outStream, filter);
    NodeIDSizeT resultNodes = result.first;
    size_t resultEdges = result.second;
    ASSERT_EQ(resultNodes, 4); // The 2 samples kept + the two nodes above them
    ASSERT_EQ(resultEdges, 3);
    outStream.close();

    // Deserialize it and check.
    grgl::IFSPointer inStream = std::make_shared<std::ifstream>(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_EQ(grg2->numEdges(), 3);
    ASSERT_EQ(grg2->numNodes(), 4);
    ASSERT_EQ(grg2->numMutations(), 2);
    Mutation loadedMut = grg2->getMutationById(0);
    ASSERT_EQ(loadedMut.getPosition(), 999);
    loadedMut = grg2->getMutationById(1);
    ASSERT_EQ(loadedMut.getPosition(), 10001);
    for (const auto& nodeMuts : grg2->getNodeMutationPairs()) {
        ASSERT_NE(nodeMuts.first, INVALID_NODE_ID);
        auto children = grg2->getDownEdges(nodeMuts.first);
        if (nodeMuts.second == 1) {
            ASSERT_EQ(children.size(), 1);
        } else {
            ASSERT_EQ(children.size(), 2);
            NodeIDList expected = {0, 1}; // Samples got renumbered.
            ASSERT_EQ(children, expected);
        }
    }
}