#include <gtest/gtest.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <unistd.h>
#include <vector>

#include "grgl/mutation.h"
#include "grgl/serialize.h"
#include "varint.h"
#include "common_grgs.h"

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

TEST(Serialization, VarInt) {
    InMemBuffer buffer(DEFAULT_BUFFER_SIZE);
    std::ostream outStream(&buffer);

    writeVarInt(123, outStream);
    writeVarInt(1234, outStream);
    writeVarInt(12345, outStream);
    writeVarInt(123456, outStream);
    writeVarInt(1234567, outStream);
    writeVarInt(12345678, outStream);
    writeVarInt(123456789, outStream);
    writeVarInt(12345678910, outStream);
    writeVarInt(1234567891011, outStream);
    writeVarInt(123456789101112, outStream);

    std::istream inStream(&buffer);
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(123, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(1234, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(12345, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(123456, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(1234567, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(12345678, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(123456789, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(12345678910, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(1234567891011, readVarInt(inStream));
    ASSERT_TRUE(inStream.good());
    ASSERT_EQ(123456789101112, readVarInt(inStream));
}

TEST(Serialization, SimpleRoundTrip) {
    const bool useVarInt = true;
    const bool allowSimplify = false;
    Mutation testMut = {999, "ACCC", "T", 19.1, NO_ORIGINAL_ID};
    const char * const testFile = "test.serialize.roundtrip.grg";
    std::ofstream outStream(testFile);
    MutableGRGPtr grg = depth3BinTree();
    grg->addPopulation("pop1");
    grg->addPopulation("pop2");
    ASSERT_TRUE(grg->numEdges() == 6);
    ASSERT_TRUE(grg->numNodes() == 7);
    grg->addMutation(testMut, 6);
    writeGrg(grg, outStream, useVarInt, allowSimplify);
    outStream.close();

    // Deserialize it and check.
    std::ifstream inStream(testFile);
    GRGPtr grg2 = readImmutableGrg(inStream);
    ASSERT_EQ(grg2->numEdges(), 6);
    ASSERT_EQ(grg2->numNodes(), 7);
    auto actualMuts = grg2->getMutationsForNode(6);
    ASSERT_EQ(actualMuts.size(), 1);
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getAllele(), testMut.getAllele());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getPosition(), testMut.getPosition());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getRefAllele(), testMut.getRefAllele());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getTime(), testMut.getTime());
    ASSERT_EQ(grg2->getMutationById(actualMuts[0]).getOriginalId(), testMut.getOriginalId());
    inStream.close();

    // Serialize/deserialize it again, disable variable-sized int encoding.
    const char * const testFile2 = "test.serialize.roundtrip2.grg";
    std::ofstream outStream2(testFile2);
    writeGrg(grg, outStream2, false, allowSimplify);
    outStream2.close();
    std::ifstream inStream2(testFile2);
    GRGPtr grg3 = readImmutableGrg(inStream2);
    ASSERT_TRUE(grg3->numEdges() == 6);
    ASSERT_TRUE(grg3->numNodes() == 7);
    inStream2.close();

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
    std::istream inStream(&buffer);
    // Buffer has random memory contents; will fail magic number check
    ASSERT_THROW(readMutableGrg(inStream), grgl::SerializationFailure);
}
