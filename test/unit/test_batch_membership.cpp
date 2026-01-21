#include <gtest/gtest.h>

#include <vector>

#include "grgl/map_mutations.h"

using namespace grgl;

TEST(BatchMembership, SetClear) {
    const size_t bitCount = 70;
    std::vector<uint32_t> wordStorage((bitCount + 31) / 32, 0xDEADBEEFu);
    BatchMembership membership(wordStorage.data(), bitCount);

    membership.setOnes(50);
    EXPECT_TRUE(membership.any());
    EXPECT_EQ(wordStorage[0], 0xFFFFFFFFu);
    EXPECT_EQ(wordStorage[1], (1u << 18) - 1u); // only first 18 bits set in second word
    EXPECT_EQ(wordStorage[2], 0u);              // remainder cleared

    membership.setOnes(0);
    EXPECT_FALSE(membership.any());
    for (auto word : wordStorage) {
        EXPECT_EQ(word, 0u);
    }
}

TEST(BatchMembership, Iterate) {
    const size_t bitCount = 70;
    std::vector<uint32_t> wordStorage((bitCount + 31) / 32, 0);
    BatchMembership membership(wordStorage.data(), bitCount);

    membership.setBit(0);
    membership.setBit(31);
    membership.setBit(32);
    membership.setBit(69); // last valid bit
    EXPECT_TRUE(membership.any());

    std::vector<size_t> setBitIndices;
    for (size_t bitIndex : membership) {
        setBitIndices.push_back(bitIndex);
    }
    EXPECT_EQ((std::vector<size_t>{0, 31, 32, 69}), setBitIndices);
}

TEST(BatchMembership, AndDiffCopy) {
    const size_t bitCount = 64;
    std::vector<uint32_t> leftStorage((bitCount + 31) / 32, 0);
    std::vector<uint32_t> rightStorage((bitCount + 31) / 32, 0);
    BatchMembership leftMembership(leftStorage.data(), bitCount);
    BatchMembership rightMembership(rightStorage.data(), bitCount);

    leftMembership.setBit(1);
    leftMembership.setBit(5);
    leftMembership.setBit(63);
    rightMembership.setBit(5);
    rightMembership.setBit(9);
    rightMembership.setBit(63);

    leftMembership.andAssign(rightMembership);
    std::vector<size_t> intersectionBits;
    for (size_t bitIndex : leftMembership) {
        intersectionBits.push_back(bitIndex);
    }
    EXPECT_EQ((std::vector<size_t>{5, 63}), intersectionBits);

    std::vector<uint32_t> copyStorage((bitCount + 31) / 32, 0);
    BatchMembership copiedMembership(copyStorage.data(), bitCount);
    copiedMembership.copyFrom(leftMembership);
    std::vector<size_t> copiedBitIndices;
    for (size_t bitIndex : copiedMembership) {
        copiedBitIndices.push_back(bitIndex);
    }
    EXPECT_EQ(intersectionBits, copiedBitIndices);
}

TEST(BatchMembership, Diff) {
    const size_t bitCount = 64;
    std::vector<uint32_t> baseStorage((bitCount + 31) / 32, 0);
    std::vector<uint32_t> otherStorage((bitCount + 31) / 32, 0);
    BatchMembership baseMembership(baseStorage.data(), bitCount);
    BatchMembership otherMembership(otherStorage.data(), bitCount);

    baseMembership.setBit(1);
    baseMembership.setBit(5);
    baseMembership.setBit(8);
    baseMembership.setBit(63);
    otherMembership.setBit(5);
    otherMembership.setBit(9);
    otherMembership.setBit(63);

    std::vector<size_t> diffBits;
    for (auto diffs = baseMembership.beginDiff(otherMembership); diffs != baseMembership.endDiff(otherMembership); ++diffs) {
        diffBits.push_back(*diffs);
    }
    EXPECT_EQ((std::vector<size_t>{1, 8, 9}), diffBits);
}

