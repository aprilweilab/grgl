#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

#include "grgl/grg.h"

TEST(Mult, vectorAddDouble) {
    std::vector<double> v1(8, 0);
    std::vector<double> v2 = {9.0, 1.111, 3.14159, 100000, 88, 100000000.1, 99, 2.7};
    // Only add the first four values.
    grgl::vectorAdd<double, false>(v1.data(), v2.data(), 0, 0, 4);
    for (size_t i = 4; i < 8; i++) {
        ASSERT_EQ(v1[i], 0.0);
    }
    for (size_t i = 0; i < 4; i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double, false>(v1.data(), v2.data(), 0, 0, 5);
    ASSERT_EQ(v1[3], v2[3]);
    ASSERT_EQ(v1[4], v2[4]);
    ASSERT_EQ(v1[5], 0.0);

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double, false>(v1.data(), v2.data(), 0, 0, 6);
    ASSERT_EQ(v1[4], v2[4]);
    ASSERT_EQ(v1[5], v2[5]);
    ASSERT_EQ(v1[6], 0.0);

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double, false>(v1.data(), v2.data(), 0, 0, 8);
    grgl::vectorAdd<double, false>(v1.data(), v2.data(), 0, 0, 8);
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_NEAR(v1[i], v2[i]*2, 0.0001);
    }
}

TEST(Mult, vectorAddFloat) {
    std::vector<float> v1(16, 0);
    std::vector<float> v2 = {9.0, 1.111, 3.14159, std::numeric_limits<float>::max(), 88, 1000.1, 99, 2.7,
                             1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<float, false>(v1.data(), v2.data(), 0, 0, v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }

    grgl::vectorAdd<float, false>(v1.data(), v2.data(), 0, 0, 1);
    ASSERT_NEAR(v1[0], 2*v2[0], 0.0001);
    ASSERT_EQ(v1[1], v2[1]);
}

TEST(Mult, vectorAddInt64) {
    std::vector<int64_t> v1(15, 0);
    std::vector<int64_t> v2 = {16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<int64_t, false>(v1.data(), v2.data(), 0, 0, v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }
}

TEST(Mult, vectorAddInt32) {
    std::vector<int32_t> v1(15, 0);
    std::vector<int32_t> v2 = {16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<int32_t, false>(v1.data(), v2.data(), 0, 0, v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }
}

TEST(Mult, matmulPerformIOAddition) {
    bool testInput[20] = {false};
    for (size_t i = 0; i < 20; i++) {
        ASSERT_EQ(testInput[i], false);
    }
    // 20 * 8 == 160 bits.
    uint8_t bitvector[20] = {0};

    // Bitwise IO test: round-trip a single bit flip. We set the 5'th booled input, and then map
    // it to the 99th bit. Then we map the 99th
    testInput[5] = true;
    grgl::matmulPerformIOAddition<uint8_t, bool, true>(&bitvector[0], 99, &testInput[0], 5);
    grgl::matmulPerformIOAddition<bool, uint8_t, true>(&testInput[0], 15, &bitvector[0], 99);
    ASSERT_EQ(testInput[5], true);
    ASSERT_EQ(testInput[15], true);
    ASSERT_EQ(testInput[19], false);
}

TEST(Mult, vectorAddBitvector) {
    uint8_t bitvector1[20] = {0};
    uint8_t bitvector2[20] = {0};

    // Add (xor) bv1[16-23] and bv2[24-31] together.
    bitvector1[2] = 0xFF;
    bitvector2[3] = 0x00;
    grgl::vectorAdd<uint8_t, true>(&bitvector1[0], &bitvector2[0], 16, 24, 8);
    ASSERT_EQ(bitvector1[2], 0xFF);

    bitvector1[2] = 0xFF;
    bitvector2[4] = 0xFF;
    grgl::vectorAdd<uint8_t, true>(&bitvector1[0], &bitvector2[0], 16, 32, 8);
    ASSERT_EQ(bitvector1[2], 0x00);

    bitvector1[2] = 0xF0;
    bitvector2[9] = 0x0F;
    grgl::vectorAdd<uint8_t, true>(&bitvector1[0], &bitvector2[0], 16, 9*8, 8);
    ASSERT_EQ(bitvector1[2], 0xFF);

    bitvector1[9] = 0xFF;
    bitvector2[3] = 0x0F;
    grgl::vectorAdd<uint8_t, true>(&bitvector1[0], &bitvector2[0], 9*8, 24, 8);
    ASSERT_EQ(bitvector1[9], 0xF0);
}