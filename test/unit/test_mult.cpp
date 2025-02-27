#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

#include "grgl/grg.h"

TEST(Mult, vectorAddDouble) {
    std::vector<double> v1(8, 0);
    std::vector<double> v2 = {9.0, 1.111, 3.14159, 100000, 88, 100000000.1, 99, 2.7};
    // Only add the first four values.
    grgl::vectorAdd<double>(v1.data(), v2.data(), 4);
    for (size_t i = 4; i < 8; i++) {
        ASSERT_EQ(v1[i], 0.0);
    }
    for (size_t i = 0; i < 4; i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double>(v1.data(), v2.data(), 5);
    ASSERT_EQ(v1[3], v2[3]);
    ASSERT_EQ(v1[4], v2[4]);
    ASSERT_EQ(v1[5], 0.0);

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double>(v1.data(), v2.data(), 6);
    ASSERT_EQ(v1[4], v2[4]);
    ASSERT_EQ(v1[5], v2[5]);
    ASSERT_EQ(v1[6], 0.0);

    v1 = std::move(std::vector<double>(8, 0));
    grgl::vectorAdd<double>(v1.data(), v2.data(), 8);
    grgl::vectorAdd<double>(v1.data(), v2.data(), 8);
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_NEAR(v1[i], v2[i]*2, 0.0001);
    }
}

TEST(Mult, vectorAddFloat) {
    std::vector<float> v1(16, 0);
    std::vector<float> v2 = {9.0, 1.111, 3.14159, std::numeric_limits<float>::max(), 88, 1000.1, 99, 2.7,
                             1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<float>(v1.data(), v2.data(), v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }

    grgl::vectorAdd<float>(v1.data(), v2.data(), 1);
    ASSERT_NEAR(v1[0], 2*v2[0], 0.0001);
    ASSERT_EQ(v1[1], v2[1]);
}

TEST(Mult, vectorAddInt64) {
    std::vector<int64_t> v1(15, 0);
    std::vector<int64_t> v2 = {16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<int64_t>(v1.data(), v2.data(), v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }
}

TEST(Mult, vectorAddInt32) {
    std::vector<int32_t> v1(15, 0);
    std::vector<int32_t> v2 = {16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2};
    ASSERT_EQ(v1.size(), v2.size());
    grgl::vectorAdd<int32_t>(v1.data(), v2.data(), v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        ASSERT_EQ(v1[i], v2[i]);
    }
}