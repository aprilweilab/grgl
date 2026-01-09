#include <gtest/gtest.h>
#include <limits>

#include "grg_helpers.h"
#include "grgl/common.h"

using namespace grgl;

TEST(Common, FloatRange) {
    IntRange denorm;
    FloatRange orig;
    ASSERT_TRUE(orig.isUnspecified());
    ASSERT_EQ(orig, orig.normalized(0, 100));
    denorm = orig.denormalized(0, 100);
    ASSERT_EQ(0, denorm.start());
    ASSERT_EQ((size_t)denorm.end(), std::numeric_limits<size_t>::max());

    orig = {0, 1000};
    ASSERT_EQ(orig.start(), 0);
    ASSERT_EQ(orig.end(), 1000);
    denorm = orig.normalized(0, 1000).denormalized(0, 1000);
    ASSERT_EQ((size_t)orig.start(), denorm.start());
    ASSERT_EQ((size_t)orig.end(), denorm.end());

    FloatRange firstHalf = {0, 500};
    FloatRange secondHalf = {500, 1000};
    ASSERT_FALSE(firstHalf.contains(secondHalf.start()));
    FloatRange firstHalfNorm = {0, 0.5};
    FloatRange secondHalfNorm = {0.5, 1.0};
    ASSERT_FALSE(firstHalfNorm.contains(secondHalfNorm.start()));
    IntRange firstHalfDenorm = firstHalfNorm.denormalized(0, 1000);
    IntRange secondHalfDenorm = secondHalfNorm.denormalized(0, 1000);
    ASSERT_FALSE(firstHalfDenorm.contains(secondHalfDenorm.start()));
}

// A range that is not evenly integer divisible should cover every positons
// once and only once.
TEST(Common, UnevenRange) {
    std::vector<FloatRange> normRanges = {
        {0, 0.333}, {0.333, 0.666}, {0.666, 1.0}
    };
    std::vector<IntRange> denormRanges;
    for (const auto& range : normRanges) {
        denormRanges.push_back(range.denormalized(0, 1000));
    }
    for (size_t i = 0; i < 1000; i++) {
        size_t foundTimes = 0;
        for (const auto& range : denormRanges) {
            if (range.contains(i)) {
                foundTimes++;
            }
        }
        ASSERT_EQ(foundTimes, 1);
    }
}
