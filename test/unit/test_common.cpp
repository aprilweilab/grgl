#include <gtest/gtest.h>

#include "grgl/common.h"

using namespace grgl;

TEST(Common, FloatRange) {
    FloatRange orig;
    ASSERT_TRUE(orig.isUnspecified());
    ASSERT_EQ(orig, orig.normalized(0, 100));
    ASSERT_EQ(orig, orig.denormalized(0, 100));

    orig = {0, 1000};
    ASSERT_EQ(orig.start(), 0);
    ASSERT_EQ(orig.end(), 1000);
    ASSERT_EQ(orig, orig.normalized(0, 1000).denormalized(0, 1000));
    FloatRange firstHalf = {0, 500};
    FloatRange secondHalf = {500, 1000};
    ASSERT_FALSE(firstHalf.contains(secondHalf.start()));
    FloatRange firstHalfNorm = {0, 0.5};
    FloatRange secondHalfNorm = {0.5, 1.0};
    ASSERT_FALSE(firstHalfNorm.contains(secondHalfNorm.start()));
    FloatRange firstHalfDenorm = firstHalfNorm.denormalized(0, 1000);
    FloatRange secondHalfDenorm = secondHalfNorm.denormalized(0, 1000);
    ASSERT_FALSE(firstHalfDenorm.contains(secondHalfDenorm.start()));
}
