#include <gtest/gtest.h>

#include "grgl/mutation.h"

using namespace grgl;

TEST(Mutation, Allelees) {
    {
        Mutation m1;
        ASSERT_TRUE(m1.isEmpty());
        ASSERT_EQ(m1.getAllele(), "");
        ASSERT_EQ(m1.getRefAllele(), "");
    }

    {
        Mutation m2(1, "A");
        ASSERT_TRUE(!m2.isEmpty());
        ASSERT_EQ(m2.getAllele(), "A");
        ASSERT_EQ(m2.getRefAllele(), "");
        ASSERT_EQ(m2.getPosition(), 1);
    }

    {
        Mutation m3(0xFFFFFFFF, "ACGTGTT");
        ASSERT_TRUE(!m3.isEmpty());
        ASSERT_EQ(m3.getAllele(), "ACGTGTT");
        ASSERT_EQ(m3.getRefAllele(), "");
        ASSERT_EQ(m3.getPosition(), 0xFFFFFFFF);
    }

    {
        Mutation m4(0xFFFFFFFFFFFFFFFF-1, "ACGTGTTAAAAA");
        ASSERT_TRUE(!m4.isEmpty());
        ASSERT_EQ(m4.getAllele(), "ACGTGTTAAAAA");
        ASSERT_EQ(m4.getRefAllele(), "");
        ASSERT_EQ(m4.getPosition(), 0xFFFFFFFFFFFFFFFF-1);
    }

    {
        // Yep, this was a bug when I used float...
        ASSERT_TRUE(((uint64_t)(float)18002039) != 18002039);
        Mutation m5(18002039, "ACGTGTTAAAAA", "AAAAAAAAAAAAAA");
        ASSERT_TRUE(!m5.isEmpty());
        ASSERT_EQ(m5.getAllele(), "ACGTGTTAAAAA");
        ASSERT_EQ(m5.getRefAllele(), "AAAAAAAAAAAAAA");
        ASSERT_EQ(m5.getPosition(), 18002039);
    }

    {
        Mutation m6(9, "A", "AAAAAAAAAAAAAG");
        ASSERT_TRUE(!m6.isEmpty());
        ASSERT_EQ(m6.getAllele(), "A");
        ASSERT_EQ(m6.getRefAllele(), "AAAAAAAAAAAAAG");
        ASSERT_EQ(m6.getPosition(), 9);
    }
}
