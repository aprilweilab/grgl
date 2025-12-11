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

TEST(Mutation, Positions) {
    ASSERT_TRUE(bpInRange({0, 5000}, 4999));
    ASSERT_TRUE(bpInRange({0, 5000}, 0));
    ASSERT_FALSE(bpInRange({0, 5000}, 5000));

    ASSERT_FALSE(bpInRangeExceptEdges({0, 5000}, 4999));
    ASSERT_FALSE(bpInRangeExceptEdges({0, 5000}, 0));
    ASSERT_FALSE(bpInRangeExceptEdges({0, 5000}, 5000));
    ASSERT_TRUE(bpInRangeExceptEdges({0, 5000}, 1));
    ASSERT_TRUE(bpInRangeExceptEdges({0, 5000}, 4998));

    ASSERT_TRUE(bpOverlap({50, 100}, {98, 1000}));
    ASSERT_TRUE(bpOverlap({98, 1000}, {50, 100}));
    ASSERT_FALSE(bpOverlap({50, 100}, {99, 1000}));
    ASSERT_FALSE(bpOverlap({99, 1000}, {50, 100}));
}