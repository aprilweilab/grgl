#include <gtest/gtest.h>
#include <stdexcept>

#include "util.h"
#include "hap_helpers.h"
#include "testing_utilities.h"

TEST(Util, MapFromTSVGood) {
    std::string testString = 
        "sample\tpop\tsuper_pop\trandom\n"
        "A\tPOP1\tSPOP5\t2342830492\n"
        "B\tPOP2\tSPOP3\tksdjflksjdf\n"
    ;
    auto filename = writeTempFile(testString);
    auto resultMap = loadMapFromTSV(filename, "sample", "pop");
    ASSERT_EQ(resultMap.size(), 2);
    ASSERT_NE(resultMap.find("A"), resultMap.end());
    ASSERT_NE(resultMap.find("B"), resultMap.end());
    ASSERT_EQ(resultMap.find("A")->second, "POP1");
    ASSERT_EQ(resultMap.find("B")->second, "POP2");

    auto otherResultMap = loadMapFromTSV(filename, "random", "super_pop");
    ASSERT_EQ(otherResultMap.size(), 2);
    ASSERT_EQ(otherResultMap.find("ksdjflksjdf")->second, "SPOP3");
    ASSERT_EQ(otherResultMap.find("2342830492")->second, "SPOP5");

	remove_file(filename);
}

TEST(Util, MapFromTSVBad) {
    std::string testString = 
        "sample\tpop\tsuper_pop\trandom\n"
        "A\tPOP1\tSPOP5\t2342830492\n"
        "B\tPOP2\tksdjflksjdf\n"
    ;
    auto filename = writeTempFile(testString);
    EXPECT_THROW(loadMapFromTSV(filename, "sample", "pop"), std::runtime_error);
	remove_file(filename);
}

TEST(HapHelpers, getBitsAsList) {
    grgl::HaplotypeVector hapVect(10);
    grgl::setBit(hapVect, 10);
    grgl::setBit(hapVect, 100);
    grgl::setBit(hapVect, 99);
    grgl::setBit(hapVect, 33);
    grgl::setBit(hapVect, (10 * 8 * 4) - 1);
    auto result = grgl::getBitsAsList(hapVect);
    std::vector<size_t> expected = {10, 33, 99, 100, 319};
    ASSERT_EQ(result, expected);
}


TEST(HapHelpers, bitwiseSubtract) {
    grgl::HaplotypeVector hap1(10);
    for (size_t bit : {10, 100, 99, 33}) {
        grgl::setBit(hap1, bit);
    }
    grgl::HaplotypeVector hap2(10);
    for (size_t bit : {10, 77, 99, 11, 4, 81}) {
        grgl::setBit(hap2, bit);
    }

    size_t count = grgl::bitwiseSubtract(hap2, hap1);
    ASSERT_EQ(count, 4);
    auto bitsSet = grgl::getBitsAsList(hap2);
    std::vector<size_t> expected = {4, 11, 77, 81};
    ASSERT_EQ(bitsSet, expected);
}

TEST(Util, vectReverse) {
    grgl::NodeIDList theVect;

    // Edge case: empty
    vectReverse(theVect);
    ASSERT_TRUE(theVect.empty());

    // Edge case: 1 element
    theVect = {91919191};
    vectReverse(theVect);
    ASSERT_EQ(theVect[0], 91919191);

    // Edge case: 2 elements
    theVect = {1001, 999};
    vectReverse(theVect);
    ASSERT_EQ(theVect[0], 999);
    ASSERT_EQ(theVect[1], 1001);

    // Odd # of elements
    theVect = {1, 9, 2};
    vectReverse(theVect);
    ASSERT_EQ(theVect[0], 2);
    ASSERT_EQ(theVect[1], 9);
    ASSERT_EQ(theVect[2], 1);

    // Even # of elements
    theVect = {100, 66, 87, 33};
    vectReverse(theVect);
    ASSERT_EQ(theVect[0], 33);
    ASSERT_EQ(theVect[1], 87);
    ASSERT_EQ(theVect[2], 66);
    ASSERT_EQ(theVect[3], 100);

}