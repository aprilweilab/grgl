#include <gtest/gtest.h>
#include <stdexcept>

#include "util.h"
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
