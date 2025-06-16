#include <gtest/gtest.h>
#include <stdexcept>
#include <random>

#include "grgl/grg.h"
#include "build_shape.h"
#include "testing_utilities.h"

using namespace grgl;

inline std::string randomVcfLine(size_t position, size_t numIndivs) {
    static std::random_device randDevice;
    static std::mt19937 generator(randDevice());

    std::uniform_int_distribution<size_t> fnSampler(0, 3);
    std::stringstream ss;
    ss << "1\t" << position << "\tV1\tA\tC\t.\tPASS\t.\tGT";
    for (size_t i = 0; i < numIndivs; i++) {
        size_t allele = fnSampler(generator);
        ss << "\t" << (allele & 0x1U) << "|" << ((allele & 0x2U) >> 1U);
    }
    ss << "\n";
    return ss.str();
}

TEST(Construct, WithPopIds) {
    std::string vcfHeaderString = 
        "##fileformat=VCFv4.2\n"
        "##source=TEST\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX1\tX2\tX3\tA1\tA2\tB4\tZ1\n"
    ;
    const size_t numIndividuals = 7;
    std::stringstream vcfFileSS;
    vcfFileSS << vcfHeaderString;
    vcfFileSS << randomVcfLine(60000, numIndividuals);
    vcfFileSS << randomVcfLine(60001, numIndividuals);
    vcfFileSS << randomVcfLine(60002, numIndividuals);
    vcfFileSS << randomVcfLine(60003, numIndividuals);
    vcfFileSS << randomVcfLine(60004, numIndividuals);
    vcfFileSS << randomVcfLine(60005, numIndividuals);
    vcfFileSS << randomVcfLine(60006, numIndividuals);
    auto filename = writeTempFile(vcfFileSS.str(), ".vcf");
    std::map<std::string, std::string> indivIdToPop;
    FloatRange fullRange;

    // Test1: Incomplete individual -> population map
    indivIdToPop.emplace("Z1", "Population2");
    indivIdToPop.emplace("X1", "Population4");
    EXPECT_THROW(createEmptyGRGFromSamples(filename, fullRange, 8, 0x0U, 0.0,
                                           indivIdToPop), std::runtime_error);

    // Test2: Complete individual -> population map
    indivIdToPop.emplace("A2", "Population3");
    indivIdToPop.emplace("X2", "Population1");
    indivIdToPop.emplace("B4", "Population3");
    indivIdToPop.emplace("A1", "Population1");
    indivIdToPop.emplace("X3", "Population4");
    auto grg = createEmptyGRGFromSamples(filename, fullRange, 8, 0x0U, 0.0,
                                         indivIdToPop);
    ASSERT_EQ(grg->numSamples(), numIndividuals*2);
    auto popDescriptions = grg->getPopulations();
    ASSERT_EQ(popDescriptions.size(), 4);
    ASSERT_EQ(popDescriptions[0], "Population4"); // Based on order of individuals in file
    // Test a bunch of sample nodes
    ASSERT_EQ(popDescriptions.at(grg->getPopulationId(0)), "Population4");
    ASSERT_EQ(popDescriptions.at(grg->getPopulationId(12)), "Population2");
    ASSERT_EQ(popDescriptions.at(grg->getPopulationId(13)), "Population2");
    ASSERT_EQ(popDescriptions.at(grg->getPopulationId(8)), "Population3");

	remove_file(filename);
}
