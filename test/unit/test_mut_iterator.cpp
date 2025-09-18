#include <gtest/gtest.h>

#include "grgl/common.h"
#include "grgl/mut_iterator.h"
#include "picovcf.hpp"
#include "testing_utilities.h"

using namespace grgl;

TEST(MutationIterator, VCFMultipleContigs) {
    std::stringstream ssVCF;
    ssVCF << "##fileformat=VCFv4.2\n";
    ssVCF << "##source=Testing\n";
    ssVCF << "##contig=<ID=20,length=10000>\n";
    ssVCF << "##contig=<ID=blargh,length=20000>\n";
    ssVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n";
    ssVCF << "20\t4000\t.\tG\tA\t.\tPASS\t.GT\t0|0\t1|1\t0|1\n";
    ssVCF << "20\t4500\t.\tT\tC\t.\tPASS\t.GT\t1|1\t1|1\t1|1\n";
    ssVCF << "20\t7777\t.\tA\tC\t.\tPASS\t.GT\t1|1\t1|1\t1|1\n";
    ssVCF << "blargh\t14000\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|1\n";
    ssVCF << "blargh\t14500\t.\tT\tC\t.\tPASS\t.\tGT\t1|0\t1|1\t1|1\n";
    ssVCF << "blargh\t17777\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1\n";
    const auto filename = writeTempFile(ssVCF.str(), ".vcf");
    // If we don't FORCE, then a plain VCF will fail
    ASSERT_THROW(auto it = makeMutationIterator(filename, {}, MIT_FLAG_EMPTY),
                 ApiMisuseFailure);
    // This should fail because of multiple contigs, but only after we scan the file.
    auto it = makeMutationIterator(filename, {}, MIT_FLAG_FORCE);
    ASSERT_THROW(it->countMutations(), picovcf::ApiMisuse);
}

TEST(MutationIterator, VCFValid) {
    std::stringstream ssVCF;
    ssVCF << "##fileformat=VCFv4.2\n";
    ssVCF << "##source=Testing\n";
    ssVCF << "##contig=<ID=20,length=10000>\n";
    ssVCF << "##contig=<ID=blargh,length=20000>\n";
    ssVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n";
    ssVCF << "20\t4000\t.\tG\tA\t.\tPASS\t.\tGT\t0|0\t1|1\t0|1\n";
    ssVCF << "20\t4500\t.\tT\tC\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\n";
    ssVCF << "20\t7777\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\n";
    ssVCF << "20\t14000\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|1\n";
    ssVCF << "20\t14500\t.\tT\tC\t.\tPASS\t.\tGT\t1|0\t1|1\t1|1\n";
    ssVCF << "20\t17777\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1";
    const auto filename = writeTempFile(ssVCF.str(), ".vcf");
    auto it = makeMutationIterator(filename, {}, MIT_FLAG_FORCE);
    ASSERT_EQ(6, it->countMutations());
    size_t total = 0;
    std::vector<size_t> counts;
    MutationAndSamples item;
    while (it->next(item, total)) {
        ASSERT_GE(item.mutation.getPosition(), 4000);
        ASSERT_LE(item.mutation.getPosition(), 17777);
        counts.push_back(item.samples.size());
    }
    std::vector<size_t> expected = {3, 6, 6, 3, 5, 3};
    ASSERT_EQ(expected, counts);

    auto it2 = makeMutationIterator(filename, {4500, 14501}, MIT_FLAG_FORCE);
    ASSERT_EQ(4, it2->countMutations());
}

TEST(MutationIterator, VCFWithMissing) {
    std::stringstream ssVCF;
    ssVCF << "##fileformat=VCFv4.2\n";
    ssVCF << "##source=Testing\n";
    ssVCF << "##contig=<ID=20,length=10000>\n";
    ssVCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n";
    ssVCF << "20\t4000\t.\tG\tA\t.\tPASS\t.\tGT\t0|0\t1|1\t0|.\n";
    ssVCF << "20\t4500\t.\tT\tC\t.\tPASS\t.\tGT\t1|1\t1|.\t.|1\n";
    ssVCF << "20\t7777\t.\tA\tC\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\n";
    ssVCF << "20\t14000\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\t0|1\t0|1\n";
    ssVCF << "20\t14500\t.\tT\tC\t.\tPASS\t.\tGT\t1|0\t1|1\t1|1\n";
    const auto filename = writeTempFile(ssVCF.str(), ".vcf");
    auto it = makeMutationIterator(filename, {}, MIT_FLAG_FORCE);
    ASSERT_EQ(5, it->countMutations());
    MutationAndSamples item;
    size_t total = 0;
    size_t lastPosition = 0;
    size_t numMissing = 0;
    while (it->next(item, total)) {
        if (item.mutation.isMissing()) {
            numMissing++;
            // This asserts that missing data _always comes first_ in the mutation order, by position
            ASSERT_GT(item.mutation.getPosition(), lastPosition);
        }
        ASSERT_GE(item.mutation.getPosition(), 4000);
        ASSERT_LE(item.mutation.getPosition(), 14500);
        lastPosition = item.mutation.getPosition();
    }

    auto it2 = makeMutationIterator(filename, {4500, 14001}, MIT_FLAG_FORCE);
    ASSERT_EQ(3, it2->countMutations());
}
