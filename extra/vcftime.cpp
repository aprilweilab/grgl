/* Just parse all genotypes in the VCF, don't do anything with the data. Time the file scan.
 *
 * Usage:
 *  vcftime <file>
 */
#include <iostream>

#include "picovcf.hpp"

using namespace picovcf;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Please pass in an input file" << std::endl;
        return 1;
    }

    const std::string filename(argv[1]);

    VCFFile vcf(filename);

    size_t totalSamples = 0;
    vcf.seekBeforeVariants();
    while (vcf.hasNextVariant()) {
        vcf.nextVariant();
        VCFVariantView variant = vcf.currentVariant();
        IndividualIteratorGT iterator = variant.getIndividualIterator();
        while (iterator.hasNext()) {
            VariantT allele1 = 0;
            VariantT allele2 = 0;
            bool isPhased = iterator.getAlleles(allele1, allele2);
            if (!isPhased) {
                std::cerr << "Cannot create a matrix for unphased data" << std::endl;
                return 2;
            }
            totalSamples++;
            if (allele2 != NOT_DIPLOID) {
                totalSamples++;
            }
        }
    }
    std::cout << "Scanned " << totalSamples << std::endl;
    return 0;
}
