/* Utility for timing the loading of genotype data.
 *
 * Usage:
 *  gtime <infile name>
 */
#include <args.hxx>
#include <chrono>
#include <iostream>
#include <limits>
#include <random>

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mut_iterator.h"
#include "picovcf.hpp"
#include "util.h"

static void mutationIteratorScan(const std::string& inFilename,
                                 bool verbose = true) {
    constexpr size_t EMIT_EVERY = 10000;

    std::shared_ptr<grgl::MutationIterator> iterator = grgl::makeMutationIterator(
        inFilename,
        /*genomeRange=*/{},
        /*binaryMutations=*/false,
        /*emitMissingData=*/false,
        /*flipRefMajor=*/false);
    grgl::MutationAndSamples mutAndSamples;
    size_t variantCount = 0;
    size_t _unused = 0;
    while (iterator->next(mutAndSamples, _unused)) {
        variantCount++;
        if (verbose && (variantCount % EMIT_EVERY == 0)) {
            std::cout << "Variants completed: " << variantCount << std::endl;
        }
    }
}


int main(int argc, char *argv[]) {
    args::ArgumentParser parser("Convert between input file formats.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "infile", "The input file (.vcf or .bgen)");
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (!infile) {
        std::cout << parser;
        return 0;
    }

    mutationIteratorScan(*infile);
    return 0;
}
