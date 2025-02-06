/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2024 April Wei
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <args.hxx>
#include <chrono>
#include <iostream>
#include <string>

#include "calculations.h"
#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "util.h"

int main(int argc, char** argv) {
    auto operationStartTime = std::chrono::high_resolution_clock::now();
#define START_TIMING_OPERATION() operationStartTime = std::chrono::high_resolution_clock::now();
#define EMIT_TIMING_MESSAGE(msg)                                                                                       \
    do {                                                                                                               \
        std::cerr << msg                                                                                               \
                  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - \
                                                                           operationStartTime)                         \
                         .count()                                                                                      \
                  << " ms" << std::endl;                                                                               \
    } while (0)

    args::ArgumentParser parser("Process GRG files.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "grg_file", "The input GRG file");
    args::Flag showStats(parser, "show-stats", "Show statistics about the GRG", {'s', "show-stats"});
    args::Flag alleleFrequency(parser, "allele-frequency", "Calculate allele frequencies", {'f', "freq"});
    args::Flag beta(parser, "association-study", "Conduct GWAS association study", {'a', "association-study"});
    args::ValueFlag<std::string> phenotype(
        parser, "phenotype", "The input phenotype file for association study", {'p', "phenotype"});
    args::Flag onlyCalcBeta(
        parser, "only-beta", "Only calculate the beta value, not all GWAS stats", {'b', "only-beta"});
    args::Flag zygosityInfo(
        parser, "zygosity-info", "Calculate hetero/homozygosity ratios for each allele", {'z', "zygosity-info"});
    args::ValueFlag<std::string> compareGrg(
        parser, "compare", "Compare the input GRG to the given GRG", {'c', "compare"});
    args::ValueFlag<std::string> sampleSubset(
        parser, "sample-subset", "Use only the given file to restrict the subset of samples", {'i', "sample-subset"});
    args::ValueFlag<std::string> region(
        parser, "region", "Use only the given region of the genome, lower-upper.", {'r', "region"});
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

    std::cout << std::fixed << std::setprecision(4);

    grgl::GRGPtr theGRG;
    START_TIMING_OPERATION();
    if (ends_with(*infile, ".grg")) {
        // We only need up edges if we're calculating subsets via samples.
        const bool loadUpEdges = (bool)sampleSubset;
        theGRG = grgl::loadImmutableGRG(*infile, loadUpEdges);
        if (!theGRG) {
            std::cerr << "Failed to load " << *infile << std::endl;
            return 2;
        }
    } else {
        std::cerr << "Unsupported/undetected filetype for " << *infile << std::endl;
        std::cerr << "Only .grg files are supported." << std::endl;
        return 1;
    }
    EMIT_TIMING_MESSAGE("Construction took ");

    std::pair<uint32_t, uint32_t> bpRange = {0, 0};
    if (region) {
        auto tokens = split(*region, '-');
        release_assert(tokens.size() == 2);
        if (!parseExactUint32(tokens[0], bpRange.first)) {
            std::cerr << "Could not parse range arg as int value." << std::endl;
            return 1;
        }
        if (!parseExactUint32(tokens[1], bpRange.second)) {
            std::cerr << "Could not parse range arg as int value." << std::endl;
            return 1;
        }
        release_assert(bpRange.second > bpRange.first);
    }

    grgl::NodeIDList onlySamples;
    if (sampleSubset) {
        onlySamples = loadNodeIDs(*sampleSubset);
    }

    if (showStats) {
        dumpStats(theGRG);
    }

    if (alleleFrequency) {
        START_TIMING_OPERATION();
        emitAlleleFrequency(theGRG, std::cout, bpRange, onlySamples);
        EMIT_TIMING_MESSAGE("Allele frequency calculation took");
    }

    if (zygosityInfo) {
        START_TIMING_OPERATION();
        emitZygosityInfo(theGRG, std::cout, bpRange, onlySamples);
        EMIT_TIMING_MESSAGE("Zygosity information calculation took");
    }

    if (beta) {
        if (!phenotype) {
            std::cerr << "Using random phenotype since no phenotype file specified" << std::endl;
        }
        START_TIMING_OPERATION();
        emitBeta(theGRG, phenotype ? *phenotype : USE_RANDOM_PHENOTYPE, std::cout, onlyCalcBeta);
        EMIT_TIMING_MESSAGE("Association took");
        return 0;
    }

    if (compareGrg) {
        START_TIMING_OPERATION();
        const grgl::GRGPtr grg2 = grgl::loadImmutableGRG(*compareGrg);
        if (!grg2) {
            return 1;
        }
        std::string disagreeReason;
        if (grgl::equivalentGRGs(theGRG, grg2, disagreeReason, /*skipRecurrent=*/true)) {
            std::cout << "COMPARE: Equivalent to " << *compareGrg << std::endl;
        } else {
            std::cout << "COMPARE: Not equivalent to " << *compareGrg << std::endl;
            std::cout << "    Reason: " << disagreeReason << std::endl;
        }
        EMIT_TIMING_MESSAGE("Comparison took ");
    }

    return 0;
}
