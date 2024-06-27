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
 * should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <args.hxx>
#include <chrono>
#include <iostream>
#include <string>

#include "calculations.h"
#include "grg_helpers.h"
#include "grgl/grg.h"
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
    args::ValueFlag<std::string> compareGrg(
        parser, "compare", "Compare the input GRG to the given GRG", {'c', "compare"});
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
        theGRG = grgl::loadImmutableGRG(*infile);
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

    if (showStats) {
        dumpStats(theGRG);
    }

    if (alleleFrequency) {
        START_TIMING_OPERATION();
        emitAlleleFrequency(theGRG, std::cout);
        EMIT_TIMING_MESSAGE("Top-down DFS (allele frequency) took");
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
