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
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tskit.h>

#include "calculations.h"
#include "fast_build.h"
#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/map_mutations.h"
#include "grgl/mut_iterator.h"
#include "grgl/serialize.h"
#include "grgl/ts2grg.h"
#include "grgl/windowing.h"
#include "pooled_jobs.h"
#include "tskit_util.h"
#include "util.h"

#include "grgl/version.h"

inline bool supportedInputFormat(const std::string& filename) {
    return ends_with(filename, ".vcf") || ends_with(filename, ".vcf.gz") || ends_with(filename, ".igd") ||
           ends_with(filename, ".bgen");
}

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

    args::ArgumentParser parser("Genotype represenation graphs.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "infile", "The input file (must be .trees or .grg)");
    args::ValueFlag<std::string> outfile(
        parser, "outfile", "The file to write the GRG to (optional)", {'o', "outfile"});
    args::Flag showStats(parser, "show-stats", "Show statistics about the GRG", {'s', "show-stats"});
    args::Flag verbose(parser, "verbose", "Show verbose details of the GRG", {'v', "verbose"});
    args::ValueFlag<std::string> mapMutations(
        parser, "map-muts", "Map the mutations from the provided file", {'m', "map-muts"});
    args::Flag binaryMutations(parser,
                               "binary-muts",
                               "Do not store the allele with the mutation, only that a mutation occurred",
                               {'b', "binary-muts"});
    args::ValueFlag<std::string> genomeRange(
        parser,
        "genomeRange",
        "Only construct GRG for the given genome range: 'x:y' means [x, y) (x inclusive, y exclusive)",
        {'r', "range"});
    args::Flag noSimplify(parser, "no-simplify", "Compare the results to the given GRG", {'l', "no-simplify"});
    args::Flag MAFFlip(
        parser, "maf-flip", "Switch the reference allele with the major allele when they differ", {"maf-flip"});
    args::Flag showVersion(parser, "version", "Show version and exit", {"version"});
    args::Flag ignoreMissing(
        parser, "ignore-missing", "Ignore missing data (treat missing alleles as REF)", {"ignore-missing"});
    args::ValueFlag<std::string> windowedSplit(
        parser,
        "windowedSplit",
        "Split input graph into multiple output graphs. Use --outfile to specify a directory to be created "
        "that will hold the output GRG files. There are two ways to specify splitting the input graph: "
        "(1) a single number means split into GRGs where each covers the given number of BP or cM. "
        "(2) a filename containing a list of ranges, one per line, where each range is \"lower upper\" in "
        " either BP or cM (depends on whether a recombination map was specified)"
        "If the option is prefixed by a hapmap-style recombination map filename (e.g. \"filename:value\")"
        " then all positions will be interpreted as cM instead of base pair.",
        {"split"});
    args::ValueFlag<size_t> jobsArg(
        parser,
        "jobs",
        "Use this many threads for the given task. Currently only applies to the --split command",
        {'j', "jobs"});

    ///// Build tree (shape) related arguments /////
    args::ValueFlag<double> lfNoTree(
        parser,
        "lf-no-tree",
        "Ignore variants with frequency less than this threshold, for tree building. If >= 1.0, it is a count."
        " If <1.0 then it is a frequency. Default: 1.0",
        {'f', "lf-no-tree"});
    args::Flag noIndividualIds(
        parser, "no-indiv-ids", "Do not store individual string identifiers in the GRG", {"no-indiv-ids"});
    args::Flag countVariants(
        parser, "count-variants", "Just count the number of variants in the input file.", {"count-variants"});
    args::Flag noTreeMap(parser, "no-tree-map", "Don't map mutations during tree building.", {"no-tree-map"});
    args::ValueFlag<std::string> treesArg(
        parser,
        "trees",
        "How many trees to use to span the region. Either an integer between 1-K (where K is number of variants "
        " in the region) or a string value of 'optimal' (best for compression), 'faster1' (less optimal), "
        "'faster2' (even less optimal).",
        {"trees"});
    args::Flag force(parser, "force", "Ignore any warnings that normally terminate execution.", {"force"});

    ///// Tree-sequence related arguments /////
    args::Flag tsNodeTimes(parser,
                           "ts-node-times",
                           "When converting tree-seq, use node times instead of mutation times",
                           {"ts-node-times"});
    args::Flag maintainTopo(
        parser,
        "maintain-topo",
        "When converting tree-seq, maintain all topology below mutations (at the cost of a larger graph)",
        {"maintain-topo"});
    args::Flag tsComputeCoals(
        parser, "ts-coals", "When converting tree-seq, compute node individual coalescences.", {"ts-coals"});
    args::ValueFlag<std::string> populationIds(parser,
                                               "population-ids",
                                               "Format: \"filename:fieldname\". Read population ids from the given "
                                               "tab-separate file, using the given fieldname.",
                                               {"population-ids"});

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
    if (showVersion) {
        std::cout << "GRGL Version " << GRGL_MAJOR_VERSION << "." << GRGL_MINOR_VERSION << std::endl;
        return 0;
    }
    if (!infile) {
        std::cout << parser;
        return 0;
    }

    std::cout << std::fixed << std::setprecision(4);

    // Default values for parameters.
    grgl::MutationIteratorFlags itFlags = grgl::MIT_FLAG_EMPTY;
    if (binaryMutations) {
        itFlags |= grgl::MIT_FLAG_BINARY_MUTATIONS;
    }
    if (!ignoreMissing) {
        itFlags |= grgl::MIT_FLAG_EMIT_MISSING_DATA;
    }
    if (MAFFlip) {
        itFlags |= grgl::MIT_FLAG_FLIP_REF_MAJOR;
    }
    if (force) {
        itFlags |= grgl::MIT_FLAG_FORCE;
    }

    grgl::FloatRange restrictRange;
    if (genomeRange) {
        const std::string basePairSuffix = "b";
        const std::string variantsSuffix = "v";
        if (ends_with(*genomeRange, basePairSuffix) || ends_with(*genomeRange, variantsSuffix)) {
            if (ends_with(*genomeRange, variantsSuffix)) {
                itFlags |= grgl::MIT_FLAG_USE_VARIANT_RANGE;
            }
            genomeRange->resize(genomeRange->size() - 1);
        }

        auto tokens = split(*genomeRange, ':');
        release_assert(tokens.size() == 2);
        double gStart = 0.0;
        if (!parseExactDouble(tokens[0], gStart)) {
            std::cerr << "Could not parse range arg as double value." << std::endl;
            return 1;
        }
        double gEnd = 0.0;
        if (!parseExactDouble(tokens[1], gEnd)) {
            std::cerr << "Could not parse range arg as double value." << std::endl;
            return 1;
        }
        release_assert(gEnd > gStart);
        restrictRange = grgl::FloatRange(gStart, gEnd);
    }

    std::map<std::string, std::string> indivIdToPop;
    if (populationIds) {
        std::vector<std::string> parts = split(*populationIds, ':');
        if (parts.size() != 2) {
            std::cerr << "Must specify \"filename:fieldname\" for --population-ids" << std::endl;
            return 1;
        }
        indivIdToPop = loadMapFromTSV(parts[0], "sample", parts[1]);
    }

#define UNSUPPORTED_FOR_INPUT(parameter, parameterName)                                                                \
    do {                                                                                                               \
        if ((bool)(parameter)) {                                                                                       \
            std::cerr << "Parameter " << (parameterName) << " is not supported for the given input type" << std::endl; \
            return 1;                                                                                                  \
        }                                                                                                              \
    } while (0)

    grgl::GRGPtr theGRG;
    START_TIMING_OPERATION();
    if (ends_with(*infile, ".trees")) {
        UNSUPPORTED_FOR_INPUT(countVariants, "--count-variants");
        tsk_treeseq_t treeSeq;
        TSKIT_OK_OR_EXIT(tsk_treeseq_load(&treeSeq, infile->c_str(), 0), "Failed to load tree-seq");

        try {
            theGRG = grgl::convertTreeSeqToGRG(&treeSeq, binaryMutations, tsNodeTimes, maintainTopo, tsComputeCoals);
        } catch (grgl::TskitApiFailure& e) {
            std::cerr << e.what();
            return 2;
        }
    } else if (ends_with(*infile, ".grg")) {
        const bool mutableNeeded = (bool)mapMutations;
        if (mutableNeeded) {
            theGRG = grgl::loadMutableGRG(*infile);
        } else {
            theGRG = grgl::loadImmutableGRG(*infile);
        }
        if (!theGRG) {
            std::cerr << "Failed to load " << *infile << std::endl;
            return 2;
        }
        if (countVariants) {
            std::cout << theGRG->numMutations() << std::endl;
            return 0;
        }
    } else if (supportedInputFormat(*infile)) {
        if (countVariants) {
            if (ends_with(*infile, ".vcf") || ends_with(*infile, ".vcf.gz")) {
                std::cerr << "Will not count variants in VCF files (too slow)" << std::endl;
                return 1;
            }
            std::shared_ptr<grgl::MutationIterator> mutIt = makeMutationIterator(*infile, restrictRange, itFlags);
            if (!mutIt) {
                std::cerr << "Could not load input file " << *infile << std::endl;
                return 1;
            }

            std::cout << mutIt->totalFileVariants() << std::endl;
            return 0;
        }

        uint64_t buildFlags = verbose ? grgl::GBF_VERBOSE_OUTPUT : grgl::GBF_EMPTY;
        if (noIndividualIds) {
            buildFlags |= grgl::GBF_NO_INDIVIDUAL_IDS;
        }
        if (noTreeMap) {
            buildFlags |= grgl::GBF_NO_TREE_MAP;
        }
        size_t treeCount = 0;
        if (treesArg) {
            if (*treesArg == "optimal") {
            } else if (*treesArg == "faster1") {
                buildFlags |= grgl::GBF_TREES_FASTER1;
            } else if (*treesArg == "faster2") {
                buildFlags |= grgl::GBF_TREES_FASTER2;
            } else {
                uint32_t numTrees = 0;
                if (!parseExactUint32(*treesArg, numTrees)) {
                    throw std::runtime_error("Cannot parse --trees arguments");
                }
                treeCount = (size_t)numTrees;
                if (treeCount == 0) {
                    throw std::runtime_error("Must provide non-zero --trees arguments");
                }
            }
        }

        theGRG = grgl::fastGRGFromSamples(outfile ? *outfile : *infile,
                                          *infile,
                                          restrictRange,
                                          buildFlags,
                                          itFlags,
                                          treeCount,
                                          lfNoTree ? *lfNoTree : 0.0,
                                          indivIdToPop);
    } else {
        std::cerr << "Unsupported/undetected filetype for " << *infile << std::endl;
        std::cerr << "Only .trees and .grg files are supported currently." << std::endl;
        return 1;
    }
    if (verbose) {
        EMIT_TIMING_MESSAGE("Construction took ");
    }

    if (mapMutations) {
        if (lfNoTree) {
            std::cerr << "lf-no-tree not supported for mutation mapping" << std::endl;
            abort();
        }
        START_TIMING_OPERATION();
        std::shared_ptr<grgl::MutationIterator> unmappedMutations =
            makeMutationIterator(*mapMutations, restrictRange, itFlags);
        if (!unmappedMutations) {
            std::cerr << "Could not load mutations file " << *mapMutations << std::endl;
            return 1;
        }
        grgl::MutationMappingStats stats;
        stats = grgl::mapMutations(std::dynamic_pointer_cast<grgl::MutableGRG>(theGRG), *unmappedMutations, verbose);
        if (verbose) {
            EMIT_TIMING_MESSAGE("Mapping mutations took");
            std::cerr << std::endl;
            std::cerr << "=== Stats ===" << std::endl;
            stats.print(std::cerr);
            std::cerr << "Final node count: " << theGRG->numNodes() << std::endl;
            std::cerr << "Final edge count: " << theGRG->numEdges() << std::endl;
            std::cerr << std::endl;
        }
    }

    if (showStats) {
        dumpStats(theGRG);
    }

    if (windowedSplit) {
        std::stringstream splitOutPrefix;
        if (outfile) {
            splitOutPrefix << *outfile;
        } else {
            splitOutPrefix << *infile << ".split";
        }
        if (pathExists(splitOutPrefix.str())) {
            std::cerr << "Split output directory " << splitOutPrefix.str()
                      << " already exists; remove and try again (or specify a different directory to create)"
                      << std::endl;
            return 2;
        }
        makeDir(splitOutPrefix.str());
        splitOutPrefix << "/" << removeExt(basename(*infile)) << ".split";

        START_TIMING_OPERATION();
        std::string mapFile;
        const auto tokens = split(*windowedSplit, ':');
        if (tokens.size() == 2) {
            mapFile = tokens[0];
        } else if (tokens.size() != 1) {
            std::cerr << "Invalid split argument: \"" << *windowedSplit << "\"" << std::endl;
            return 2;
        }

        const size_t jobs = jobsArg ? *jobsArg : 1;
        grgl::WindowList windows;
        double perWindowAmt = 0;
        const std::string& argValue = tokens[tokens.size() - 1];
        if (!parseExactDouble(argValue, perWindowAmt)) {
            if (pathExists(argValue)) {
                std::vector<std::pair<std::string, std::string>> rangeStrings;
                try {
                    rangeStrings =
                        loadMapFromTSV<std::vector<std::pair<std::string, std::string>>>(argValue, "start", "end", ' ');
                } catch (const std::exception& e) {
                    std::cerr << "Failed to properly read file " << argValue << ": " << e.what() << std::endl;
                    return 2;
                }
                if (!mapFile.empty()) {
                    std::cerr << "TODO: Implement range file support with cM (recombination map)" << std::endl;
                    throw std::runtime_error("TODO: Implement range file support with cM (recombination map)");
                }
                for (size_t i = 0; i < rangeStrings.size(); i++) {
                    const auto& pair = rangeStrings[i];
                    bool success = true;
                    if (mapFile.empty()) {
                        uint64_t lower = 0;
                        uint64_t upper = 0;
                        success &= (bool)parseExactUint64(pair.first, lower);
                        success &= (bool)parseExactUint64(pair.second, upper);
                        windows.push_back({lower, upper});
                    } else {
                        double lower = 0.0;
                        double upper = 0.0;
                        success &= (bool)parseExactDouble(pair.first, lower);
                        success &= (bool)parseExactDouble(pair.second, upper);
                        // TODO: map these to base-pair positions via the recombination map.
                    }
                    if (!success) {
                        std::cerr << "Invalid number at line " << i + 2 << " of file " << argValue << std::endl;
                        return 2;
                    }
                }
                if (rangeStrings.empty()) {
                    std::cerr << "No range data in file " << argValue << std::endl;
                    return 2;
                }
            } else {
                std::cerr << "Bad split value: \"" << argValue << "\": invalid number and/or filename" << std::endl;
                return 2;
            }
        } else {
            const size_t overlap = 1; // TODO: add support for overlapping windows.
            release_assert(perWindowAmt > 0);
            if (mapFile.empty()) {
                windows = grgl::windowByBP(theGRG->getBPRange(), (size_t)perWindowAmt, 1);
            } else {
                windows = grgl::windowByCM(theGRG->getBPRange(), mapFile, perWindowAmt, 1);
            }
        }
        release_assert(!windows.empty());

        // Worker pool for splitting the GRG.
        class SplitJobs : public grgl::PooledJobs<grgl::Window> {
        public:
            explicit SplitJobs(grgl::GRGPtr grg, std::string filenamePrefix)
                : m_grg(std::move(grg)),
                  m_filenamePrefix(std::move(filenamePrefix)) {}

        protected:
            void processItem(grgl::Window window) override {
                grgl::NodeIDList mutations = mutationsForWindow(m_grg, window);
                std::stringstream filename;
                filename << m_filenamePrefix << "_" << window.begin << ".grg";
                saveGRGSubset(m_grg,
                              filename.str(),
                              grgl::TraversalDirection::DIRECTION_DOWN,
                              mutations,
                              {window.begin, window.end});
            }

            grgl::GRGPtr m_grg;
            std::string m_filenamePrefix;
        };

        SplitJobs workers(theGRG, splitOutPrefix.str());
        for (const auto& window : windows) {
            workers.addWork(window);
        }
        workers.doAllWork(jobs);
        EMIT_TIMING_MESSAGE("Split GRG into " << windows.size() << " parts in ");
    } else if (outfile) {
        START_TIMING_OPERATION();
        auto counts = saveGRG(theGRG, *outfile, !noSimplify);
        if (verbose) {
            std::cout << "Wrote simplified GRG with:" << std::endl;
            std::cout << "  Nodes: " << counts.first << std::endl;
            std::cout << "  Edges: " << counts.second << std::endl;
            EMIT_TIMING_MESSAGE("Wrote GRG to " << *outfile << " in ");
        }
    }

    return 0;
}
