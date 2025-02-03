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

/* Convert other formats to IGD.
 *
 * Usage:
 *  gconverter <infile name> <outfile name>
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

static grgl::NodeIDList trimIndividuals(grgl::NodeIDList sampleList,
                                        const size_t ploidy,
                                        const std::unordered_map<grgl::NodeID, grgl::NodeID>& keepIndividuals) {
    if (!keepIndividuals.empty()) {
        grgl::NodeIDList newList;
        newList.reserve(keepIndividuals.size() * ploidy);
        for (auto sampleId : sampleList) {
            const auto extra = sampleId % ploidy;
            const auto indivId = sampleId / ploidy;
            const auto& findIt = keepIndividuals.find(indivId);
            if (findIt != keepIndividuals.end()) {
                grgl::NodeID newID = (findIt->second * ploidy) + extra;
                newList.emplace_back(newID);
            }
        }
        return std::move(newList);
    }
    return std::move(sampleList);
}

static void trim(grgl::NodeIDList& sampleList, size_t newNumSamples) {
    if (0 != newNumSamples) {
        size_t trimTo = sampleList.size();
        for (size_t i = 0; i < sampleList.size(); i++) {
            if (sampleList[i] >= newNumSamples) {
                trimTo = i;
                break;
            }
        }
        sampleList.resize(trimTo);
    }
}

void addNoise(grgl::NodeIDList& sampleList,
              const size_t numSamples,
              const size_t fpThisVariant,
              const size_t fnThisVariant) {
    static std::random_device randDevice;
    static std::mt19937 generator(randDevice());
    grgl::NodeIDList dual;
    dual.reserve(numSamples - sampleList.size());
    size_t last = 0;
    for (size_t i = 0; i < sampleList.size(); i++) {
        for (size_t j = last; j < sampleList[i]; j++) {
            dual.push_back(j);
        }
        last = sampleList[i] + 1;
    }
    std::uniform_int_distribution<size_t> fnSampler(0, sampleList.size() - 1);
    const size_t newSize = (sampleList.size() > fnThisVariant) ? (sampleList.size() - fnThisVariant) : 0;
    for (size_t i = sampleList.size(); i > newSize; i--) {
        const size_t j = fnSampler(generator);
        grgl::NodeID tmp = sampleList[i - 1];
        sampleList[i - 1] = sampleList[j];
        sampleList[j] = tmp;
    }
    sampleList.resize(newSize);

    std::uniform_int_distribution<size_t> fpSampler(0, dual.size() - 1);
    const size_t additional = std::min<size_t>(fpThisVariant, dual.size());
    for (size_t i = 0; i < additional; i++) {
        const size_t j = fpSampler(generator);
        grgl::NodeID tmp = dual[i];
        dual[i] = dual[j];
        dual[j] = tmp;
    }
    for (size_t i = 0; i < additional; i++) {
        sampleList.push_back(dual[i]);
    }
    std::sort(sampleList.begin(), sampleList.end());
}

static void mutationIteratorToIGD(const std::string& inFilename,
                                  const std::string& outFilename,
                                  const grgl::FloatRange& restrictRange,
                                  const double fpPerVariant = 0,
                                  const double fnPerVariant = 0,
                                  const size_t trimToSamples = 0,
                                  std::string keepIndividualFile = "",
                                  bool removeEmpty = false,
                                  bool verbose = true) {
    constexpr size_t EMIT_EVERY = 10000;

    std::cout << "Converting " << inFilename << " to " << outFilename << std::endl;
    std::shared_ptr<grgl::MutationIterator> iterator =
        grgl::makeMutationIterator(inFilename, restrictRange, false, true, false);
    size_t ploidy = 0;
    size_t numIndividuals = 0;
    bool isPhased = false;
    iterator->getMetadata(ploidy, numIndividuals, isPhased);
    std::cout << "Loaded file with: [ploidy=" << ploidy << ", individuals=" << numIndividuals
              << ", isPhased=" << isPhased << "]" << std::endl;
    if (ploidy == 0) {
        throw grgl::BadInputFileFailure("Unexpected failure processing input file");
    }
    if (numIndividuals == 0) {
        throw grgl::BadInputFileFailure("No individuals in input file");
    }
    if (!isPhased) {
        throw grgl::BadInputFileFailure("Can only process phased input files");
    }

    const size_t numSamples = numIndividuals * ploidy;
    release_assert(trimToSamples == 0 || (trimToSamples % ploidy == 0));
    size_t effectiveSamples = (0 == trimToSamples) ? numSamples : trimToSamples;

    std::unordered_map<grgl::NodeID, grgl::NodeID> keepIndividuals;
    auto individualIds = iterator->getIndividualIds();
    if (!keepIndividualFile.empty()) {
        if (individualIds.empty()) {
            throw grgl::BadInputFileFailure("Individual filtering requires the input file to contain individual IDs");
        }
        std::map<std::string, grgl::NodeID> idToNodeId;
        for (size_t i = 0; i < individualIds.size(); i++) {
            idToNodeId.emplace(individualIds[i], i);
        }
        std::vector<std::string> newIndividualIds;
        std::ifstream filterFile(keepIndividualFile);
        std::string indivId;
        size_t count = 0;
        while (std::getline(filterFile, indivId)) {
            release_assert(!indivId.empty());
            const auto& findIt = idToNodeId.find(indivId);
            if (findIt == idToNodeId.end()) {
                std::stringstream ssErr;
                ssErr << "Could not find individual with id " << indivId;
                throw grgl::BadInputFileFailure(ssErr.str().c_str());
            }
            keepIndividuals.emplace(findIt->second, count++);
            newIndividualIds.push_back(indivId);
        }
        individualIds = std::move(newIndividualIds);
        size_t newSamples = individualIds.size() * ploidy;
        if (newSamples < effectiveSamples) {
            effectiveSamples = newSamples;
        }
    }

    std::ofstream outFile(outFilename, std::ios::binary);
    picovcf::IGDWriter writer(ploidy, effectiveSamples / ploidy, isPhased);
    writer.writeHeader(outFile, inFilename, "");

    size_t totalSamples = 0;
    grgl::MutationAndSamples mutAndSamples;
    size_t variantCount = 0;
    size_t missingCount = 0;
    double fpLeftovers = 0.0;
    double fnLeftovers = 0.0;
    while (iterator->next(mutAndSamples, totalSamples)) {
        const size_t position = (size_t)mutAndSamples.mutation.getPosition();
        const auto& refAllele = mutAndSamples.mutation.getRefAllele();
        const auto& altAllele = mutAndSamples.mutation.getAllele();
        if (altAllele == grgl::Mutation::ALLELE_MISSING) {
            missingCount += mutAndSamples.samples.size();
        }
        trim(mutAndSamples.samples, trimToSamples);
        mutAndSamples.samples = trimIndividuals(std::move(mutAndSamples.samples), ploidy, keepIndividuals);
        if (removeEmpty && mutAndSamples.samples.empty()) {
            continue;
        }
        if (fpPerVariant + fnPerVariant > 0) {
            const double fp = fpPerVariant + fpLeftovers;
            const size_t fpThisVariant = (size_t)fp;
            if ((double)fpThisVariant != fp) {
                fpLeftovers = fp - (double)fpThisVariant;
            }
            const double fn = fpPerVariant + fpLeftovers;
            const size_t fnThisVariant = (size_t)fn;
            if ((double)fnThisVariant != fn) {
                fnLeftovers = fn - (double)fnThisVariant;
            }
            addNoise(mutAndSamples.samples, effectiveSamples, fpThisVariant, fnThisVariant);
        }
        writer.writeVariantSamples(outFile, (uint64_t)position, refAllele, altAllele, mutAndSamples.samples);
        variantCount++;
        if (verbose && (variantCount % EMIT_EVERY == 0)) {
            std::cout << "Variants completed: " << variantCount << " (" << missingCount << " missing items)"
                      << std::endl;
        }
    }
    writer.writeIndex(outFile);
    writer.writeVariantInfo(outFile);
    writer.writeIndividualIds(outFile, individualIds);
    outFile.seekp(0);
    writer.writeHeader(outFile, inFilename, "");
}

int main(int argc, char* argv[]) {
    args::ArgumentParser parser("Convert between input file formats.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "infile", "The input file (.vcf or .bgen)");
    args::Positional<std::string> outfile(parser, "outfile", "The output file (.igd)");
    args::ValueFlag<size_t> trimTo(parser,
                                   "trim-to",
                                   "Trim number of samples to this value, must be even divisor of current sample count",
                                   {'t', "trim-to"});
    args::ValueFlag<double> falseNegPerVariant(
        parser, "fn-per-variant", "Number of false negatives per variant", {'n', "false-neg"});
    args::ValueFlag<double> falsePosPerVariant(
        parser, "fp-per-variant", "Number of false positives per variant", {'p', "false-pos"});
    args::ValueFlag<std::string> genomeRange(
        parser,
        "genomeRange",
        "Only convert for the given genome range: 'x:y' means [x, y) (x inclusive, y exclusive)",
        {'r', "range"});
    args::ValueFlag<std::string> keepIndivs(
        parser,
        "keepIndivs",
        "Only retain the individuals with the IDs given in this file (one ID per line).",
        {'i', "keep-indivs"});
    args::Flag removeEmpty(parser, "removeEmpty", "Remove empty alternate alleles.", {'e', "remove-empty"});
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
    if (!infile || !outfile) {
        std::cout << parser;
        return 0;
    }

    if (falseNegPerVariant) {
        std::cout << "Adding " << *falseNegPerVariant << " false negatives per variant" << std::endl;
    }
    if (falsePosPerVariant) {
        std::cout << "Adding " << *falsePosPerVariant << " false positives per variant" << std::endl;
    }
    if (trimTo) {
        std::cout << "Trimming each variant to " << *trimTo << " samples" << std::endl;
    }

    grgl::FloatRange restrictRange;
    if (genomeRange) {
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

    mutationIteratorToIGD(*infile,
                          *outfile,
                          restrictRange,
                          falseNegPerVariant ? *falseNegPerVariant : 0,
                          falsePosPerVariant ? *falsePosPerVariant : 0,
                          trimTo ? *trimTo : 0,
                          keepIndivs ? *keepIndivs : "",
                          removeEmpty,
                          /*verbose=*/true);
    return 0;
}
