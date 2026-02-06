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
#include "grgl/mut_iterator.h"

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "picovcf.hpp"
#include "util.h"

#include <algorithm>
#include <limits>
#include <random>
#include <sstream>
#include <sys/stat.h>

#if BGEN_ENABLED
extern "C" {
#include <bgen/bgen.h>
}
#endif

template <> struct std::hash<std::pair<std::string, std::string>> {
    size_t operator()(const std::pair<std::string, std::string>& pair) const noexcept {
        return grgl::hash_combine(std::hash<std::string>{}(pair.first), std::hash<std::string>{}(pair.second));
    }
};

namespace grgl {

inline bool fileExists(const std::string& filename) {
    struct stat statBuffer;
    return stat(filename.c_str(), &statBuffer) == 0;
}

void MutationIterator::reset() {
    m_flippedAlleles = 0;
    m_alreadyLoaded.clear();
    reset_specific();
}

bool MutationIterator::inRange(size_t variantIndex, size_t position) const {
    if (useVariantRange()) {
        return m_genomeRange.contains(variantIndex);
    }
    return m_genomeRange.contains(position);
}

bool MutationIterator::next(MutationAndSamples& mutAndSamples, size_t& totalSamples) {
    totalSamples = m_totalSamples;
    if (!m_alreadyLoaded.empty()) {
        mutAndSamples.mutation = std::move(m_alreadyLoaded.back().mutation);
        mutAndSamples.samples = std::move(m_alreadyLoaded.back().samples);
        m_alreadyLoaded.pop_back();
        return true;
    }
    mutAndSamples.samples.clear();

    // If we have an alt allele with > 50% occurrence, make it the new reference allele.
    buffer_next(m_totalSamples);
    const size_t mafThreshold = std::max<size_t>(m_totalSamples / 2, 1);
    Mutation newRefAllele;
    if (flipRefMajor()) {
        for (const auto& loadedMut : m_alreadyLoaded) {
            if (loadedMut.samples.size() > mafThreshold) {
                newRefAllele = loadedMut.mutation;
            }
        }
    }
    // We found a mutation that should be switched to become the reference allele.
    if (!newRefAllele.isEmpty()) {
        assert(flipRefMajor());
        m_flippedAlleles++;
        std::vector<bool> oldAltAlleles(m_totalSamples);
        for (auto& mutAndSamples : m_alreadyLoaded) {
            for (size_t sampleId : mutAndSamples.samples) {
                oldAltAlleles.at(sampleId) = true;
            }
        }
        NodeIDList newMutSamples;
        for (size_t sampleId = 0; sampleId < oldAltAlleles.size(); sampleId++) {
            if (!oldAltAlleles[sampleId]) {
                newMutSamples.emplace_back(sampleId);
            }
        }
        for (auto& mutAndSamples : m_alreadyLoaded) {
            // Replace the old mutation with the new one (which used to be the ref allele)
            if (mutAndSamples.mutation == newRefAllele) {
                mutAndSamples.mutation =
                    Mutation(newRefAllele.getPosition(), newRefAllele.getRefAllele(), newRefAllele.getAllele());
                mutAndSamples.samples = std::move(newMutSamples);
            } else {
                mutAndSamples.mutation = Mutation(
                    mutAndSamples.mutation.getPosition(), mutAndSamples.mutation.getAllele(), newRefAllele.getAllele());
            }
        }
    }
    if (binaryMutations() && !m_alreadyLoaded.empty()) {
        auto& actualResult = m_alreadyLoaded.front();
        for (auto& mutAndSamples : m_alreadyLoaded) {
            if (actualResult.mutation == mutAndSamples.mutation) {
                actualResult.mutation = Mutation(
                    mutAndSamples.mutation.getPosition(), Mutation::ALLELE_1, mutAndSamples.mutation.getRefAllele());
            } else {
                for (auto sampleId : mutAndSamples.samples) {
                    actualResult.samples.push_back(sampleId);
                }
            }
        }
        m_alreadyLoaded.resize(1);
    }

    const bool foundMutations = !m_alreadyLoaded.empty();
    if (!m_alreadyLoaded.empty()) {
        mutAndSamples.mutation = std::move(m_alreadyLoaded.back().mutation);
        mutAndSamples.samples = std::move(m_alreadyLoaded.back().samples);
        m_alreadyLoaded.pop_back();
        totalSamples = m_totalSamples;
    }
    return foundMutations;
}

VCFMutationIterator::VCFMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags)
    : MutationIterator(genomeRange, flags),
      m_vcf(new picovcf::VCFFile(filename, PVCF_VCFFILE_CONTIG_REQUIRE_ONE)) {
    if (!m_vcf->isUsingIndex()) {
        const char* errMsg = "WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.";
        std::cerr << errMsg << std::endl;
        if (!(bool)(flags & MIT_FLAG_FORCE)) {
            throw ApiMisuseFailure(errMsg);
        }
    }
    // If we got a normalized range, we have to denormalize it before use.
    if (genomeRange.isNormalized()) {
        if (useVariantRange()) {
            const char* errMsg = "ERROR: Using variant ranges with VCF files is not supported.";
            std::cerr << errMsg << std::endl;
            throw ApiMisuseFailure(errMsg);
        }
        const auto vcfGenomeRange = m_vcf->getGenomeRange();
        m_genomeRange = genomeRange.denormalized(vcfGenomeRange.first, vcfGenomeRange.second);
    }
    getMetadataHelper();

    // Tell it to scan to the first variant.
    m_needsReset = true;
}
void VCFMutationIterator::getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) {
    ploidy = m_ploidy;
    numIndividuals = m_numIndividuals;
    isPhased = m_isPhased;
}

void VCFMutationIterator::getMetadataHelper() {
    const auto originalPosition = m_vcf->getFilePosition();
    m_numIndividuals = m_vcf->numIndividuals();
    m_ploidy = 0;
    m_isPhased = false;
    m_vcf->seekBeforeVariants();
    // We only look at the first individual of the first variant, because subsequent code in
    // this iterator ensures that we have consistent ploidy and phasedness for all data.
    if (m_vcf->nextVariant()) {
        picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        api_exc_check(variant.hasGenotypeData(), "No GT data in VCF file");
        variant.getGenotypeArray();
        m_ploidy = variant.getMaxPloidy();
        api_exc_check(variant.getPhasedness() != picovcf::PVCFP_MIXED, "Mixed phased VCF files not supported");
        m_isPhased = (variant.getPhasedness() == picovcf::PVCFP_PHASED);
    }
    m_vcf->setFilePosition(originalPosition);
}

size_t VCFMutationIterator::countMutations() const {
    size_t mutations = 0;
    const auto originalPosition = m_vcf->getFilePosition();
    release_assert(!useVariantRange());
    bool alreadyHave = m_vcf->lowerBoundPosition(m_genomeRange.start());
    while (alreadyHave || m_vcf->nextVariant()) {
        alreadyHave = false;
        picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        if (m_genomeRange.contains(variant.getPosition())) {
            mutations += variant.getAltAlleles().size();
        } else {
            break;
        }
    }
    m_vcf->setFilePosition(originalPosition);
    return mutations;
}

size_t VCFMutationIterator::totalFileVariants() const { return m_vcf->numVariants(); }

std::vector<std::string> VCFMutationIterator::getIndividualIds() { return m_vcf->getIndividualLabels(); }

void VCFMutationIterator::buffer_next(size_t& totalSamples) {
    bool foundMutations = false;
    bool alreadyHave = false;
    if (m_needsReset) {
        alreadyHave = m_vcf->lowerBoundPosition(m_genomeRange.start());
        m_needsReset = false;
    }
    while (!foundMutations && m_alreadyLoaded.empty() && (alreadyHave || m_vcf->nextVariant())) {
        alreadyHave = false;
        picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        const size_t position = variant.getPosition();
        if (!m_genomeRange.contains(position)) {
            break;
        }
        const std::string& refAllele = variant.getRefAllele();
        const std::vector<std::string> altAlleles = variant.getAltAlleles();
        std::vector<NodeIDList> altAlleleToSamples(altAlleles.size());
        NodeIDList missingSamples;

        const std::vector<picovcf::AlleleT> gtArray = variant.getGenotypeArray();
        const size_t variantPloidy = variant.getMaxPloidy();
        if (m_ploidy == 0) {
            m_ploidy = variantPloidy;
        } else if (variantPloidy != m_ploidy) {
            throw ApiMisuseFailure("VCF files with different ploidy values is not supported");
        }
        for (size_t sampleId = 0; sampleId < gtArray.size(); sampleId++) {
            const picovcf::AlleleT& allele = gtArray[sampleId];
            switch (allele) {
            case picovcf::MIXED_PLOIDY:
                throw ApiMisuseFailure("VCF files with different ploidy values is not supported");
                break;
            case picovcf::MISSING_VALUE:
                if (emitMissingData()) {
                    missingSamples.push_back(sampleId);
                }
                break;
            default:
                if (allele > 0) {
                    altAlleleToSamples.at(allele - 1).push_back(sampleId);
                }
                break;
            }
        }
        totalSamples = m_vcf->numIndividuals() / variantPloidy;

        // Convert the above information into (Mutation, NodeIDSet) pairs.
        for (size_t altAllele = 0; altAllele < altAlleleToSamples.size(); altAllele++) {
            if (binaryMutations() && !m_alreadyLoaded.empty()) {
                for (NodeID sampleId : altAlleleToSamples[altAllele]) {
                    m_alreadyLoaded.back().samples.push_back(sampleId);
                }
            } else {
                m_alreadyLoaded.push_back(
                    {Mutation(position, binaryMutations() ? Mutation::ALLELE_1 : altAlleles.at(altAllele), refAllele),
                     std::move(altAlleleToSamples[altAllele])});
            }
        }

        // If we have missing data, always emit it first (m_alreadyLoaded is a stack, so the end will emit first).
        if (emitMissingData() && !missingSamples.empty()) {
            m_alreadyLoaded.push_back(
                {Mutation(position, Mutation::ALLELE_MISSING, refAllele), std::move(missingSamples)});
        }
    }
}

void VCFMutationIterator::reset_specific() { m_needsReset = true; }

IGDMutationIterator::IGDMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags)
    : MutationIterator(genomeRange, flags),
      m_igd(new picovcf::IGDData(filename)) {
    // If we got a normalized range, we have to denormalize it before use.
    if (genomeRange.isNormalized()) {
        if (useVariantRange()) {
            m_genomeRange = genomeRange.denormalized(0, m_igd->numVariants());
        } else {
            const auto range = m_igd->getGenomeRange();
            m_genomeRange = genomeRange.denormalized(range.first, range.second);
        }
    }
    // Ick. This is a bit complex. Because of the way IGD stores variants (one
    // row per allele), we want to always start on a new position boundary. The
    // code in buffer_next() always ends on a position boundary.
    if (useVariantRange()) {
        m_startVariant = m_genomeRange.start();
        release_assert(inRange(m_startVariant, 0));
        if (m_startVariant > 0) {
            const size_t previousPosition = m_igd->getPosition(m_startVariant - 1);
            size_t position = m_igd->getPosition(m_startVariant);
            while (previousPosition == position) {
                position = m_igd->getPosition(++m_startVariant);
            }
        }
        size_t endVariant = m_genomeRange.end();
        if (endVariant < m_igd->numVariants()) {
            release_assert(endVariant > m_startVariant);
            const size_t origPosition = m_igd->getPosition(endVariant - 1);
            while (endVariant < m_igd->numVariants() && m_igd->getPosition(endVariant) == origPosition) {
                endVariant++;
            }
        }
        // Update the range.
        m_genomeRange = IntRange(m_startVariant, endVariant);
    } else {
        m_startVariant = m_igd->lowerBoundPosition(m_genomeRange.start());
        release_assert(inRange(m_startVariant, m_igd->getPosition(m_startVariant)));
    }
    m_currentVariant = m_startVariant;
}

void IGDMutationIterator::getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) {
    ploidy = m_igd->getPloidy();
    numIndividuals = m_igd->numIndividuals();
    isPhased = m_igd->isPhased();
}

size_t IGDMutationIterator::countMutations() const {
    if (useVariantRange()) {
        return m_genomeRange.span();
    }
    size_t mutations = 0;
    for (size_t i = m_startVariant; i < m_igd->numVariants(); i++) {
        if (m_genomeRange.contains(m_igd->getPosition(i))) {
            mutations++;
        } else {
            break;
        }
    }
    return mutations;
}

size_t IGDMutationIterator::totalFileVariants() const { return m_igd->numVariants(); }

std::vector<std::string> IGDMutationIterator::getIndividualIds() { return std::move(m_igd->getIndividualIds()); }

// Template for mutations that represent missing data.
using AllelePair = std::pair<std::string, std::string>;
static const AllelePair AP_FOR_MISSING = {Mutation::ALLELE_MISSING, Mutation::ALLELE_MISSING};

void IGDMutationIterator::buffer_next(size_t& totalSamples) {
    static std::random_device randDevice;
    static std::mt19937 generator(randDevice());
    std::uniform_int_distribution<size_t> copySampler(0, m_igd->getPloidy() - 1);

    release_assert(m_currentVariant >= m_startVariant);
    totalSamples = m_igd->numSamples();

    std::unordered_map<AllelePair, size_t> seenMap;
    // If we have an alt allele with > 50% occurrence, make it the new reference allele.
    while (m_currentVariant < m_igd->numVariants() && m_alreadyLoaded.empty()) {
        uint8_t numCopies = 0;
        bool isMissingDataRow = false;
        const size_t position = m_igd->getPosition(m_currentVariant, isMissingDataRow, numCopies);
        if (inRange(m_currentVariant, position)) {
            NodeIDList missingSamples;
            while (m_currentVariant < m_igd->numVariants() &&
                   position == m_igd->getPosition(m_currentVariant, isMissingDataRow, numCopies)) {
                if (isMissingDataRow) {
                    if (emitMissingData()) {
                        api_exc_check(numCopies != 1,
                                      "Unphased IGD files do not support partial missingness at the individual level");
                        const bool needsSort = !missingSamples.empty();
                        auto sampleSet = m_igd->getSamplesWithAlt(m_currentVariant);
                        for (auto sampleId : sampleSet) {
                            if (m_igd->isPhased()) {
                                missingSamples.emplace_back(sampleId);
                            } else {
                                const size_t firstCopy = sampleId * m_igd->getPloidy();
                                for (size_t j = 0; j < m_igd->getPloidy(); j++) {
                                    missingSamples.emplace_back(firstCopy + j);
                                }
                            }
                        }
                        if (needsSort) {
                            std::sort(missingSamples.begin(), missingSamples.end());
                        }
                    }
                } else {
                    AllelePair alleles = {m_igd->getRefAllele(m_currentVariant), m_igd->getAltAllele(m_currentVariant)};

                    // We collection sample sets by alleles
                    bool needsSort = false;
                    size_t loadedIndex = std::numeric_limits<size_t>::max();
                    auto findIt = seenMap.find(alleles);
                    if (findIt != seenMap.end()) {
                        loadedIndex = findIt->second;
                        needsSort = !m_alreadyLoaded[loadedIndex].samples.empty();
                    }
                    if (loadedIndex > m_alreadyLoaded.size()) {
                        loadedIndex = m_alreadyLoaded.size();
                        m_alreadyLoaded.push_back({Mutation(position, alleles.second, alleles.first), {}});
                        seenMap.emplace(std::move(alleles), loadedIndex);
                    }

                    auto sampleSet = m_igd->getSamplesWithAlt(m_currentVariant);
                    for (auto sampleId : sampleSet) {
                        if (numCopies >= 1) {
                            const size_t firstCopy = sampleId * m_igd->getPloidy();
                            for (size_t j = 0; j < numCopies; j++) {
                                m_alreadyLoaded[loadedIndex].samples.push_back(firstCopy + j);
                            }
                        } else {
                            m_alreadyLoaded[loadedIndex].samples.push_back(sampleId);
                        }
                    }
                    if (needsSort) {
                        std::sort(m_alreadyLoaded[loadedIndex].samples.begin(),
                                  m_alreadyLoaded[loadedIndex].samples.end());
                    }
                }
                m_currentVariant++;
            }
            // If we have missing data, always emit it first (m_alreadyLoaded is a stack, so the end will emit first).
            if (!missingSamples.empty()) {
                Mutation mutForMissing = {Mutation(position, Mutation::ALLELE_MISSING, Mutation::ALLELE_MISSING)};
                m_alreadyLoaded.push_back({std::move(mutForMissing), std::move(missingSamples)});
            }
        } else {
            m_currentVariant++;
            break;
        }
    }
}

void IGDMutationIterator::reset_specific() { m_currentVariant = m_startVariant; }

#if BGEN_ENABLED

#define BAD_BGEN(msgStream)                                                                                            \
    do {                                                                                                               \
        std::stringstream errMsgBG;                                                                                    \
        errMsgBG << msgStream;                                                                                         \
        throw BadInputFileFailure(errMsgBG.str().c_str());                                                             \
    } while (0)

BGENMutationIterator::BGENMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags)
    : MutationIterator(genomeRange, flags),
      m_file(bgen_file_open(filename)),
      m_partition(nullptr) {
    std::stringstream metaNameStream;
    metaNameStream << filename << ".meta";
    std::string metafilename = metaNameStream.str();

    if (nullptr == m_file) {
        BAD_BGEN("Failed to open BGEN file " << filename);
    }

    bgen_metafile* metafile = nullptr;
    if (!fileExists(metafilename)) {
        BAD_BGEN("BGEN must be indexed first. " << metafilename << " does not exist");
    }
    metafile = bgen_metafile_open(metafilename.c_str());
    release_assert(metafile != nullptr);

    const uint32_t numParts = bgen_metafile_npartitions(metafile);
    if (bgen_metafile_npartitions(metafile) != 1) {
        BAD_BGEN(metafilename << " has more than 1 partition; it is not properly indexed");
    }
    m_partition = bgen_metafile_read_partition(metafile, 0);
    release_assert(m_partition != nullptr);
    bgen_metafile_close(metafile);
    metafile = nullptr;

    const size_t numVariants = bgen_partition_nvariants(m_partition);
    if (useVariantRange()) {
        m_genomeRange = genomeRange.denormalized(0, numVariants);
    } else {
        for (size_t i = 0; i < numVariants; i++) {
            const bgen_variant* const variant = bgen_partition_get_variant(m_partition, i);
            if (variant->position < minPosition) {
                minPosition = variant->position;
            }
            if (variant->position > maxPosition) {
                maxPosition = variant->position;
            }
        }
        m_genomeRange = genomeRange.denormalized(minPosition, maxPosition);
    }
}

BGENMutationIterator::~BGENMutationIterator() {
    if (nullptr != m_file) {
        bgen_file_close(m_file);
        m_file = nullptr;
    }
    if (nullptr != m_partition) {
        bgen_partition_destroy(m_partition);
        m_partition = nullptr;
    }
}

void BGENMutationIterator::getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) {
    numIndividuals = bgen_file_nsamples(m_file);
    ploidy = 0;
    isPhased = false;
    if (bgen_partition_nvariants(m_partition) > 0) {
        const bgen_variant* const variant = bgen_partition_get_variant(m_partition, 0);
        bgen_genotype* genotype = bgen_file_open_genotype(m_file, variant->genotype_offset);
        release_assert(genotype != nullptr);
        const uint8_t minPloidy = bgen_genotype_min_ploidy(genotype);
        const uint8_t maxPloidy = bgen_genotype_max_ploidy(genotype);
        if (minPloidy != maxPloidy) {
            BAD_BGEN("BGEN file must have consistent ploidy for all samples");
        }
        ploidy = static_cast<size_t>(maxPloidy);
        isPhased = bgen_genotype_phased(genotype);

        bgen_genotype_close(genotype);
    }
    // TODO: BGEN is only supported for phased data, currently.
    release_assert(isPhased);
}

size_t BGENMutationIterator::countMutations() const {
    size_t mutations = 0;
    for (size_t i = 0; i < bgen_partition_nvariants(m_partition); i++) {
        const bgen_variant* const variant = bgen_partition_get_variant(m_partition, i);
        if (!inRange(i, variant->position)) {
            continue;
        }
        mutations += (variant->nalleles - 1);
    }
    return mutations;
}

size_t BGENMutationIterator::totalFileVariants() const { return bgen_partition_nvariants(m_partition); }

std::vector<std::string> BGENMutationIterator::getIndividualIds() {
    if (bgen_file_contain_samples(m_file)) {
        std::vector<std::string> result;
        struct bgen_samples* samples = bgen_file_read_samples(m_file);
        release_assert(samples != nullptr);
        for (size_t i = 0; i < bgen_file_nsamples(m_file); i++) {
            const bgen_string* bgStr = bgen_samples_get(samples, i);
            release_assert(bgStr != nullptr);
            result.emplace_back(bgen_string_data(bgStr));
        }
        bgen_samples_destroy(samples);
        return std::move(result);
    }
    return {};
}

void BGENMutationIterator::buffer_next(size_t& totalSamples) {
    bool foundMutations = false;
    while (m_currentVariant < bgen_partition_nvariants(m_partition) && m_alreadyLoaded.empty()) {
        const bgen_variant* const variant = bgen_partition_get_variant(m_partition, m_currentVariant);
        if (!inRange(m_currentVariant, variant->position)) {
            m_currentVariant++;
            continue;
        }

        std::vector<NodeIDList> alleleToSamples(variant->nalleles);
        NodeIDList missingSamples;

        bgen_genotype* genotype = bgen_file_open_genotype(m_file, variant->genotype_offset);
        release_assert(genotype != nullptr);
        const uint8_t minPloidy = bgen_genotype_min_ploidy(genotype);
        const uint8_t maxPloidy = bgen_genotype_max_ploidy(genotype);
        if (minPloidy != maxPloidy) {
            BAD_BGEN("BGEN file must have consistent ploidy for all samples");
        }
        if (!bgen_genotype_phased(genotype)) {
            BAD_BGEN("BGEN file must contain phased data");
        }

        const size_t numIndividuals = bgen_file_nsamples(m_file);
        const size_t probsPerIndividual = bgen_genotype_ncombs(genotype);
        const size_t numProbabilities = probsPerIndividual * numIndividuals;
        release_assert(numProbabilities >= numIndividuals);
        std::vector<double> probabilities(numProbabilities);
        if (0 != bgen_genotype_read(genotype, probabilities.data())) {
            BAD_BGEN("Failed to read genotype for variant");
        }

        const size_t numAlleles = static_cast<size_t>(bgen_genotype_nalleles(genotype));
        if (numAlleles == 0) {
            BAD_BGEN("No alleles for BGEN genotype");
        }
        const size_t numSamples = numIndividuals * maxPloidy;
        for (size_t sampleId = 0; sampleId < numSamples; sampleId++) {
            const size_t baseOffset = (numAlleles * sampleId);
            double maxProb = 0.0;
            size_t maxAllele = 0;
            double sumProb = 0.0;

            if (!bgen_genotype_missing(genotype, sampleId / maxPloidy)) {
                for (size_t k = 0; k < numAlleles; k++) {
                    const size_t idx = baseOffset + k;
                    const double prob = probabilities.at(idx);
                    sumProb += prob;
                    if (prob > maxProb) {
                        maxProb = prob;
                        maxAllele = k;
                    }
                }
                // The BGEN spec leaves out the k'th allele, but the library we use recovers it, so
                // this should be true.
                release_assert(sumProb == 1.0);
            } else if (emitMissingData()) {
                missingSamples.push_back(sampleId);
            }

            // We assume that the allele at position 0 is the reference allele.
            if (maxAllele != 0) {
                alleleToSamples.at(maxAllele).push_back(sampleId);
            }
        }

        bgen_genotype_close(genotype);
        genotype = nullptr;

        totalSamples = numSamples;

        // Convert the above information into (Mutation, NodeIDSet) pairs.
        std::string refAllele = bgen_string_data(variant->allele_ids[0]);
        for (size_t i = 0; i < alleleToSamples.size(); i++) {
            if (i == 0) {
                release_assert(alleleToSamples[i].empty());
                continue;
            }
            if (binaryMutations() && !m_alreadyLoaded.empty()) {
                for (NodeID sampleId : alleleToSamples[i]) {
                    m_alreadyLoaded.back().samples.push_back(sampleId);
                }
            } else {
                m_alreadyLoaded.push_back(
                    {Mutation(variant->position,
                              binaryMutations() ? Mutation::ALLELE_1 : bgen_string_data(variant->allele_ids[i]),
                              refAllele),
                     std::move(alleleToSamples[i])});
            }
        }

        // If we have missing data, always emit it first (m_alreadyLoaded is a stack, so the end will emit first).
        if (emitMissingData() && !missingSamples.empty()) {
            m_alreadyLoaded.push_back(
                {Mutation(variant->position, Mutation::ALLELE_MISSING, refAllele), std::move(missingSamples)});
        }

        m_currentVariant++;
    }
}

void BGENMutationIterator::reset_specific() { m_currentVariant = 0; }
#endif

std::shared_ptr<grgl::MutationIterator>
makeMutationIterator(const std::string& filename, FloatRange genomeRange, MutationIteratorFlags flags) {
    std::shared_ptr<grgl::MutationIterator> mutationIterator;
    if (ends_with(filename, ".vcf") || ends_with(filename, ".vcf.gz")) {
        mutationIterator.reset(new grgl::VCFMutationIterator(filename.c_str(), genomeRange, flags));
    } else if (ends_with(filename, ".igd")) {
        mutationIterator.reset(new grgl::IGDMutationIterator(filename.c_str(), genomeRange, flags));
#if BGEN_ENABLED
    } else if (ends_with(filename, ".bgen")) {
        mutationIterator.reset(new grgl::BGENMutationIterator(filename.c_str(), genomeRange, flags));
#endif
    } else {
        abort();
    }
    return mutationIterator;
}

} // namespace grgl
