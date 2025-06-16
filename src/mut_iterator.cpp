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
#include "picovcf.hpp"
#include "util.h"

#include <limits>
#include <sstream>
#include <sys/stat.h>

#if BGEN_ENABLED
extern "C" {
#include <bgen/bgen.h>
}
#endif

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
    if (m_flipRefMajor) {
        for (const auto& loadedMut : m_alreadyLoaded) {
            if (loadedMut.samples.size() > mafThreshold) {
                newRefAllele = loadedMut.mutation;
            }
        }
    }
    // We found a mutation that should be switched to become the reference allele.
    if (!newRefAllele.isEmpty()) {
        assert(m_flipRefMajor);
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
    if (m_binaryMutations && !m_alreadyLoaded.empty()) {
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

VCFMutationIterator::VCFMutationIterator(
    const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor)
    : MutationIterator(genomeRange, binaryMutations, emitMissingData, flipRefMajor),
      m_vcf(new picovcf::VCFFile(filename)) {
    // If we got a normalized range, we have to denormalize it before use.
    if (m_genomeRange.isNormalized()) {
        const auto vcfGenomeRange = m_vcf->getGenomeRange();
        m_genomeRange = m_genomeRange.denormalized(vcfGenomeRange.first, vcfGenomeRange.second);
    }
    // Scan to the first variant.
    m_vcf->seekBeforeVariants();
}

void VCFMutationIterator::getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) {
    const auto originalPosition = m_vcf->getFilePosition();
    numIndividuals = m_vcf->numIndividuals();
    ploidy = 0;
    isPhased = false;
    m_vcf->seekBeforeVariants();
    // We only look at the first individual of the first variant, because subsequent code in
    // this iterator ensures that we have consistent ploidy and phasedness for all data.
    if (m_vcf->hasNextVariant()) {
        m_vcf->nextVariant();
        const picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        if (variant.hasGenotypeData()) {
            picovcf::IndividualIteratorGT individuals = variant.getIndividualIterator();
            if (individuals.hasNext()) {
                picovcf::VariantT sample1Allele = 0;
                picovcf::VariantT sample2Allele = 0;
                isPhased = individuals.getAlleles(sample1Allele, sample2Allele);
                if (sample2Allele == picovcf::NOT_DIPLOID) {
                    ploidy = 1;
                } else {
                    ploidy = 2;
                }
            }
        }
    }
    m_vcf->setFilePosition(originalPosition);
}

size_t VCFMutationIterator::countMutations() const {
    size_t mutations = 0;
    const auto originalPosition = m_vcf->getFilePosition();
    m_vcf->seekBeforeVariants();
    while (m_vcf->hasNextVariant()) {
        m_vcf->nextVariant();
        const picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        if (m_genomeRange.contains((double)variant.getPosition())) {
            mutations += variant.getAltAlleles().size();
        }
    }
    m_vcf->setFilePosition(originalPosition);
    return mutations;
}

std::vector<std::string> VCFMutationIterator::getIndividualIds() { return m_vcf->getIndividualLabels(); }

void VCFMutationIterator::buffer_next(size_t& totalSamples) {
    bool foundMutations = false;
    while (!foundMutations && m_vcf->hasNextVariant() && m_alreadyLoaded.empty()) {
        m_vcf->nextVariant();
        const picovcf::VCFVariantView& variant = m_vcf->currentVariant();
        if (!m_genomeRange.contains((double)variant.getPosition())) {
            continue;
        }
        const std::string& refAllele = variant.getRefAllele();
        const std::vector<std::string> altAlleles = variant.getAltAlleles();
        std::vector<NodeIDList> altAlleleToSamples(altAlleles.size());
        NodeIDList missingSamples;

        // Collect the map from alternative allele index to sample sets.
        picovcf::IndividualIteratorGT individuals = variant.getIndividualIterator();
        NodeID sampleId = 0;
        while (individuals.hasNext()) {
            picovcf::VariantT sample1Allele = 0;
            picovcf::VariantT sample2Allele = 0;
            individuals.getAlleles(sample1Allele, sample2Allele);
            if (sample1Allele != 0) {
                if (sample1Allele != picovcf::MISSING_VALUE) {
                    altAlleleToSamples.at(sample1Allele - 1).push_back(sampleId);
                } else if (m_emitMissingData) {
                    missingSamples.push_back(sampleId);
                }
            }
            sampleId++;
            if (sample2Allele != picovcf::NOT_DIPLOID) {
                if (sample2Allele != 0) {
                    if (sample2Allele != picovcf::MISSING_VALUE) {
                        altAlleleToSamples.at(sample2Allele - 1).push_back(sampleId);
                    } else if (m_emitMissingData) {
                        missingSamples.push_back(sampleId);
                    }
                }
                sampleId++;
            }
        }
        totalSamples = sampleId;

        // Convert the above information into (Mutation, NodeIDSet) pairs.
        for (size_t altAllele = 0; altAllele < altAlleleToSamples.size(); altAllele++) {
            if (m_binaryMutations && !m_alreadyLoaded.empty()) {
                for (NodeID sampleId : altAlleleToSamples[altAllele]) {
                    m_alreadyLoaded.back().samples.push_back(sampleId);
                }
            } else {
                m_alreadyLoaded.push_back({Mutation(variant.getPosition(),
                                                    m_binaryMutations ? Mutation::ALLELE_1 : altAlleles.at(altAllele),
                                                    refAllele),
                                           std::move(altAlleleToSamples[altAllele])});
            }
        }
        if (m_emitMissingData && !missingSamples.empty()) {
            m_alreadyLoaded.push_back(
                {Mutation(variant.getPosition(), Mutation::ALLELE_MISSING, refAllele), std::move(missingSamples)});
        }
    }
}

void VCFMutationIterator::reset_specific() { m_vcf->seekBeforeVariants(); }

IGDMutationIterator::IGDMutationIterator(
    const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor)
    : MutationIterator(genomeRange, binaryMutations, emitMissingData, flipRefMajor),
      m_igd(new picovcf::IGDData(filename)),
      m_currentVariant(0) {
    // If we got a normalized range, we have to denormalize it before use.
    if (m_genomeRange.isNormalized()) {
        const auto range = m_igd->getGenomeRange();
        m_genomeRange = m_genomeRange.denormalized(range.first, range.second);
    }
}

void IGDMutationIterator::getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) {
    ploidy = m_igd->getPloidy();
    numIndividuals = m_igd->numIndividuals();
    isPhased = m_igd->isPhased();
}

size_t IGDMutationIterator::countMutations() const {
    size_t mutations = 0;
    for (size_t i = 0; i < m_igd->numVariants(); i++) {
        if (m_genomeRange.contains((double)m_igd->getPosition(i))) {
            mutations++;
        }
    }
    return mutations;
}

std::vector<std::string> IGDMutationIterator::getIndividualIds() { return std::move(m_igd->getIndividualIds()); }

void IGDMutationIterator::buffer_next(size_t& totalSamples) {
    totalSamples = m_igd->numSamples();

    // If we have an alt allele with > 50% occurrence, make it the new reference allele.
    while (m_currentVariant < m_igd->numVariants() && m_alreadyLoaded.empty()) {
        bool isMissingDataRow = false;
        const size_t position = m_igd->getPosition(m_currentVariant, isMissingDataRow);
        if (m_genomeRange.contains((double)position)) {
            while (m_currentVariant < m_igd->numVariants() &&
                   position == m_igd->getPosition(m_currentVariant, isMissingDataRow)) {
                if (!isMissingDataRow || m_emitMissingData) {
                    const std::string& refAllele = m_igd->getRefAllele(m_currentVariant);
                    auto allele = m_igd->getAltAllele(m_currentVariant);
                    NodeIDList samples;
                    NodeIDList samplesHomozyg;
                    auto sampleSet = m_igd->getSamplesWithAlt(m_currentVariant);
                    samples.reserve(sampleSet.size());
                    for (auto sampleId : sampleSet) {
                        samples.push_back(sampleId);
                    }
                    if (isMissingDataRow) {
                        allele = Mutation::ALLELE_MISSING;
                    }
                    m_alreadyLoaded.push_back({Mutation(position, allele, refAllele), std::move(samples)});
                }
                m_currentVariant++;
            }
        } else {
            m_currentVariant++;
        }
    }
}

void IGDMutationIterator::reset_specific() { m_currentVariant = 0; }

#if BGEN_ENABLED

#define BAD_BGEN(msgStream)                                                                                            \
    do {                                                                                                               \
        std::stringstream errMsgBG;                                                                                    \
        errMsgBG << msgStream;                                                                                         \
        throw BadInputFileFailure(errMsgBG.str().c_str());                                                             \
    } while (0)

BGENMutationIterator::BGENMutationIterator(
    const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor)
    : MutationIterator(genomeRange, binaryMutations, emitMissingData, flipRefMajor),
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

    // If we got a normalized range, we have to denormalize it before use.
    if (m_genomeRange.isNormalized()) {
        size_t minPosition = std::numeric_limits<size_t>::max();
        size_t maxPosition = 0;
        for (size_t i = 0; i < bgen_partition_nvariants(m_partition); i++) {
            const bgen_variant* const variant = bgen_partition_get_variant(m_partition, i);
            if (variant->position < minPosition) {
                minPosition = variant->position;
            }
            if (variant->position > maxPosition) {
                maxPosition = variant->position;
            }
        }
        m_genomeRange = m_genomeRange.denormalized(minPosition, maxPosition);
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
}

size_t BGENMutationIterator::countMutations() const {
    size_t mutations = 0;
    for (size_t i = 0; i < bgen_partition_nvariants(m_partition); i++) {
        const bgen_variant* const variant = bgen_partition_get_variant(m_partition, i);
        if (!m_genomeRange.contains((double)variant->position)) {
            continue;
        }
        mutations += (variant->nalleles - 1);
    }
    return mutations;
}

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
        m_currentVariant++;
        if (!m_genomeRange.contains((double)variant->position)) {
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
            } else if (m_emitMissingData) {
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
            if (m_binaryMutations && !m_alreadyLoaded.empty()) {
                for (NodeID sampleId : alleleToSamples[i]) {
                    m_alreadyLoaded.back().samples.push_back(sampleId);
                }
            } else {
                m_alreadyLoaded.push_back(
                    {Mutation(variant->position,
                              m_binaryMutations ? Mutation::ALLELE_1 : bgen_string_data(variant->allele_ids[i]),
                              refAllele),
                     std::move(alleleToSamples[i])});
            }
        }
        if (m_emitMissingData && !missingSamples.empty()) {
            m_alreadyLoaded.push_back(
                {Mutation(variant->position, Mutation::ALLELE_MISSING, refAllele), std::move(missingSamples)});
        }
    }
}

void BGENMutationIterator::reset_specific() { m_currentVariant = 0; }
#endif

std::shared_ptr<grgl::MutationIterator> makeMutationIterator(const std::string& filename,
                                                             FloatRange genomeRange,
                                                             bool binaryMutations,
                                                             bool emitMissingData,
                                                             bool flipRefMajor) {
    std::shared_ptr<grgl::MutationIterator> mutationIterator;
    if (ends_with(filename, ".vcf") || ends_with(filename, ".vcf.gz")) {
        mutationIterator.reset(new grgl::VCFMutationIterator(
            filename.c_str(), genomeRange, binaryMutations, emitMissingData, flipRefMajor));
    } else if (ends_with(filename, ".igd")) {
        mutationIterator.reset(new grgl::IGDMutationIterator(
            filename.c_str(), genomeRange, binaryMutations, emitMissingData, flipRefMajor));
#if BGEN_ENABLED
    } else if (ends_with(filename, ".bgen")) {
        mutationIterator.reset(new grgl::BGENMutationIterator(
            filename.c_str(), genomeRange, binaryMutations, emitMissingData, flipRefMajor));
#endif
    } else {
        abort();
    }
    return mutationIterator;
}

} // namespace grgl
