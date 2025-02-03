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
#ifndef GRGL_MUT_ITERATOR_H
#define GRGL_MUT_ITERATOR_H

#include "grgl/common.h"
#include "grgl/grgnode.h"
#include "util.h"
#include <utility>

// Forward declarations
namespace picovcf {
class VCFFile;
class IGDData;
} // namespace picovcf

#if BGEN_ENABLED
struct bgen_file;
struct bgen_partition;
#endif

namespace grgl {

class Mutation;

struct MutationAndSamples {
    Mutation mutation;
    NodeIDList samples;
};

/**
 * A MutationIterator provides sequential access to the _mutations_ (not variants) in an input
 * file. The difference between "mutation" and "variant" here is that the former only has one
 * alternate allele, and the latter can have many -- mutation is like a normalization.
 *
 * @param[in] binaryMutations If true, all alternate alleles will be combined into a single mutation
 *      with allele value Mutation::ALLELE_1.
 * @param[in] emitMissingData If true, an "extra" Mutation will be emitted for each variant that has
 *      missing data, with allele value Mutation::ALLELE_MISSING. Otherwise missing data will be
 *      omitted entirely (i.e., indistinguishable from the reference allele).
 * @param[in] flipRefMajor If true, flip the major allele to become the reference allele. I.e., if
 *      an alternate allele affects more than 50% of the samples then make that the reference allele
 *      and make the old reference allele an alternate.
 */
class MutationIterator {
public:
    explicit MutationIterator(FloatRange genomeRange,
                              bool binaryMutations,
                              bool emitMissingData,
                              bool flipRefMajor = false)
        : m_genomeRange(genomeRange),
          m_binaryMutations(binaryMutations),
          m_emitMissingData(emitMissingData),
          m_flipRefMajor(flipRefMajor) {}
    virtual ~MutationIterator() = default;

    MutationIterator(const MutationIterator& rhs) = delete;
    MutationIterator& operator=(const MutationIterator& rhs) = delete;
    MutationIterator(const MutationIterator&& rhs) = delete;
    MutationIterator& operator=(const MutationIterator&& rhs) = delete;

    virtual void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) = 0;

    virtual std::vector<std::string> getIndividualIds() = 0;

    virtual size_t countMutations() const = 0;

    bool next(MutationAndSamples& mutAndSamples, size_t& totalSamples);
    void reset();

    size_t numFlippedAlleles() const { return m_flippedAlleles; }

    FloatRange getBpRange() const { return m_genomeRange; }

protected:
    virtual void buffer_next(size_t& totalSamples) = 0;
    virtual void reset_specific() = 0;

    // A buffer of loaded mutations from the current position/variant.
    std::list<MutationAndSamples> m_alreadyLoaded;

    // Range to use.
    FloatRange m_genomeRange;

    // How many alleles were flipped due to m_flipRefMajor?
    size_t m_flippedAlleles{};
    // Remember how many samples we saw last "row"
    size_t m_totalSamples{};
    bool m_binaryMutations;
    bool m_emitMissingData;
    bool m_flipRefMajor;
};

/**
 * Iterate the mutations in a VCF file.
 *
 * TODO: if we keep this around, should consider adding an index file that maps genome position to
 * file position.
 */
class VCFMutationIterator : public MutationIterator {
public:
    explicit VCFMutationIterator(
        const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor);

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    std::vector<std::string> getIndividualIds() override;

protected:
    void buffer_next(size_t& totalSamples) override;
    void reset_specific() override;

private:
    std::unique_ptr<picovcf::VCFFile> m_vcf;
};

/**
 * Iterate the mutations in an IGD file.
 */
class IGDMutationIterator : public MutationIterator {
public:
    explicit IGDMutationIterator(
        const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor);

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    std::vector<std::string> getIndividualIds() override;

protected:
    void buffer_next(size_t& totalSamples) override;
    void reset_specific() override;

private:
    std::unique_ptr<picovcf::IGDData> m_igd;
    size_t m_currentVariant;
};

#if BGEN_ENABLED
/**
 * Iterate the mutations in an BGEN (v1.2 or v1.3) file.
 *
 * This assumes that the reference allele is always the 0th allele.
 */
class BGENMutationIterator : public MutationIterator {
public:
    explicit BGENMutationIterator(
        const char* filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor);
    ~BGENMutationIterator() override;

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    std::vector<std::string> getIndividualIds() override;

protected:
    void buffer_next(size_t& totalSamples) override;
    void reset_specific() override;

private:
    bgen_file* m_file;
    const bgen_partition* m_partition;
    size_t m_currentVariant{};
};
#endif

std::shared_ptr<grgl::MutationIterator> makeMutationIterator(
    const std::string& filename, FloatRange genomeRange, bool binaryMutations, bool emitMissingData, bool flipRefMajor);

} // namespace grgl

#endif /* GRGL_MUT_ITERATOR_H */
