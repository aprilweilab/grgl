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

using MutationIteratorFlags = uint64_t;
enum {
    ///< Empty flags.
    MIT_FLAG_EMPTY = 0x0,
    ///< Convert variants to be biallelic, collapsing alt alleles into a single one.
    MIT_FLAG_BINARY_MUTATIONS = 0x1,
    ///< If there is missing data, emit it as a separate Mutation, and ALWAYS BEFORE other Mutations at the same site.
    MIT_FLAG_EMIT_MISSING_DATA = 0x2,
    ///< Flip the major allele to be the reference allele, where necessary.
    MIT_FLAG_FLIP_REF_MAJOR = 0x4,
    ///< Skip mutations with no samples.
    MIT_FLAG_SKIP_EMPTY = 0x8,
    ///< Treat range as a range across _numbered variants_ instead of base pairs.
    MIT_FLAG_USE_VARIANT_RANGE = 0x10,
};

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
    explicit MutationIterator(FloatRange genomeRange, MutationIteratorFlags flags)
        : m_genomeRange(genomeRange.toIntRange()),
          m_flags(flags) {}
    virtual ~MutationIterator() = default;

    MutationIterator(const MutationIterator& rhs) = delete;
    MutationIterator& operator=(const MutationIterator& rhs) = delete;
    MutationIterator(const MutationIterator&& rhs) = delete;
    MutationIterator& operator=(const MutationIterator&& rhs) = delete;

    virtual void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) = 0;

    virtual std::vector<std::string> getIndividualIds() = 0;

    virtual size_t countMutations() const = 0;
    virtual size_t totalFileVariants() const = 0;

    bool inRange(size_t variantIndex, size_t position) const;

    /**
     * Get the next Mutation. Properties:
     * * If emitMissingData() is true, then any mutations that represent missing data will be emitted
     *   _before_ other Mutations that are at the same position (site).
     *
     * @param[out] mutAndSamples The Mutation and a list of sample identifiers that have it.
     * @param[out] totalSamples The total number of samples for the dataset.
     * @return True if a Mutation was found, false if we are at the end of the iterator.
     */
    bool next(MutationAndSamples& mutAndSamples, size_t& totalSamples);
    void reset();

    size_t numFlippedAlleles() const { return m_flippedAlleles; }

    IntRange getBpRange() const { return m_genomeRange; }

    bool binaryMutations() const { return (bool)(m_flags & MIT_FLAG_BINARY_MUTATIONS); }

    bool emitMissingData() const { return (bool)(m_flags & MIT_FLAG_EMIT_MISSING_DATA); }

    bool flipRefMajor() const { return (bool)(m_flags & MIT_FLAG_FLIP_REF_MAJOR); }

    bool skipEmpty() const { return (bool)(m_flags & MIT_FLAG_SKIP_EMPTY); }

    bool useVariantRange() const { return (bool)(m_flags & MIT_FLAG_USE_VARIANT_RANGE); }

protected:
    virtual void buffer_next(size_t& totalSamples) = 0;
    virtual void reset_specific() = 0;

    // A buffer of loaded mutations from the current position/variant.
    std::vector<MutationAndSamples> m_alreadyLoaded;

    // Range to use.
    IntRange m_genomeRange;

    // Current variant (index)
    size_t m_currentVariant{};

    // How many alleles were flipped due to m_flipRefMajor?
    size_t m_flippedAlleles{};
    // Remember how many samples we saw last "row"
    size_t m_totalSamples{};
    MutationIteratorFlags m_flags;
};

/**
 * Iterate the mutations in a VCF file.
 *
 * NOTES:
 * 1. This does not support fast random access to VCF, thus can be extremely slow for generating a GRG.
 * 2. This only supports VCFs that have all alleles on the same row, not one allele per row.
 */
class VCFMutationIterator : public MutationIterator {
public:
    explicit VCFMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags);

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    size_t totalFileVariants() const override;
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
    explicit IGDMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags);

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    size_t totalFileVariants() const override;
    std::vector<std::string> getIndividualIds() override;

protected:
    void buffer_next(size_t& totalSamples) override;
    void reset_specific() override;

private:
    std::unique_ptr<picovcf::IGDData> m_igd;
    size_t m_startVariant{};
};

#if BGEN_ENABLED
/**
 * Iterate the mutations in an BGEN (v1.2 or v1.3) file.
 *
 * This assumes that the reference allele is always the 0th allele.
 */
class BGENMutationIterator : public MutationIterator {
public:
    explicit BGENMutationIterator(const char* filename, FloatRange genomeRange, MutationIteratorFlags flags);
    ~BGENMutationIterator() override;

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override;
    size_t countMutations() const override;
    size_t totalFileVariants() const override;
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

std::shared_ptr<grgl::MutationIterator>
makeMutationIterator(const std::string& filename, FloatRange genomeRange, MutationIteratorFlags flags);

} // namespace grgl

#endif /* GRGL_MUT_ITERATOR_H */
