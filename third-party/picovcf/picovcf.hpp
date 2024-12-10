#ifndef PICOVCF_HPP
#define PICOVCF_HPP

#include <array>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fcntl.h>

#if VCF_GZ_SUPPORT
#include <zlib.h>
#endif

namespace picovcf {

// With these type sizes, we support up to 4 billion samples, and trillions of
// variants.

/** Represents the index of an individual or sample. */
using SampleT = uint32_t;
static constexpr size_t MAX_SAMPLES = std::numeric_limits<SampleT>::max() - 1;
static constexpr SampleT SAMPLE_INDEX_NOT_SET = std::numeric_limits<SampleT>::max();

/** Represents the index of a variant. */
using VariantT = uint64_t;

/** Pair of integers that define a range */
using RangePair = std::pair<VariantT, VariantT>;

/** Pair of integers that define a non-reference allele (mutation) value:
 * variant row index and alt allele index (0-based, so 0 is the first
 * non-reference allele value). */
using MutationPair = std::pair<VariantT, VariantT>;

/** When processing pairs of alleles, this is used for the second item for
 * haploids */
static constexpr VariantT NOT_DIPLOID = std::numeric_limits<VariantT>::max();

/** Represents a missing allele value (e.g., "." in VCF nomenclature) */
static constexpr VariantT MISSING_VALUE = std::numeric_limits<VariantT>::max() - 1;

static constexpr size_t INTERNAL_VALUE_NOT_SET = std::numeric_limits<size_t>::max();

/** IGD can only store up to ploidy = 8 */
static constexpr size_t MAX_PLOIDY = 8;

/**
 * Exception thrown when there is a problem reading the underlying file.
 */
class FileReadError : public std::runtime_error {
public:
    explicit FileReadError(char const* const message)
        : std::runtime_error(message) {}
};

/**
 * Exception thrown when the VCF file is not well formed according to the spec.
 */
class MalformedFile : public std::runtime_error {
public:
    explicit MalformedFile(char const* const message)
        : std::runtime_error(message) {}
};

/**
 * Exception thrown when the API is misused (bad arguments, using iterators
 * incorrectly).
 */
class ApiMisuse : public std::runtime_error {
public:
    explicit ApiMisuse(char const* const message)
        : std::runtime_error(message) {}
};

#if FUZZING
#define PICOVCF_THROW_ERROR(excType, msgOp)                                                                            \
    do {                                                                                                               \
        std::cout << msgOp;                                                                                            \
        exit(0);                                                                                                       \
    } while (0)

#define PICOVCF_ASSERT_OR_MALFORMED(condition, msgOp)                                                                  \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cout << msgOp;                                                                                        \
            exit(0);                                                                                                   \
        }                                                                                                              \
    } while (0)
#else
#define PICOVCF_THROW_ERROR(excType, msgOp)                                                                            \
    do {                                                                                                               \
        std::stringstream _ssErrMsg;                                                                                   \
        _ssErrMsg << msgOp;                                                                                            \
        throw excType(_ssErrMsg.str().c_str());                                                                        \
    } while (0)

#define PICOVCF_ASSERT_OR_MALFORMED(condition, msgOp)                                                                  \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::stringstream _ssErrMsg;                                                                               \
            _ssErrMsg << msgOp;                                                                                        \
            throw MalformedFile(_ssErrMsg.str().c_str());                                                              \
        }                                                                                                              \
    } while (0)

#endif

#define PICOVCF_RELEASE_ASSERT(condition)                                                                              \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cerr << "PICOVCF_RELEASE_ASSERT(" #condition ") failed at " << __FILE__ << ":" << __LINE__            \
                      << std::endl;                                                                                    \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)

#define PICOVCF_GOOD_OR_READ_ERROR(fstream, filename)                                                                  \
    do {                                                                                                               \
        if (!(fstream).good()) {                                                                                       \
            PICOVCF_THROW_ERROR(FileReadError, "Failure to read input file " << (filename));                           \
        }                                                                                                              \
    } while (0)

#define PICOVCF_GOOD_OR_MALFORMED_FILE(fstream)                                                                        \
    do {                                                                                                               \
        if (!(fstream).good()) {                                                                                       \
            PICOVCF_THROW_ERROR(MalformedFile,                                                                         \
                                "Unexpected file error reading input file (position " << (fstream).tellg() << ")");    \
        }                                                                                                              \
    } while (0)

#define PICOVCF_GOOD_OR_API_MISUSE(fstream)                                                                            \
    do {                                                                                                               \
        if (!(fstream).good()) {                                                                                       \
            PICOVCF_THROW_ERROR(ApiMisuse,                                                                             \
                                "File error when seeking to user-supplied location " << (fstream).tellg() << ")");     \
        }                                                                                                              \
    } while (0)

// Divide by TO and take the ceiling. E.g., for converting bits to the number of
// bytes that will be needed to store those bits (ceiling because bits may not
// be divisible by 8).
template <typename T, size_t TO> inline T picovcf_div_ceiling(T value) { return (value + (TO - 1)) / TO; }

inline void picovcf_split(const std::string& theString, const char token, std::vector<std::string>& result) {
    size_t position = 0;
    do {
        const size_t nextPosition = theString.find_first_of(token, position);
        if (nextPosition == std::string::npos) {
            result.push_back(theString.substr(position));
            position = nextPosition;
        } else {
            result.push_back(theString.substr(position, nextPosition - position));
            position = nextPosition + 1;
        }
    } while (position != std::string::npos);
}

inline std::string picovcf_getKey(const std::string& line, const size_t start) {
    size_t position = line.find_first_of('=', start);
    if (position != std::string::npos) {
        return line.substr(start, position - start);
    }
    return {};
}

inline std::string picovcf_getValue(const std::string& line, const size_t start) {
    size_t position = line.find_first_of('=', start);
    if (position != std::string::npos) {
        return line.substr(position + 1);
    }
    return {};
}

using FileOffset = std::pair<size_t, size_t>;

/**
 * PRIVATE: File reader that reads a block at a time under-the-hood.
 */
class BufferedReader {
public:
    BufferedReader(const std::string& filename, size_t bufferedByteCt)
        : m_file(std::fopen(filename.c_str(), "rb")),
          m_capacity(bufferedByteCt),
          m_position(0) {
#if HAVE_POSIX_FADVISE
        posix_fadvise(fileno(m_file), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }

    virtual ~BufferedReader() { fclose(m_file); }

    virtual void seek(const FileOffset& offset) {
        std::fseek(m_file, offset.first, SEEK_SET);
        readBuffered();
    }

    virtual FileOffset tell() { return {std::ftell(m_file) - (m_buffer.size() - m_position), 0}; }

    size_t readline(std::string& buffer) {
        if (m_buffer.empty()) {
            readBuffered();
        }
        if (m_buffer.empty()) {
            buffer.resize(0);
            return 0;
        }
        assert(m_position < m_buffer.size());
        bool found = false;
        size_t lineBytes = 0;
        while (!found) {
            const uint8_t* bufPtr = &m_buffer.data()[m_position];
            const size_t bufAvail = m_buffer.size() - m_position;
            char* p = (char*)memchr(bufPtr, '\n', bufAvail);
            size_t bytesFound;
            if (nullptr == p) {
                bytesFound = bufAvail;
                m_position = m_buffer.size();
            } else {
                assert(p <= (const char*)(m_buffer.data() + m_buffer.size()));
                assert(*p == '\n');
                bytesFound = p - (const char*)bufPtr;
                m_position += (bytesFound + 1); // +1 to skip newline
                found = true;
            }
            if (buffer.size() - lineBytes < bytesFound) {
                buffer.resize(lineBytes + bytesFound);
            }
            std::memcpy(&((char*)buffer.c_str())[lineBytes], bufPtr, bytesFound);
            lineBytes += bytesFound;
            if (m_position == m_buffer.size()) {
                if (!readBuffered()) {
                    found = true; // EOF
                }
            } else {
                assert(m_position < m_buffer.size());
                found = true;
            }
        }
        buffer.resize(lineBytes);
        return lineBytes;
    }

    virtual bool eof() {
        if (m_buffer.empty()) {
            readBuffered();
        }
        return std::feof(m_file) && (m_position == m_buffer.size());
    }

    virtual bool peek_eof() {
        if (m_buffer.empty()) {
            readBuffered();
        }
        if (m_position == m_buffer.size()) {
            return readBuffered();
        }
        return false;
    }

protected:
    // Return false if EOF
    virtual bool readBuffered() {
        m_position = 0;
        if (std::feof(m_file)) {
            m_buffer.resize(0);
            return false;
        } else if (m_buffer.size() < m_capacity) {
            m_buffer.clear();
            m_buffer.resize(m_capacity);
        }
        size_t bytesRead = std::fread(m_buffer.data(), 1, m_buffer.size(), m_file);
        if (std::ferror(m_file)) {
            PICOVCF_THROW_ERROR(FileReadError, "Failed to read input file");
        }
        if (bytesRead < m_buffer.size()) {
            m_buffer.resize(bytesRead);
        }
        return true;
    }

    FILE* m_file;
    std::vector<uint8_t> m_buffer;
    size_t m_capacity;
    size_t m_position;
};

#if VCF_GZ_SUPPORT
/**
 * PRIVATE: Buffered file reader that decompresses input via zlib.
 */
class ZBufferedReader : public BufferedReader {
public:
    ZBufferedReader(const std::string& filename, size_t bufferedByteCt)
        : BufferedReader(filename, bufferedByteCt),
          m_compressed(bufferedByteCt),
          m_lastFileRead(0) {
        constexpr int GZIP = 31;
        m_zlibStream.zalloc = Z_NULL;
        m_zlibStream.zfree = Z_NULL;
        m_zlibStream.opaque = Z_NULL;
        m_zlibStream.avail_in = 0;
        m_zlibStream.next_in = Z_NULL;
        if (Z_OK != inflateInit2(&m_zlibStream, GZIP)) {
            PICOVCF_THROW_ERROR(FileReadError, "Failed to init zlib stream.");
        }
    }

    virtual ~ZBufferedReader() { inflateEnd(&m_zlibStream); }

    void resetStream() { PICOVCF_RELEASE_ASSERT(Z_OK == inflateReset(&m_zlibStream)); }

    // Seek is not a constant operation in .gz files. We could do something like
    // the zran.c example from zlib, but I'd like to avoid that complexity unless
    // we really need it. Generally, seeking around a VCF file is going to be
    // slow, so avoid it!
    void seek(const FileOffset& offset) override {
        PICOVCF_RELEASE_ASSERT(offset.first % m_compressed.size() == 0);
        std::fseek(m_file, 0, SEEK_SET);
        resetStream();
        for (size_t i = 0; i <= offset.first; i += m_compressed.size()) {
            readBuffered();
        }
        PICOVCF_RELEASE_ASSERT(m_lastFileRead == offset.first);
        m_position = offset.second;
    }

    FileOffset tell() override {
        if (m_buffer.empty()) {
            readBuffered();
        }
        // The file offset is the block that we decompressed, and we also provide
        // the position within that block (_after_ decompression).
        PICOVCF_RELEASE_ASSERT(m_lastFileRead % m_compressed.size() == 0);
        return {m_lastFileRead, m_position};
    }

    bool eof() override {
        if (m_buffer.empty()) {
            readBuffered();
        }
        return std::feof(m_file) && (m_position == m_buffer.size());
    }

    bool peek_eof() override {
        if (m_buffer.empty()) {
            readBuffered();
        }
        if (m_position == m_buffer.size()) {
            return readBuffered();
        }
        return false;
    }

private:
    static constexpr size_t EST_COMPRESSION_FACTOR = 20;

    // Return false if EOF
    bool readBuffered() override {
        m_position = 0;
        if (std::feof(m_file)) {
            m_buffer.resize(0);
            return false;
        } else if (m_buffer.size() < m_capacity) {
            m_buffer.clear();
            m_buffer.resize(m_capacity);
        }
        const size_t bytesRead = zReadBuffered();
        if (bytesRead < m_buffer.size()) {
            m_buffer.resize(bytesRead);
        }
        return true;
    }

    // Read m_compressed.size() bytes from the file stream, and then decompress
    // that into the regular buffer. This means that the regular buffer is now
    // variable-sized, and will be as big as necessary to hold all of the
    // uncompressed data from that block.
    size_t zReadBuffered() {
        m_buffer.clear();
        m_buffer.resize(m_compressed.size() * EST_COMPRESSION_FACTOR);

        m_lastFileRead = std::ftell(m_file);
        m_zlibStream.avail_in = std::fread(m_compressed.data(), 1, m_compressed.size(), m_file);
        if (std::ferror(m_file)) {
            PICOVCF_THROW_ERROR(FileReadError, "Failed to read input file @ " << m_lastFileRead);
        }
        if (m_zlibStream.avail_in == 0)
            return 0;
        m_zlibStream.next_in = m_compressed.data();
        size_t readOffset = 0;
        // Consume all of our input data and decompress it, resizing our output
        // buffer as necessary. Note this assumes that ALL the data is compressed,
        // it does not support trailing data in the file that is uncompressed.
        do {
            PICOVCF_RELEASE_ASSERT(readOffset < m_buffer.size());
            m_zlibStream.avail_out = (m_buffer.size() - readOffset);
            m_zlibStream.next_out = (m_buffer.data() + readOffset);
            const size_t startOut = m_zlibStream.avail_out;
            int ret = inflate(&m_zlibStream, Z_NO_FLUSH);
            if (ret == Z_STREAM_END) {
                resetStream();
            } else {
                PICOVCF_ASSERT_OR_MALFORMED(ret == Z_OK, "inflate() error: " << ret);
            }
            const size_t have = startOut - m_zlibStream.avail_out;
            readOffset += have;
            if (m_zlibStream.avail_out == 0) {
                m_buffer.resize(2 * m_buffer.size());
            }
        } while (m_zlibStream.avail_in > 0);

        return readOffset;
    }

    std::vector<uint8_t> m_compressed;
    size_t m_lastFileRead;
    z_stream m_zlibStream;
};
#endif

// See if the current or next character we read will be EOF
inline bool picovcf_peek_eof(BufferedReader* instream) { return instream->eof() || instream->peek_eof(); }

inline bool isMutation(VariantT alleleIndex) { return (alleleIndex != MISSING_VALUE) && (alleleIndex != 0); }

/**
 * Iterate over individuals' genotype data only (skipping other data).
 *
 * This is a common use-case that needs to be fast and incremental.
 */
class IndividualIteratorGT {
public:
    /** The VCF diploid separator that indicates the data is phased */
    static constexpr const char PHASED_SEPARATOR = '|';

    /**
     * Are there more individuals?
     *
     * @returns true If there is another individual -- if this returns false, then
     * you cannot call getAlleles() (it will throw an exception).
     */
    bool hasNext() const { return m_currentPosition != std::string::npos; }

    /**
     * Get the allele values associated with the current individual, and (by
     * default) increment the iterator to the next individual.
     *
     * @param[out] allele1 The first (perhaps only) allele index, populated as an
     * output parameter.
     * @param[out] allele2 The second allele index, populated as an output
     * parameter. Will be set to NOT_DIPLOID as appropriate.
     * @param[in] moveNext Pass false to stay at the current individual; by
     * default it will move to the next one.
     * @returns true if the aleles are phased, false if they are not (or it is a
     * haploid).
     */
    bool getAlleles(VariantT& allele1, VariantT& allele2, bool moveNext = true) {
        if (m_currentPosition == std::string::npos) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Iterator is at end of individuals");
        }
        size_t stopAt = m_currentLine.find_first_of(":\t", m_currentPosition);
        size_t length = 0;
        if (stopAt == std::string::npos) {
            length = m_currentLine.size() - m_currentPosition;
        } else {
            length = stopAt - m_currentPosition;
        }
        // These are the two fast paths (haploid + diploid)
        bool isPhased = true;
        if (length == 1) {
            allele1 = singleCharToVariantT(m_currentLine[m_currentPosition]);
            allele2 = NOT_DIPLOID;
        } else if (length == 3) {
            allele1 = singleCharToVariantT(m_currentLine[m_currentPosition]);
            allele2 = singleCharToVariantT(m_currentLine[m_currentPosition + 2]);
            ;
            isPhased = m_currentLine[m_currentPosition + 1] == PHASED_SEPARATOR;
        } else {
            std::string alleles = m_currentLine.substr(m_currentPosition, length);
            size_t splitPos = alleles.find_first_of("|/", 0);
            PICOVCF_ASSERT_OR_MALFORMED(splitPos != std::string::npos,
                                        "Invalid allele string: " << alleles << "\nAt line: " << m_currentLine);
            char splitChar = alleles[splitPos];
            std::vector<std::string> alleleStrings;
            picovcf_split(alleles, splitChar, alleleStrings);
            PICOVCF_ASSERT_OR_MALFORMED(alleleStrings.size() == 2,
                                        "Invalid allele string: " << alleles << "\nAt line: " << m_currentLine);
            allele1 = stringToVariantT(alleleStrings[0]);
            allele2 = stringToVariantT(alleleStrings[1]);
            isPhased = (splitChar == PHASED_SEPARATOR);
        }
        if (moveNext) {
            next();
        }
        return isPhased;
    }

    /**
     * Move this iterator to the next individual.
     */
    void next() {
        if (m_currentPosition == std::string::npos) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Iterator is at end of individuals");
        }
        m_individualIndex++;
        m_currentPosition = m_currentLine.find_first_of('\t', m_currentPosition + 1);
        if (m_currentPosition != std::string::npos) {
            m_currentPosition++; // Move past separator.
        }
    }

private:
    IndividualIteratorGT(const std::string& currentLine, const size_t currentPosition)
        : m_currentLine(currentLine),
          m_currentPosition(currentPosition),
          m_individualIndex(0) {}

    static inline VariantT singleCharToVariantT(char c) {
        switch (c) {
        case '0': return 0;
        case '1': return 1;
        case '2': return 2;
        case '3': return 3;
        case '4': return 4;
        case '5': return 5;
        case '6': return 6;
        case '7': return 7;
        case '8': return 8;
        case '9': return 9;
        case '.': return MISSING_VALUE;
        }
        PICOVCF_THROW_ERROR(MalformedFile, "Invalid allele value: " << c);
    }

    static inline VariantT stringToVariantT(const std::string& str) {
        if (str == ".") {
            return MISSING_VALUE;
        }
        char* endPtr = nullptr;
        VariantT result = static_cast<uint32_t>(std::strtoull(str.c_str(), &endPtr, 10));
        if (endPtr != (str.c_str() + str.size())) {
            PICOVCF_THROW_ERROR(MalformedFile, "Invalid allele value: " << str);
        }
        return result;
    }

    // Position in the line of the first individual.
    const std::string& m_currentLine;
    size_t m_currentPosition;
    // The 0-based index of the current individual.
    SampleT m_individualIndex;

    friend class VCFVariantView;
};

/**
 * Pre-parsed information from a Variant, minus any individual data.
 */
struct VCFVariantInfo {
    std::string chromosome;
    size_t position;
    std::string identifier;
    std::string referenceAllele;
    std::vector<std::string> alternativeAlleles;
    double quality;
    std::string filter;
    std::unordered_map<std::string, std::string> information;
    std::vector<std::string> format;
};

/**
 * A subset of the VCF that has been rotated to go from sample (haplotype index)
 * to variant/mutation, excluding samples that map to the reference allele.
 *
 * VCF is ordered by variant, making it easy to go from variant/mutation to a
 * list of individuals/samples affected. This datastructure captures the inverse
 * relationship: individuals/samples mapped to the variant/mutation.
 *
 * This data structure probably doesn't make much sense for diploid unphased
 * datasets, since it is rotating per-sample (however, the resulting data could
 * still be useful, if you iterate the samples in pairs and collect all
 * mutations from both samples per pair).
 */
struct VCFRotatedWindow {
    /** The variants associated with this window */
    std::vector<VCFVariantInfo> parsedVariants;
    /** The 0-based index of the first individual in the window */
    size_t firstIndividual;
    /** The map from individual (index 0 is firstIndividual) to sets
     * of indexes into the parseVariants vector */
    std::vector<std::set<MutationPair>> sampleToMutation;
    /** The file position _after_ the last variant that was parsed - can be used
     * as an "iterator" */
    FileOffset posAfterLastVariant;
    /** Is the data diploid */
    bool isDiploid;
};

/**
 * A class that lazily interprets the data for a single variant.
 *
 * This does not "parse" the whole row associated with the variant, it locates
 * the positions of required fields and then parses those fields on demand. The
 * individual genotype data is accessed through another lazy view (iterator).
 */
class VCFVariantView {
public:
    /** Genotype format string identifier. */
    static constexpr const char* const FORMAT_GT = "GT";

    /**
     * Parse and make a copy of the non-genotype data in this variant row. This
     * can be expensive, especially without basicInfoOnly set, but does allow you
     * to capture the information from this view and not lose it when you move to
     * the next variant.
     *
     * @param[in] basicInfoOnly True by default, this only populates the
     * chromosome, position, id, ref allele, and alt allele fields. Set to false
     * to get the additional fields.
     */
    inline VCFVariantInfo parseToVariantInfo(bool basicInfoOnly = true) {
        VCFVariantInfo result = {this->getChrom(),
                                 this->getPosition(),
                                 this->getID(),
                                 this->getRefAllele(),
                                 this->getAltAlleles(),
                                 std::numeric_limits<double>::quiet_NaN()};
        if (!basicInfoOnly) {
            result.quality = this->getQuality();
            result.filter = this->getFilter();
            result.information = this->getInfo();
            result.format = this->getFormat();
        }
        return result;
    }

    /**
     * Get the chromosome identifier.
     * @returns String of the chromosome identifier.
     */
    std::string getChrom() const { return m_currentLine.substr(0, m_nonGTPositions[POS_CHROM_END]); }

    /**
     * Get the genome position.
     * @returns Integer of the genome position.
     */
    size_t getPosition() const {
        const size_t posSize = m_nonGTPositions[POS_POS_END] - m_nonGTPositions[POS_CHROM_END];
        assert(posSize > 0);
        std::string posStr = m_currentLine.substr(m_nonGTPositions[POS_CHROM_END] + 1, posSize - 1);
        char* endPtr = nullptr;
        auto result = static_cast<size_t>(strtoull(posStr.c_str(), &endPtr, 10));
        if (endPtr != (posStr.c_str() + posStr.size())) {
            PICOVCF_THROW_ERROR(MalformedFile, "Invalid position (cannot parse): " << posStr);
        }
        return result;
    }

    /**
     * Get the ID for this variant.
     * @returns String of the ID.
     */
    std::string getID() const { return stringForPosition(POS_ID_END); }

    /**
     * Get the reference allele for this variant.
     * @returns String of the reference allele.
     */
    std::string getRefAllele() const { return stringForPosition(POS_REF_END); }

    /**
     * Get the alternative alleles for this variant.
     * @returns Vector of strings for all the alternative alleles. The allele
     * index associated with each individual can be used to lookup the actual
     * allele in this vector.
     */
    std::vector<std::string> getAltAlleles() const {
        const std::string alleleStr = stringForPosition(POS_ALT_END);
        // Optimize for a common case for genotype data.
        if (alleleStr.size() == 1) {
            if (alleleStr == ".") {
                return {};
            }
            return {alleleStr};
        }
        std::vector<std::string> result;
        picovcf_split(alleleStr, ',', result);
        return result;
    }

    /**
     * Get the quality value.
     * @returns The numeric value for the quality.
     */
    double getQuality() const {
        std::string qualString = stringForPosition(POS_QUAL_END);
        char* endPtr = nullptr;
        double result = strtod(qualString.c_str(), &endPtr);
        if (endPtr != (qualString.c_str() + qualString.size())) {
            PICOVCF_THROW_ERROR(MalformedFile, "Invalid quality number (cannot parse): " << qualString);
        }
        return result;
    }

    /**
     * Get the filter value.
     * @returns The string value for the filter.
     */
    std::string getFilter() const { return stringForPosition(POS_FILTER_END); }

    /**
     * Get the info key/value pairs.
     * @returns The information as a map from key to value.
     */
    std::unordered_map<std::string, std::string> getInfo() const {
        std::string infoString = stringForPosition(POS_INFO_END);
        std::vector<std::string> infoPairs;
        picovcf_split(infoString, ';', infoPairs);
        std::unordered_map<std::string, std::string> result;
        for (auto pair : infoPairs) {
            result.emplace(picovcf_getKey(pair, 0), picovcf_getValue(pair, 0));
        }
        return result;
    }

    /**
     * @returns true if this VCF file contains genotype data.
     */
    bool hasGenotypeData() const { return INTERNAL_VALUE_NOT_SET != m_nonGTPositions[POS_FORMAT_END]; }

    /**
     * Gets the list of formats.
     *
     * Enforces that "GT" must be the first FORMAT. If the resulting vector is
     * empty then there is no FORMAT and thus there is no genotype data.
     *
     * @returns A vector the format strings.
     */
    std::vector<std::string> getFormat() const {
        if (!hasGenotypeData()) {
            return {};
        }
        std::string formatString = stringForPosition(POS_FORMAT_END);
        if (formatString.size() == 2) {
            if (formatString != FORMAT_GT) {
                PICOVCF_THROW_ERROR(MalformedFile, "The first item in FORMAT _must_ be " << FORMAT_GT);
            }
            return {formatString};
        }
        std::vector<std::string> result;
        picovcf_split(formatString, ':', result);
        if (result.empty() || result[0] != FORMAT_GT) {
            PICOVCF_THROW_ERROR(MalformedFile, "The first item in FORMAT _must_ be " << FORMAT_GT);
        }
        return result;
    }

    /**
     * Get an iterator for traversing over the individual genotype data.
     * @returns An IndividualIteratorGT for efficiently accessing the genotype
     * data.
     */
    IndividualIteratorGT getIndividualIterator() const {
        if (!hasGenotypeData()) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Cannot iterate individuals when there is no genotype data");
        }
        return IndividualIteratorGT(m_currentLine, m_nonGTPositions[POS_FORMAT_END] + 1);
    }

    // TODO we should have another (more general) iterator that can be used for
    // non-GT data.
private:
    enum {
        POS_CHROM_END = 0,
        POS_POS_END = 1,
        POS_ID_END = 2,
        POS_REF_END = 3,
        POS_ALT_END = 4,
        POS_QUAL_END = 5,
        POS_FILTER_END = 6,
        POS_INFO_END = 7,
        POS_FORMAT_END = 8,
    };

    inline std::string stringForPosition(const size_t curEnd) const {
        assert(curEnd > 0);
        const size_t prevEnd = curEnd - 1;
        const size_t strSize = m_nonGTPositions[curEnd] - m_nonGTPositions[prevEnd];
        return m_currentLine.substr(m_nonGTPositions[prevEnd] + 1, strSize - 1);
    }

    void parseNonGenotypePositions() {
        size_t position = 0;
        for (size_t i = 0; i < m_nonGTPositions.size(); i++) {
            position = m_currentLine.find_first_of("\t", position);
            if (position == std::string::npos) {
                if (i < REQUIRED_FIELDS) {
                    PICOVCF_THROW_ERROR(MalformedFile,
                                        "Invalid line (missing required fields) at " << m_currentLine.substr(0, 100));
                } else {
                    m_nonGTPositions[i] = INTERNAL_VALUE_NOT_SET;
                }
            } else {
                m_nonGTPositions[i] = position;
            }
            position++;
        }
    }

    explicit VCFVariantView(const std::string& currentLine)
        : m_currentLine(currentLine) {}

    void reset() { parseNonGenotypePositions(); }

    // FORMAT and genotype data is not required.
    static constexpr size_t REQUIRED_FIELDS = 8;

    const std::string& m_currentLine;
    std::array<size_t, 9> m_nonGTPositions;

    friend class VCFFile;
};

/**
 * A lazy parser for VCF files.
 *
 * Loads the data into memory a line at a time. Does not parse the entire line,
 * as users often only need certain pieces of information.
 */
class VCFFile {
public:
    static constexpr const char* const SUPPORTED_PREFIX = "VCFv4";
    static constexpr const char* const META_FILE_FORMAT = "fileformat";

    // 128kb for reading compressed data from the file, note this is likely to
    // result in an uncompressed buffer of a few megabytes.
    static constexpr size_t COMPRESSED_BUFFER_SIZE = 128 * 1024;
    // 1MB buffer for reading uncompressed data. Larger is faster, to a point.
    static constexpr size_t UNCOMPRESSED_BUFFER_SIZE = 1024 * 1024;

    explicit VCFFile(const std::string& filename)
        : m_variants(INTERNAL_VALUE_NOT_SET),
          m_genomeRange({INTERNAL_VALUE_NOT_SET, INTERNAL_VALUE_NOT_SET}),
          m_posVariants({INTERNAL_VALUE_NOT_SET, INTERNAL_VALUE_NOT_SET}),
          m_currentVariant(m_currentLine) {
        if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
#if VCF_GZ_SUPPORT
            m_infile = std::unique_ptr<BufferedReader>(new ZBufferedReader(filename, COMPRESSED_BUFFER_SIZE));
#else
            PICOVCF_THROW_ERROR(ApiMisuse, "picovcf was not compiled with zlib support (VCF_GZ_SUPPORT)");
#endif
        } else {
            m_infile = std::unique_ptr<BufferedReader>(new BufferedReader(filename, UNCOMPRESSED_BUFFER_SIZE));
        }
        parseHeader();

        std::string version = getMetaInfo(VCFFile::META_FILE_FORMAT);
        if (version.substr(0, 5) == SUPPORTED_PREFIX) {
            if (version.size() > 6 && !(version[6] != '0' || version[6] == '1' || version[6] == '2')) {
                PICOVCF_THROW_ERROR(MalformedFile, "Unsupported VCF version: " << version);
            }
        } else {
            PICOVCF_THROW_ERROR(MalformedFile, "Unsupported VCF version: " << version);
        }
    }

    /**
     * Compute the number of variants in the file: expensive!
     * This is not a constant time operation, it involves scanning the entire
     * file.
     *
     * @return The number of variants.
     */
    size_t numVariants() {
        if (m_variants == INTERNAL_VALUE_NOT_SET) {
            scanVariants();
        }
        return m_variants;
    }

    /**
     * Compute the range of positions for variants in the file.
     * This is not a constant time operation, it involves scanning the entire
     * file.
     *
     * @return A pair of the minimum and maximum variant positions present in the
     * file.
     */
    RangePair getGenomeRange() {
        if (m_genomeRange.first == INTERNAL_VALUE_NOT_SET) {
            scanVariants();
        }
        return m_genomeRange;
    }

    /**
     * Get the number of individuals with labels in the VCF file.
     * @return number of individuals.
     */
    size_t numIndividuals() { return m_individualLabels.size(); }

    /**
     * Get a list of the labels for the individuals in the VCF file.
     * @return vector of strings, where the 0th is the 0th individuals label, etc.
     */
    std::vector<std::string>& getIndividualLabels() { return m_individualLabels; }

    /**
     * Get a metadata value from the VCF header rows.
     * @param[in] key The metadata key name.
     * @return the string associated with the given key, or empty string.
     */
    std::string getMetaInfo(const char* const key) const {
        const auto metaIt = m_metaInformation.find(key);
        if (metaIt != m_metaInformation.end()) {
            return metaIt->second;
        }
        return {};
    }

    /**
     * Get an opaque handle describing the current file position of the parser.
     * @return Current FileOffset.
     */
    FileOffset getFilePosition() { return m_infile->tell(); }

    /**
     * Use an opaque handle to return to a previously-recorded file position.
     * @param[in] position A FileOffset saved via getFilePosition().
     */
    void setFilePosition(const FileOffset& position) { m_infile->seek(position); }

    /**
     * Change the parser position to be immediately _before_ the first variant.
     */
    void seekBeforeVariants() {
        PICOVCF_ASSERT_OR_MALFORMED(m_posVariants.first != INTERNAL_VALUE_NOT_SET, "File has no variant data");
        setFilePosition(m_posVariants);
    }

    /**
     * Is there a variant at the current file position?
     * @returns true if calling nextVariant() will place us at a valid variant.
     */
    bool hasNextVariant() { return !picovcf_peek_eof(m_infile.get()); }

    /**
     * Read the variant at the current file position and move the file position to
     * the following variant.
     */
    void nextVariant() {
        if (picovcf_peek_eof(m_infile.get())) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Tried to move to next variant when there aren't any");
        }
        m_infile->readline(m_currentLine);
        m_currentVariant.reset();
    }

    /**
     * Get a parseable view of the variant that we last encountered with
     * nextVariant().
     * @return A VCFVariantView that can be queried for variant information.
     */
    VCFVariantView& currentVariant() { return m_currentVariant; }

    /**
     * Get a VCFRotatedWindow for the given rectange defined by individualRange
     * and genomeRange.
     *
     * The most efficient way to use this function is from left-to-right
     * (individuals) first, and then top-to-bottom (variants). E.g. if you have a
     * VCF laid out like this: xyzw XYZW Then ask for the windows in this order:
     * x,y,z,w,X,Y,Z,W. This allows use of posAfterLastVariant for efficiency.
     *
     * Assumptions:
     * 1. The VCF rows are ordered by genome position (ascending).
     * 2. The genotype data is uniform in ploidy and phased-ness.
     *
     * @param[in] individualRange The range [start, end) of individual indexes to
     * include. For example, [0, 10) will only include individuals 0-9.
     * @param[in] genomeRange The range [start, end) of genome positions to
     * include. For example, [100, 1000) will include any variants that have
     * position 100-999.
     */
    void getRotatedWindow(RangePair individualRange,
                          RangePair genomeRange,
                          VCFRotatedWindow& result,
                          FileOffset posAfterLastVariant = {0, 0}) {
        const size_t totalIndividuals = this->numIndividuals();
        if (individualRange.second <= individualRange.first || individualRange.first >= totalIndividuals) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Invalid individualRange argument");
        }
        if (genomeRange.second <= genomeRange.first) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Invalid genomeRange argument");
        }
        if (posAfterLastVariant != FileOffset(0, 0)) {
            this->setFilePosition(posAfterLastVariant);
        } else {
            this->seekBeforeVariants();
        }
        auto prevFilePosition = this->getFilePosition();
        const size_t numIndividuals =
            std::min<size_t>(totalIndividuals, individualRange.second) - individualRange.first;
        result.firstIndividual = individualRange.first;
        result.sampleToMutation.clear();
        result.parsedVariants.clear();
        result.isDiploid = false;
        while (this->hasNextVariant()) {
            prevFilePosition = this->getFilePosition();
            this->nextVariant();
            VCFVariantView& variant = this->currentVariant();
            const size_t currentVariantIndex = result.parsedVariants.size();
            const double position = variant.getPosition();
            if (position >= genomeRange.first && position < genomeRange.second) {
                IndividualIteratorGT iterator = variant.getIndividualIterator();
                // Scan for the start position.
                for (size_t i = 0; i < individualRange.first && iterator.hasNext(); i++, iterator.next())
                    ;

                for (size_t i = 0; i < (individualRange.second - individualRange.first) && iterator.hasNext();
                     i++, iterator.next()) {
                    VariantT allele1 = 0, allele2 = 0;
                    iterator.getAlleles(allele1, allele2, /*moveNext=*/false);
                    if (result.sampleToMutation.empty()) {
                        result.isDiploid = (allele2 != NOT_DIPLOID);
                        result.sampleToMutation.resize(numIndividuals * (result.isDiploid ? 2 : 1));
                    }
                    const size_t multiplier = (result.isDiploid ? 2 : 1);
                    if (isMutation(allele1)) {
                        result.sampleToMutation.at(multiplier * i).emplace(currentVariantIndex, allele1);
                    }
                    // TODO should probably also check for consistent phasing.
                    if ((allele2 != NOT_DIPLOID) != result.isDiploid) {
                        PICOVCF_THROW_ERROR(ApiMisuse, "getRotatedWindow requires consistent ploidy");
                    }
                    if (isMutation(allele2)) {
                        result.sampleToMutation.at((multiplier * i) + 1).emplace(currentVariantIndex, allele2);
                    }
                }
                const VCFVariantInfo variantParsed = variant.parseToVariantInfo();
                result.parsedVariants.push_back(std::move(variantParsed));
            } else if (position >= genomeRange.second || !result.parsedVariants.empty()) {
                break;
            }
        }
        result.posAfterLastVariant = prevFilePosition;
    }

private:
    void parseHeader() {
        FileOffset prevPosition = {0, 0};
        std::string lineBuffer;
        bool sawHeader = false;
        while (!picovcf_peek_eof(m_infile.get())) {
            m_infile->readline(lineBuffer);
            if (lineBuffer.size() >= 1 && lineBuffer[0] == '#') {
                if (lineBuffer.size() >= 2 && lineBuffer[1] == '#') {
                    std::string key = picovcf_getKey(lineBuffer, 2);
                    std::string value = picovcf_getValue(lineBuffer, 2);
                    if (key.empty() || value.empty()) {
                        PICOVCF_THROW_ERROR(MalformedFile, "Invalid key/value pair at " << m_infile->tell().first);
                    }
                    m_metaInformation.emplace(key, value);
                } else {
                    std::vector<std::string> headerColumns;
                    picovcf_split(lineBuffer, '\t', headerColumns);
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns.size() >= 8, "Header line too short");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[0] == "#CHROM",
                                                "CHROM missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[1] == "POS", "POS missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[2] == "ID", "ID missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[3] == "REF", "REF missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[4] == "ALT", "ALT missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[5] == "QUAL", "QUAL missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[6] == "FILTER",
                                                "FILTER missing/out of order in header line");
                    PICOVCF_ASSERT_OR_MALFORMED(headerColumns[7] == "INFO", "INFO missing/out of order in header line");
                    if (headerColumns.size() > 8) {
                        PICOVCF_ASSERT_OR_MALFORMED(headerColumns[8] == "FORMAT",
                                                    "FORMAT missing/out of order in header line");
                        m_individualLabels.reserve(headerColumns.size() - 9);
                        for (size_t i = 9; i < headerColumns.size(); i++) {
                            m_individualLabels.push_back(std::move(headerColumns[i]));
                        }
                    }
                    sawHeader = true;
                }
            } else {
                m_posVariants = prevPosition;
                break;
            }
            prevPosition = m_infile->tell();
        }
        PICOVCF_ASSERT_OR_MALFORMED(sawHeader, "No header line (column names) found.");
    }

    void scanVariants() {
        m_variants = 0;
        const auto originalPosition = getFilePosition();
        this->seekBeforeVariants();
        m_genomeRange.second = 0;
        while (this->hasNextVariant()) {
            this->nextVariant();
            VCFVariantView variant = this->currentVariant();
            const double position = variant.getPosition();
            if (m_genomeRange.first == INTERNAL_VALUE_NOT_SET) {
                m_genomeRange.first = position;
            }
            PICOVCF_ASSERT_OR_MALFORMED(m_genomeRange.second <= position,
                                        "VCF rows must be in ascending genome position");
            m_genomeRange.second = position;
            m_variants++;
        }
        setFilePosition(originalPosition);
    }

    std::unique_ptr<BufferedReader> m_infile;
    // The number of variants represented
    size_t m_variants;
    // The min/max genome range
    RangePair m_genomeRange;
    // Starting position (offset) in the file for the variants information
    FileOffset m_posVariants;
    // Metadata dictionary
    std::unordered_map<std::string, std::string> m_metaInformation;
    // Labels for individuals
    std::vector<std::string> m_individualLabels;

    // These three members represent the current variant we're looking at.
    VCFVariantView m_currentVariant;
    std::string m_currentLine;
};

template <typename T> static inline T readScalar(std::istream& inStream) {
    T simpleValue = 0;
    inStream.read(reinterpret_cast<char*>(&simpleValue), sizeof(simpleValue));
    return simpleValue;
}

static inline size_t readStringLen(const uint64_t version, std::istream& inStream) {
    if (version == 3) {
        return readScalar<uint64_t>(inStream);
    }
    return readScalar<uint32_t>(inStream);
}

static inline std::string readString(const uint64_t version, std::istream& inStream) {
    const auto strLength = readStringLen(version, inStream);
    PICOVCF_GOOD_OR_MALFORMED_FILE(inStream);
    std::string strValue;
    strValue.resize(strLength);
    inStream.read(const_cast<char*>(strValue.c_str()), strLength);
    return std::move(strValue);
}

/** Vector of sample indexes (IDs) */
using IGDSampleList = std::vector<SampleT>;

// NOTE: in order to use this for larger than byte at a time, some more math is
// needed to adjust the sample IDs. This is more robust to uniformly distributed
// bits than other methods, in that it skips 0 bytes and skips 0-prefixes.
template <typename T> inline void samplesForBV(const T* buffer, SampleT index, IGDSampleList& result) {
    constexpr SampleT bits = sizeof(T) * 8;
    constexpr SampleT mask = 0x1 << (bits - 1);
    static_assert(sizeof(T) == 1, "More math needed to support larger types");
    T value = buffer[index];
    SampleT bit = 0;
    const SampleT sampleOffset = index * bits;
    while (value != 0) {
        const bool isSet = (bool)(value & mask);
        if (isSet) {
            const SampleT sampleId = sampleOffset + bit;
            result.push_back(sampleId);
        }
        value <<= 1U;
        bit++;
    }
}

// INTERNAL helper.
inline IGDSampleList getSamplesWithAlt(const uint8_t* buffer, const SampleT numSamples) {
    IGDSampleList result;
    for (SampleT idx = 0; idx < picovcf_div_ceiling<SampleT, 8>(numSamples); idx++) {
        samplesForBV<uint8_t>((uint8_t*)buffer, idx, result);
    }
    return std::move(result);
}

// INTERNAL CLASS
// Compact representation of allele values.
class IGDAllele {
public:
    static constexpr uint32_t IS_LONG = 0x80000000;

    IGDAllele(const std::string& stringVal, size_t numLongAlleles)
        : m_value() {
        if (stringVal.size() <= sizeof(m_value)) {
            for (size_t i = 0; i < stringVal.size(); i++) {
                ((char*)&m_value)[i] = stringVal[i];
            }
            PICOVCF_ASSERT_OR_MALFORMED(!(m_value & IS_LONG), "Invalid allele string " << stringVal);
            assert(stringVal.size() == sizeof(m_value) || ((char*)&m_value)[sizeof(m_value) - 1] == 0);
        } else {
            PICOVCF_ASSERT_OR_MALFORMED(
                !(numLongAlleles & IS_LONG || numLongAlleles > std::numeric_limits<uint32_t>::max()),
                "Too many long alleles: " << numLongAlleles);
            m_value = ((uint32_t)numLongAlleles) | IS_LONG;
        }
    }

    std::string getString() const {
        std::string result;
        for (size_t i = 0; i < sizeof(m_value); i++) {
            char c = ((char*)&m_value)[i];
            if (0 == c) {
                break;
            }
            result.push_back(c);
        }
        return std::move(result);
    }

    uint32_t getLongIndex() const { return m_value & (~IS_LONG); }

    bool isLong() const { return m_value & IS_LONG; }

private:
    uint32_t m_value;
};
static_assert(sizeof(IGDAllele) == 4, "IGDAllele size changed");

/**
 * Indexable individual genotype data.
 *
 * This class is for reading data from an IGD file. See IGDWriter for generating
 * IGD output.
 */
class IGDData {
public:
    /** Uniquely identifies this file as an IGD */
    static constexpr uint64_t IGD_MAGIC = 0x3a0c6fd7945a3481;
    /** Flag that indicates the genotype data is phased */
    static constexpr uint64_t IGD_PHASED = 0x1;
    /** The IGD file format version that this library writes */
    static constexpr uint64_t CURRENT_IGD_VERSION = 4;
    /** If fewer than NumSamples/DEFAULT_SPARSE_THRESHOLD samples have a
       particular variant, we will store it sparsely. */
    static constexpr uint32_t DEFAULT_SPARSE_THRESHOLD = 32;

    /* See IGD.FORMAT.md for detailed file format. The layout is:
     * - FixedHeader (128 bytes)
     * - Source (string: 8 byte length, followed by that many bytes)
     * - Description (string: 8 byte length, followed by that many bytes)
     * - Genotype data rows, >=1 per variant.
     *    - multi-allelic variants are expanded into multiple rows
     *    - missing data is expanded into its own row
     *    - each row may be either sparse (list of sample indices) or dense (bit
     * vector)
     * - At arbitrary locations in the file, as defined by the header:
     *    - An index of all the variants, with the genomic and file position, and
     * flags that indicate whether it is stored sparsely or not.
     *    - A list of the variant information, which contains the allele strings
     * themselves.
     *    - A list of identifiers for the individuals (samples) in the dataset
     * (optional).
     *
     * There can be many variants for each position. Each variant is a single
     * alternative. So conversion from VCF would take a single variant with N
     * alternate alleles and create N variants, each with the same reference
     * allele.
     */
#pragma pack(push, 1)
    struct FixedHeader {
        uint64_t magic;                // Magic identifier to say this is an IGD file.
        uint64_t version;              // IGD file format version.
        uint32_t ploidy;               // Ploidy of every individual.
        uint32_t sparseThreshold;      // Number of samples below which we store a
                                       // variant sparsely.
        uint64_t numVariants;          // Total number of variants.
        uint32_t numIndividuals;       // Total number of individuals.
        uint32_t unusedCenter;         // Unused
        uint64_t flags;                // Flags that indicate properties of the dataset.
        uint64_t filePosIndex;         // Byte offset in the file where the variant rows start.
        uint64_t filePosVariants;      // Byte offset in the file where the variant
                                       // allele information start.
        uint64_t filePosIndividualIds; // Byte offset in the file where the
                                       // identifiers for individuals start.
        uint64_t filePosVariantIds;    // Byte offset in the file where the identifiers
                                       // for variants start.
        uint64_t unused[6];
    };
    static_assert(sizeof(FixedHeader) == 128, "Header size is fixed");

    // Bitfields are undefined in the C standard how they are laid out in memory,
    // so we can't rely on them for a consistent serialization format. Instead we
    // "steal" the upper byte of the bpPosition to store whether or not the row at
    // the given position is sparse.
    static constexpr uint64_t BP_POS_FLAGS_MASK = 0xFF00000000000000;
    static constexpr uint64_t BP_POS_FLAGS_SPARSE = 0x0100000000000000;
    static constexpr uint64_t BP_POS_FLAGS_IS_MISSING = 0x0200000000000000;

    // We steal the next byte for indicating how many copies of the alternate allele
    // are represented by the current variant. For example, a value of 2 for alt allele
    // 'A' and ref allele 'T' means that the diploid has AA. A value of 1 means AT or
    // TA (order is irrelevant), and the absence of the sample in any list means TT.
    static constexpr uint64_t BP_POS_COPY_SHIFT = 48;
    static constexpr uint64_t BP_POS_COPY_MASK = ((uint64_t)0xFF << BP_POS_COPY_SHIFT);

    // The mask to get only the position value.
    static constexpr uint64_t BP_POS_ONLY_MASK = ~(BP_POS_FLAGS_MASK | BP_POS_COPY_MASK);

    // The maximum value of a position, (2^48)-1
    static constexpr uint64_t MAX_BP_POSITION = 0x0000FFFFFFFFFFFF;

    static_assert(MAX_BP_POSITION == BP_POS_ONLY_MASK, "Unexpected BP_POS mask");
    static_assert((BP_POS_FLAGS_MASK & BP_POS_COPY_MASK) == 0, "Overlapping masks");

    struct IndexEntry {
        uint64_t bpPosition;     // The base-pair position of the variant.
        uint64_t filePosDataRow; // Byte offset in the file where the corresponding
                                 // row is (samples containing variant).
    };
    static_assert(sizeof(IndexEntry) == 16, "IndexEntry is fixed size");
#pragma pack(pop)

    /**
     * Loads the header of an IGD file and prepares for accessing the variant and
     * genotype data. This class provides access on a per-variant level, it does
     * not load the entire data into memory.
     *
     * @param[in] filename The file path for an IGD file.
     */
    explicit IGDData(const std::string& filename)
        : m_infile(filename, std::ios::binary),
          m_header({0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    std::numeric_limits<uint64_t>::max(),
                    std::numeric_limits<uint64_t>::max(),
                    std::numeric_limits<uint64_t>::max(),
                    std::numeric_limits<uint64_t>::max()}),
          m_beforeFirstVariant(0) {
        readAndCheckHeader();
    }

    /**
     * The number of variants in this file. Note that variants are binary in IGD,
     * so this is not necessarily the same as the number of polymorphic sites (i.e.,
     * rows in a VCF file). Also for unphased data, each number of copies that is
     * stored is considered a variant, so e.g., for diploid datasets Aa and AA are
     * considered separate variants.
     *
     * @return Number of variants (rows of genotype data).
     */
    uint64_t numVariants() const { return m_header.numVariants; }

    /**
     * The number of individuals represented in the genotype data.
     * @return Number of individuals.
     */
    SampleT numIndividuals() const { return m_header.numIndividuals; }

    /**
     * IGD files have a fixed ploidy, so this is either getPloidy()*numIndividuals()
     * (for phased data) or the same as numIndividuals() (for unphased).
     *
     * @return Number of samples. Every sample index will be >= 0 and < numSamples().
     */
    SampleT numSamples() const {
        if (m_header.flags & IGD_PHASED) {
            return m_header.numIndividuals * m_header.ploidy;
        }
        return numIndividuals();
    }

    /**
     * A string describing where this file came from.
     * @return A string description.
     */
    std::string getSource() const { return m_source; }

    /**
     * Free-form description of the file contents.
     * @return A string description.
     */
    std::string getDescription() const { return m_description; }

    /**
     * The ploidy of _all_ the data contained in the file.
     * @return ploidy value >= 1.
     */
    uint64_t getPloidy() const { return m_header.ploidy; }

    /**
     * The phasedness of all the data in the IGD file.
     * @return true if phased.
     */
    bool isPhased() const { return m_header.flags & IGD_PHASED; }

    /**
     * Get the range of positions for variants in the file.
     * @return A pair of the minimum and maximum variant positions present in the
     * file.
     */
    RangePair getGenomeRange() {
        bool _ignore = false;
        VariantT firstPosition = this->getPosition(0, _ignore);
        VariantT lastPosition = this->getPosition(this->numVariants() - 1, _ignore);
        return {firstPosition, lastPosition};
    }

    /**
     * Get the position of the given variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @param[out] isMissing Will be set to true if this variant represents
     * missing data.
     * @param[out] numCopies Will be set to the number of copies of the alternate
     * allele that this variant represents.
     * @return The position of the variant on the genome.
     */
    uint64_t getPosition(VariantT variantIndex, bool& isMissing, uint8_t& numCopies) {
        const size_t offset = getVariantIndexOffset(variantIndex);
        m_infile.seekg(offset);
        PICOVCF_GOOD_OR_API_MISUSE(m_infile);
        const uint64_t position = readScalar<uint64_t>(m_infile);
        isMissing = (bool)(position & BP_POS_FLAGS_IS_MISSING);
        numCopies = (uint8_t)(position >> BP_POS_COPY_SHIFT);
        return position & BP_POS_ONLY_MASK;
    }

    /**
     * Get the position of the given variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @param[out] isMissing Will be set to true if this variant represents
     * missing data.
     * @return The position of the variant on the genome.
     */
    uint64_t getPosition(VariantT variantIndex, bool& isMissing) {
        const size_t offset = getVariantIndexOffset(variantIndex);
        m_infile.seekg(offset);
        PICOVCF_GOOD_OR_API_MISUSE(m_infile);
        const uint64_t position = readScalar<uint64_t>(m_infile);
        isMissing = (bool)(position & BP_POS_FLAGS_IS_MISSING);
        return position & BP_POS_ONLY_MASK;
    }

    /**
     * Get the position of the given variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @return The position of the variant on the genome.
     */
    uint64_t getPosition(VariantT variantIndex) {
        bool _ignore = false;
        return getPosition(variantIndex, _ignore);
    }

    /**
     * Get the single alternative allele for the given variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @return The string representing the alternative allele of the variant.
     */
    std::string getAltAllele(VariantT variantIndex) {
        // TODO: do we really want to read these all at once?
        if (m_alternateAlleles.empty()) {
            readAllAlleleInfo();
        }
        IGDAllele allele = m_alternateAlleles.at(variantIndex);
        if (allele.isLong()) {
            return m_longAlleles.at(allele.getLongIndex());
        }
        return allele.getString();
    }

    /**
     * Get the reference allele for the given variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @return The string representing the reference allele of the variant.
     */
    std::string getRefAllele(VariantT variantIndex) {
        // TODO: do we really want to read these all at once?
        if (m_referenceAlleles.empty()) {
            readAllAlleleInfo();
        }
        IGDAllele allele = m_referenceAlleles.at(variantIndex);
        if (allele.isLong()) {
            return m_longAlleles.at(allele.getLongIndex());
        }
        return allele.getString();
    }

    /**
     * Get the list of samples that have the alternate allele for the given
     * variant.
     * @param[in] variantIndex The 0-based index of the variant, i.e. the row
     * number.
     * @return A list (std::vector) of the sample indexes. Order is based on
     * individual, and then ploidy within the individual. E.g., the 0th diploid
     * individual will have sample indexes 0 and 1, the 1st will have 2 and 3,
     * etc.
     */
    IGDSampleList getSamplesWithAlt(VariantT variantIndex) {
        const size_t base = getVariantIndexOffset(variantIndex);
        m_infile.seekg(base);
        PICOVCF_GOOD_OR_API_MISUSE(m_infile);
        IndexEntry index;
        m_infile.read(reinterpret_cast<char*>(&index), sizeof(index));
        const bool isSparse = (bool)(index.bpPosition & BP_POS_FLAGS_SPARSE);

        m_infile.seekg(index.filePosDataRow);
        PICOVCF_GOOD_OR_API_MISUSE(m_infile);
        if (isSparse) {
            const SampleT numSamples = readScalar<SampleT>(m_infile);
            PICOVCF_GOOD_OR_MALFORMED_FILE(m_infile);
            IGDSampleList sampleList(numSamples);
            for (SampleT j = 0; j < numSamples; j++) {
                sampleList[j] = readScalar<SampleT>(m_infile);
            }
            return std::move(sampleList);
        } else {
            const bool phased = (bool)(m_header.flags & IGD_PHASED);
            const SampleT numSamples = phased ? (m_header.ploidy * m_header.numIndividuals) : m_header.numIndividuals;
            const SampleT readAmount = picovcf_div_ceiling<SampleT, 8>(numSamples);
            PICOVCF_RELEASE_ASSERT(readAmount > 0);
            std::unique_ptr<uint8_t> buffer(new uint8_t[readAmount]);
            if (buffer) {
                m_infile.read(reinterpret_cast<char*>(buffer.get()), readAmount);
                PICOVCF_GOOD_OR_API_MISUSE(m_infile);
                return ::picovcf::getSamplesWithAlt((const uint8_t*)buffer.get(), numSamples);
            }
        }
        return {};
    }

    /**
     * Read the (optional) list of individual identifiers from the file.
     *
     * @return A newly created std::vector<std::string> object. Individual 0 has
     * its label at position 0, and the last individual has its label at position
     * (numIndividuals()-1).
     */
    std::vector<std::string> getIndividualIds() {
        std::vector<std::string> result;
        if (0 != m_header.filePosIndividualIds) {
            m_infile.seekg(m_header.filePosIndividualIds);
            const size_t numLabels = readScalar<uint64_t>(m_infile);
            PICOVCF_ASSERT_OR_MALFORMED(numLabels == m_header.numIndividuals, "Invalid number of individual ids");
            for (size_t i = 0; i < numLabels; i++) {
                result.emplace_back(std::move(readString(m_header.version, m_infile)));
            }
        }
        return std::move(result);
    }

    /**
     * Read the (optional) list of variant identifiers from the file.
     *
     * @return A newly created std::vector<std::string> object. Variant 0 has its
     * label at position 0, and the last variant has its identifier at position
     * (numVariants()-1).
     */
    std::vector<std::string> getVariantIds() {
        std::vector<std::string> result;
        if (0 != m_header.filePosVariantIds) {
            m_infile.seekg(m_header.filePosVariantIds);
            const size_t numLabels = readScalar<uint64_t>(m_infile);
            PICOVCF_ASSERT_OR_MALFORMED(numLabels == m_header.numVariants, "Invalid number of variant ids");
            for (size_t i = 0; i < numLabels; i++) {
                result.emplace_back(std::move(readString(m_header.version, m_infile)));
            }
        }
        return std::move(result);
    }

private:
    size_t getVariantIndexOffset(VariantT variantIndex) {
        return m_header.filePosIndex + (variantIndex * sizeof(IndexEntry));
    }

    void readAndCheckHeader() {
        m_infile.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
        PICOVCF_GOOD_OR_MALFORMED_FILE(m_infile);
        PICOVCF_ASSERT_OR_MALFORMED(m_header.magic == IGD_MAGIC, "Invalid magic " << std::hex << m_header.magic);
        static_assert(CURRENT_IGD_VERSION == 4, "Remove/check backwards compatibility if incrementing version");
        PICOVCF_ASSERT_OR_MALFORMED(m_header.version == CURRENT_IGD_VERSION || m_header.version == 3,
                                    "Unsupported file format version " << m_header.version);
        m_source = readString(m_header.version, m_infile);
        m_description = readString(m_header.version, m_infile);
        m_beforeFirstVariant = m_infile.tellg();
        PICOVCF_ASSERT_OR_MALFORMED(m_header.filePosVariants > m_beforeFirstVariant,
                                    "Invalid variant info position " << m_header.filePosVariants);
        PICOVCF_ASSERT_OR_MALFORMED(m_header.ploidy <= MAX_PLOIDY,
                                    "Invalid ploidy " << m_header.ploidy << " is greater than maximum of "
                                                      << MAX_PLOIDY);
        PICOVCF_ASSERT_OR_MALFORMED(m_header.numIndividuals <= MAX_SAMPLES / m_header.ploidy,
                                    "Too many individuals to store: " << m_header.numIndividuals);
    }

    void readAllAlleleInfo() {
        m_infile.seekg(m_header.filePosVariants);
        PICOVCF_GOOD_OR_MALFORMED_FILE(m_infile);
        PICOVCF_RELEASE_ASSERT(m_referenceAlleles.empty());
        PICOVCF_RELEASE_ASSERT(m_alternateAlleles.empty());
        m_referenceAlleles.reserve(m_header.numVariants);
        m_alternateAlleles.reserve(m_header.numVariants);
        for (VariantT i = 0; i < m_header.numVariants; i++) {
            std::string reference = readString(m_header.version, m_infile);
            PICOVCF_GOOD_OR_MALFORMED_FILE(m_infile);
            IGDAllele refAllele(reference, m_longAlleles.size());
            if (refAllele.isLong()) {
                m_longAlleles.push_back(std::move(reference));
            }
            m_referenceAlleles.push_back(std::move(refAllele));

            std::string alternate = readString(m_header.version, m_infile);
            PICOVCF_GOOD_OR_MALFORMED_FILE(m_infile);
            IGDAllele altAllele(alternate, m_longAlleles.size());
            if (altAllele.isLong()) {
                m_longAlleles.push_back(std::move(alternate));
            }
            m_alternateAlleles.push_back(std::move(altAllele));
        }
        m_longAlleles.shrink_to_fit();
    }

    std::ifstream m_infile;
    FixedHeader m_header;
    std::string m_source;
    std::string m_description;
    std::streamoff m_beforeFirstVariant;

    // These are indexed by the "variant index" (0-based), and hold in memory all
    // of the reference/ alternate alleles. Even in an extremely large file, this
    // will not be much RAM.
    std::vector<IGDAllele> m_referenceAlleles;
    std::vector<IGDAllele> m_alternateAlleles;
    // Separate vector for storing the "long" alleles (larger than 4 nucleotides)
    std::vector<std::string> m_longAlleles;
};

template <typename T> static inline void writeScalar(T intValue, std::ostream& outStream) {
    outStream.write(reinterpret_cast<const char*>(&intValue), sizeof(intValue));
}

static inline void writeString(const std::string& value, std::ostream& outStream) {
    writeScalar<uint32_t>(value.size(), outStream);
    outStream.write(value.c_str(), value.size());
}

inline void writeAllelesAsOnes(std::ostream& outStream, const IGDSampleList& orderedSampleList, SampleT numSamples) {
    // XXX this can be majorly optimized.
    uint8_t buffer = 0x0;
    SampleT i = 0;
    SampleT j = 0;
    SampleT nextSample = SAMPLE_INDEX_NOT_SET;
    if (i < orderedSampleList.size()) {
        nextSample = orderedSampleList[i];
        i++;
    }
    while (j < numSamples) {
        if (j != 0 && j % 8 == 0) {
            writeScalar<uint8_t>(buffer, outStream);
        }
        buffer <<= 1;
        if (j == nextSample) {
            buffer |= 0x1;
            if (i < orderedSampleList.size()) {
                nextSample = orderedSampleList[i];
                i++;
            } else {
                nextSample = SAMPLE_INDEX_NOT_SET;
            }
        }
        j++;
    }
    SampleT leftover = numSamples % 8;
    if (leftover > 0) {
        buffer <<= (8 - leftover);
    }
    writeScalar<uint8_t>(buffer, outStream);
}

/**
 * Class for constructing IGD files on disk.
 */
class IGDWriter {
public:
    /**
     * Create an IGDWriter object, does not create the file.
     * @param[in] ploidy The fixed ploidy for all individuals.
     * @param[in] numIndividuals The number of individuals with genotype data.
     * @param[in] isPhased True if the data is phased, false otherwise.
     */
    IGDWriter(uint32_t ploidy, uint32_t numIndividuals, bool isPhased)
        : m_header({
              IGDData::IGD_MAGIC,
              IGDData::CURRENT_IGD_VERSION,
              ploidy,
              IGDData::DEFAULT_SPARSE_THRESHOLD,
              0,
              numIndividuals,
              0,
              isPhased ? IGDData::IGD_PHASED : 0x0,
              0,
              0,
              0,
              0,
          }) {
        PICOVCF_ASSERT_OR_MALFORMED(ploidy <= MAX_PLOIDY, "Ploidy exceeded maximum supported");
    }

    /**
     * Write the fixed-size IGD header to the output stream.
     * @param[in] outStream The output stream.
     * @param[in] source A description of where this data came from.
     * @param[in] description Generic description of the data.
     */
    void writeHeader(std::ostream& outStream, const std::string& source, const std::string& description) {
        outStream.write(reinterpret_cast<const char*>(&m_header), sizeof(m_header));
        writeString(source, outStream);
        writeString(description, outStream);
    }

    /**
     * Write a single variant (row) of sample/genotype data.
     *
     * This can be used for phased or unphased data: when phased the sample list always
     * consists of haploid sample indexes, when unphased the sample list consists of
     * individual indexes. For the former, indexes are 0 <= index < (ploidy*numIndividuals)
     * and for the latter 0 <= index < numIndividuals.
     *
     * @param[in] outStream The output stream.
     * @param[in] genomePosition The position of the variant on the genome.
     * @param[in] referenceAllele The string representing the reference allele.
     * @param[in] altAlleles The vector of strings representing the alternate
     *      alleles. Must be at least 1 alternate allele.
     * @param[in] sampleList A list of sample indexes. For phased diploid data, e.g.,
     *      index 0 is the 0th chromosome copy of the 0th individual, index 1
     *      is the 1st chromosome copy of the 0th individual, index 2 is the 0th
     *      chromosome copy of the 1st individual, etc. For unphased data, index 0
     *      is the 0th individual, index 1 is the 1st individual, etc.
     * @param[in] numCopes [Optional] Required only for unphased data. The number of
     *      copies of the alternate allele this variant represents, 0 < copies <= ploidy.
     * @param[in] isMissing [Optional] when set to true, this row represents all the
     *      samples that have missing data at the given polymorphic site.
     */
    void writeVariantSamples(std::ostream& outStream,
                             const uint64_t genomePosition,
                             const std::string& referenceAllele,
                             const std::string& altAllele,
                             const IGDSampleList& sampleList,
                             const bool isMissing = false,
                             const uint8_t numCopies = 0) {
        const bool phased = (bool)(m_header.flags & IGDData::IGD_PHASED);
        const SampleT numSamples = phased ? (m_header.ploidy * m_header.numIndividuals) : m_header.numIndividuals;
        const VariantT variantIndex = m_header.numVariants;
        assert(sampleList.size() <= numSamples);
        assert(phased || numCopies > 0);
        assert(phased || numCopies <= m_header.ploidy);

        m_referenceAlleles.emplace_back(referenceAllele);
        m_alternateAlleles.emplace_back(altAllele);

        const bool isSparse = (sampleList.size() <= (numSamples / IGDData::DEFAULT_SPARSE_THRESHOLD));
        m_index.push_back(std::move(makeEntry(genomePosition, isSparse, isMissing, numCopies, outStream.tellp())));
        if (isSparse) {
            writeSparse(outStream, sampleList);
            m_sparseCount++;
        } else {
            writeAllelesAsOnes(outStream, sampleList, numSamples);
        }
        m_header.numVariants++;
        m_totalCount++;
    }

    /**
     * Write the table of information about the variants. This information is
     * collected and saved by writeVariantSamples(), so this function must be
     * called _after_ that one.
     * @param[in] outStream The output stream.
     */
    void writeIndex(std::ostream& outStream) {
        m_header.filePosIndex = outStream.tellp();
        for (VariantT i = 0; i < m_index.size(); i++) {
            outStream.write((const char*)&m_index[i], sizeof(m_index[i]));
        }
    }

    /**
     * Write the table of information about the variants. This information is
     * collected and saved by writeVariantSamples(), so this function must be
     * called _after_ that one.
     * @param[in] outStream The output stream.
     */
    void writeVariantInfo(std::ostream& outStream) {
        m_header.filePosVariants = outStream.tellp();
        for (VariantT i = 0; i < m_header.numVariants; i++) {
            const std::string& ref = m_referenceAlleles.at(i);
            writeString(ref, outStream);

            const std::string& alt = m_alternateAlleles.at(i);
            writeString(alt, outStream);
        }
    }

    /**
     * Write the table of individual identifiers.
     *
     * This is an optional part of the IGD file. You can create an IGD and not
     * call this method, in which case the individual ID table will just be empty.
     *
     * @param[in] outStream The output stream.
     * @param[in] labels A list (vector) of string identifiers. One id per
     * individual.
     */
    void writeIndividualIds(std::ostream& outStream, const std::vector<std::string>& labels) {
        if (labels.empty()) {
            m_header.filePosIndividualIds = 0;
        } else {
            m_header.filePosIndividualIds = outStream.tellp();
            PICOVCF_RELEASE_ASSERT(0 != m_header.filePosIndividualIds);
            if (labels.size() != m_header.numIndividuals) {
                PICOVCF_THROW_ERROR(ApiMisuse, "Must provide one label per individual");
            }

            writeScalar<uint64_t>(labels.size(), outStream);
            for (const auto& label : labels) {
                writeString(label, outStream);
            }
        }
    }

    /**
     * Write the table of variant identifiers.
     *
     * This is an optional part of the IGD file. You can create an IGD and not
     * call this method, in which case the variant ID table will just be empty.
     *
     * @param[in] outStream The output stream.
     * @param[in] labels A list (vector) of string identifiers. One id per
     * variant.
     */
    void writeVariantIds(std::ostream& outStream, const std::vector<std::string>& labels) {
        if (labels.empty()) {
            m_header.filePosVariantIds = 0;
        } else {
            m_header.filePosVariantIds = outStream.tellp();
            PICOVCF_RELEASE_ASSERT(0 != m_header.filePosVariantIds);
            if (labels.size() != m_header.numVariants) {
                PICOVCF_THROW_ERROR(ApiMisuse, "Must provide one label per variant");
            }

            writeScalar<uint64_t>(labels.size(), outStream);
            for (const auto& label : labels) {
                writeString(label, outStream);
            }
        }
    }

    size_t m_sparseCount{};
    size_t m_totalCount{};

private:
    void writeSparse(std::ostream& outStream, const IGDSampleList& sampleList) {
        writeScalar<SampleT>(sampleList.size(), outStream);
        for (SampleT i = 0; i < sampleList.size(); i++) {
            writeScalar<SampleT>(sampleList[i], outStream);
        }
    }

    IGDData::IndexEntry makeEntry(VariantT genomePosition,
                                  const bool isSparse,
                                  const bool isMissing,
                                  uint8_t numCopies,
                                  std::streamoff filePosition) {
        uint64_t position = static_cast<uint64_t>(genomePosition);
        if (position > IGDData::MAX_BP_POSITION) {
            PICOVCF_THROW_ERROR(ApiMisuse, "Given genome position is too large: " << genomePosition);
        }
        if (numCopies > 0) {
            position |= ((uint64_t)numCopies) << IGDData::BP_POS_COPY_SHIFT;
        }
        if (isSparse) {
            position |= IGDData::BP_POS_FLAGS_SPARSE;
        }
        if (isMissing) {
            position |= IGDData::BP_POS_FLAGS_IS_MISSING;
        }
        return {position, static_cast<uint64_t>(filePosition)};
    }

    IGDData::FixedHeader m_header;
    std::vector<std::string> m_referenceAlleles;
    std::vector<std::string> m_alternateAlleles;
    std::vector<IGDData::IndexEntry> m_index;
};

/**
 * Using minimal memory, convert the given VCF file (can be gzipped) to an IGD
 * file with the given name.
 * @param[in] vcfFilename The name of the input VCF file to be converted.
 * @param[in] outFilename The name of the output IGD file to be created.
 * @param[in] description [Optional] A description of the dataset.
 * @param[in] verbose [Optional] Set to true to get statistics printed to
 * stdout.
 * @param[in] emitIndividualIds [Optional] Copy individual IDs to IGD file
 * (false by default).
 * @param[in] emitVariantIds [Optional] Copy variant IDs to IGD file (false by
 * default).
 * @param[in] forceUnphased [Optional] When true, force the result to be unphased
 * even if the input is phased (or mixed phased-ness).
 */
inline void vcfToIGD(const std::string& vcfFilename,
                     const std::string& outFilename,
                     std::string description = "",
                     bool verbose = false,
                     bool emitIndividualIds = false,
                     bool emitVariantIds = false,
                     bool forceUnphased = false) {
    VCFFile vcf(vcfFilename);
    vcf.seekBeforeVariants();
    PICOVCF_ASSERT_OR_MALFORMED(vcf.hasNextVariant(), "VCF file has no variants");
    vcf.nextVariant();
    VCFVariantView& variant1 = vcf.currentVariant();
    PICOVCF_ASSERT_OR_MALFORMED(variant1.hasGenotypeData(), "VCF file has no genotype data");

    IndividualIteratorGT firstIndividual = variant1.getIndividualIterator();
    PICOVCF_ASSERT_OR_MALFORMED(firstIndividual.hasNext(), "VCF file has no genotype data");
    VariantT allele1 = 0;
    VariantT allele2 = 0;
    const bool isPhased = forceUnphased || firstIndividual.getAlleles(allele1, allele2, /*moveNext=*/false);
    const uint64_t ploidy = (allele2 == NOT_DIPLOID) ? 1 : 2;

    std::vector<std::string> variantIds;
    std::ofstream outFile(outFilename, std::ios::binary);
    IGDWriter writer(ploidy, vcf.numIndividuals(), isPhased);
    writer.writeHeader(outFile, vcfFilename, description);
    vcf.seekBeforeVariants();
    while (vcf.hasNextVariant()) {
        vcf.nextVariant();
        VCFVariantView& variant = vcf.currentVariant();
        IndividualIteratorGT individualIt = variant.getIndividualIterator();
        auto altAlleles = variant.getAltAlleles();
        IGDSampleList missingData;
        const size_t numSampleLists = isPhased ? altAlleles.size() : (altAlleles.size() * ploidy);
        std::vector<IGDSampleList> variantGtData(numSampleLists);
        SampleT sampleIndex = 0;
        while (individualIt.hasNext()) {
            const bool isPhasedI = individualIt.getAlleles(allele1, allele2);
            PICOVCF_ASSERT_OR_MALFORMED(forceUnphased || isPhasedI == isPhased,
                                        "Cannot convert VCF with mixed phasedness, unless forceUnphased is set");
            const uint64_t ploidyI = (allele2 == NOT_DIPLOID) ? 1 : 2;
            PICOVCF_ASSERT_OR_MALFORMED(ploidyI == ploidy, "Cannot convert VCF with mixed ploidy");
            if (isPhased) {
                if (allele1 == MISSING_VALUE) {
                    missingData.push_back(sampleIndex);
                } else if (allele1 > 0) {
                    variantGtData.at(allele1 - 1).push_back(sampleIndex);
                }
                sampleIndex++;
                if (ploidy == 2) {
                    if (allele2 == MISSING_VALUE) {
                        missingData.push_back(sampleIndex);
                    } else if (allele2 > 0) {
                        variantGtData.at(allele2 - 1).push_back(sampleIndex);
                    }
                    sampleIndex++;
                }
                // sampleIndex refers to _haploid samples_ for unphased data.
            } else {
                if (allele1 == MISSING_VALUE || allele2 == MISSING_VALUE) {
                    missingData.push_back(sampleIndex);
                } else if (allele1 == allele2 && allele1 > 0) {
                    // This corresponds to a "2" (numCopies)
                    variantGtData.at((allele1 - 1) + altAlleles.size()).push_back(sampleIndex);
                } else {
                    // These correspond to a "1" (numCopies). We can have two "1"s in the case where the
                    // alt alleles are different. It is up to the downstream use-case how to handle this
                    // (e.g., filter the lower frequency one out, adjust the calculation to handle it, etc.)
                    if (allele1 > 0) {
                        variantGtData.at(allele1 - 1).push_back(sampleIndex);
                    }
                    if (allele2 > 0) {
                        variantGtData.at(allele2 - 1).push_back(sampleIndex);
                    }
                }
                sampleIndex++; // sampleIndex refers to _individuals_ for unphased data.
            }
        }
        std::string currentVariantId;
        if (emitVariantIds) {
            currentVariantId = variant.getID();
        }
        for (size_t altIndex = 0; altIndex < altAlleles.size(); altIndex++) {
            const auto position = variant.getPosition();
            const auto& ref = variant.getRefAllele();
            const auto& alt = altAlleles[altIndex];
            const uint8_t numCopies = isPhased ? 0 : 1;
            writer.writeVariantSamples(outFile, position, ref, alt, variantGtData[altIndex], false, numCopies);
            if (!isPhased && ploidy > 1) {
                const size_t copies2Index = altIndex + altAlleles.size();
                const auto& sampleList = variantGtData[copies2Index];
                if (!sampleList.empty()) {
                    writer.writeVariantSamples(outFile,
                                               position,
                                               ref,
                                               alt,
                                               sampleList,
                                               false,
                                               /*numCopies=*/2);
                    if (emitVariantIds) {
                        variantIds.emplace_back(currentVariantId);
                    }
                }
            }
            if (emitVariantIds) {
                variantIds.emplace_back(currentVariantId);
            }
        }
        if (!missingData.empty()) {
            writer.writeVariantSamples(outFile, variant.getPosition(), variant.getRefAllele(), "", missingData, true);
            if (emitVariantIds) {
                variantIds.emplace_back(currentVariantId);
            }
        }
    }
    writer.writeIndex(outFile);
    writer.writeVariantInfo(outFile);
    if (emitIndividualIds) {
        writer.writeIndividualIds(outFile, vcf.getIndividualLabels());
    }
    if (emitVariantIds) {
        writer.writeVariantIds(outFile, variantIds);
    }
    outFile.seekp(0);
    writer.writeHeader(outFile, vcfFilename, description);

    if (verbose) {
        std::cout << "Wrote " << writer.m_totalCount << " total variants" << std::endl;
        std::cout << "Of which " << writer.m_sparseCount << " were written sparsely" << std::endl;
    }
}

} // namespace picovcf

#endif /* PICOVCF_H */
