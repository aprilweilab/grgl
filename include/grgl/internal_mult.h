
template <typename T, bool useBitVector>
inline void vectorAdd(T* vect1, const T* vect2, const size_t base1, const size_t base2, const size_t count);

template <typename T, typename T2, bool useBitVector>
inline void
matmulPerformIOAddition(T* outputMatrix, const size_t outputIdx, const T2* inputMatrix, const size_t inputIdx);

template <typename T, typename T2, bool useBitVector>
inline void matmulPerformIOAddition(T* outputMatrix, const size_t outputIdx, const T2 inputValue, const size_t count);

template <typename T> class ValueSumVisitor : public GRGVisitor {
public:
    explicit ValueSumVisitor(std::vector<T>& nodeValues, const size_t numRows = 1)
        : m_nodeValues(nodeValues),
          m_numRows(numRows) {}

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        if (dfsPass != grgl::DfsPass::DFS_PASS_THERE) {
            const size_t base = nodeId * m_numRows;
            if (direction == DIRECTION_DOWN) {
                for (const auto& child : grg->getDownEdges(nodeId)) {
                    const size_t cbase = child * m_numRows;
                    vectorAdd<T, false>(&m_nodeValues[0], &m_nodeValues[0], base, cbase, m_numRows);
                }
            } else {
                for (const auto& parent : grg->getUpEdges(nodeId)) {
                    const size_t pbase = parent * m_numRows;
                    vectorAdd<T, false>(&m_nodeValues[0], &m_nodeValues[0], base, pbase, m_numRows);
                }
            }
        }
        return true;
    }

private:
    std::vector<T>& m_nodeValues;
    const size_t m_numRows;
};

// This showed no improvement in speed when using AVX, likely because of the time to load the
// data to the special registers, and then we only do a single operation at a time. Probably
// you have to do multiple operations to make it worthwhile.

template <typename T> inline void valueWiseVectorAdd(T* vect1, const T* vect2, const size_t count) {
#if USE_AVX
#error "Not implemented"
#endif
    for (size_t j = 0; j < count; j++) {
        vect1[j] += vect2[j];
    }
}

template <>
inline void vectorAdd<double, false>(
    double* vect1, const double* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void
vectorAdd<float, false>(float* vect1, const float* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<int64_t, false>(
    int64_t* vect1, const int64_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<int32_t, false>(
    int32_t* vect1, const int32_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<int16_t, false>(
    int16_t* vect1, const int16_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<int8_t, false>(
    int8_t* vect1, const int8_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<uint64_t, false>(
    uint64_t* vect1, const uint64_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<uint32_t, false>(
    uint32_t* vect1, const uint32_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<uint16_t, false>(
    uint16_t* vect1, const uint16_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<uint8_t, false>(
    uint8_t* vect1, const uint8_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    valueWiseVectorAdd(&vect1[base1], &vect2[base2], count);
}

template <>
inline void vectorAdd<uint8_t, true>(
    uint8_t* vect1, const uint8_t* vect2, const size_t base1, const size_t base2, const size_t count) {
    constexpr size_t elemBits = sizeof(uint8_t) * 8;
    // The two vectors have to be aligned to a byte. This means callers that make use of matMul() with bitvectors
    // have to ensure this. This is currently only exposed via the Python API, which guarantees this by padding out
    // the input array as needed.
    const size_t elemCount = count / elemBits;
    release_assert(count % elemBits == 0);
    const size_t base1E = base1 / elemBits;
    release_assert(base1 % elemBits == 0);
    const size_t base2E = base2 / elemBits;
    release_assert(base2 % elemBits == 0);
    size_t j = 0;
    for (; j < elemCount; j++) {
        vect1[base1E + j] ^= vect2[base2E + j];
    }
}

inline void validateMatMulInputs(const GRG* grg,
                                 const size_t inputCols,
                                 const size_t inputRows,
                                 const TraversalDirection direction,
                                 const size_t outputSize,
                                 const bool emitAllNodes,
                                 const bool byIndividual,
                                 const size_t outputCols) {
    release_assert(outputSize % inputRows == 0);
    const size_t expectSampleCols = byIndividual ? grg->numIndividuals() : grg->numSamples();
    if (direction == DIRECTION_DOWN) {
        if (inputCols != grg->numMutations()) {
            throw ApiMisuseFailure("Input vector is not size M (numMutations)");
        }
        if (!emitAllNodes && outputCols != expectSampleCols) {
            throw ApiMisuseFailure("Output vector is not size N (numSamples, or numIndividuals depending on flags)");
        }
    } else {
        if (inputCols != expectSampleCols) {
            throw ApiMisuseFailure("Input vector is not size N (numSamples, or numIndividuals depending on flags)");
        }
        if (!emitAllNodes && outputCols != grg->numMutations()) {
            throw ApiMisuseFailure("Output vector is not size M (numMutations)");
        }
    }
    if (emitAllNodes && outputCols != grg->numNodes()) {
        throw ApiMisuseFailure("Output vector is not size numNodes");
    }
}

template <typename T>
inline void
matmulPerformIOAdditionHelper(T* outputMatrix, const size_t outputIdx, const T* inputMatrix, const size_t inputIdx) {
    outputMatrix[outputIdx] += inputMatrix[inputIdx];
}

template <typename T>
inline void
matmulPerformIOAdditionHelper(T* outputMatrix, const size_t outputIdx, const T inputValue, const size_t count) {
    for (size_t i = 0; i < count; i++) {
        outputMatrix[outputIdx + i] += inputValue;
    }
}

template <>
inline void matmulPerformIOAddition<double, double, false>(double* outputMatrix,
                                                           const size_t outputIdx,
                                                           const double* inputMatrix,
                                                           const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<float, float, false>(float* outputMatrix,
                                                         const size_t outputIdx,
                                                         const float* inputMatrix,
                                                         const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<int64_t, int64_t, false>(int64_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int64_t* inputMatrix,
                                                             const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<int32_t, int32_t, false>(int32_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int32_t* inputMatrix,
                                                             const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<int16_t, int16_t, false>(int16_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int16_t* inputMatrix,
                                                             const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<int8_t, int8_t, false>(int8_t* outputMatrix,
                                                           const size_t outputIdx,
                                                           const int8_t* inputMatrix,
                                                           const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<uint64_t, uint64_t, false>(uint64_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint64_t* inputMatrix,
                                                               const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<uint32_t, uint32_t, false>(uint32_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint32_t* inputMatrix,
                                                               const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<uint16_t, uint16_t, false>(uint16_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint16_t* inputMatrix,
                                                               const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<uint8_t, uint8_t, false>(uint8_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const uint8_t* inputMatrix,
                                                             const size_t inputIdx) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputMatrix, inputIdx);
}

template <>
inline void matmulPerformIOAddition<double, double, false>(double* outputMatrix,
                                                           const size_t outputIdx,
                                                           const double inputValue,
                                                           const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<float, float, false>(float* outputMatrix,
                                                         const size_t outputIdx,
                                                         const float inputValue,
                                                         const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<int64_t, int64_t, false>(int64_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int64_t inputValue,
                                                             const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<int32_t, int32_t, false>(int32_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int32_t inputValue,
                                                             const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<int16_t, int16_t, false>(int16_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const int16_t inputValue,
                                                             const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<int8_t, int8_t, false>(int8_t* outputMatrix,
                                                           const size_t outputIdx,
                                                           const int8_t inputValue,
                                                           const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<uint64_t, uint64_t, false>(uint64_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint64_t inputValue,
                                                               const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<uint32_t, uint32_t, false>(uint32_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint32_t inputValue,
                                                               const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<uint16_t, uint16_t, false>(uint16_t* outputMatrix,
                                                               const size_t outputIdx,
                                                               const uint16_t inputValue,
                                                               const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

template <>
inline void matmulPerformIOAddition<uint8_t, uint8_t, false>(uint8_t* outputMatrix,
                                                             const size_t outputIdx,
                                                             const uint8_t inputValue,
                                                             const size_t count) {
    matmulPerformIOAdditionHelper(outputMatrix, outputIdx, inputValue, count);
}

// If the input bit is 0, leave output unchanged. If the input bit is 1, flip the output
// value to the other value.
template <typename T> inline void matmulFlipSingleBit(T* outputMatrix, const size_t outputIdx) {
    constexpr size_t elemBits = sizeof(T) * 8;
    const size_t elemOffset = outputIdx / elemBits;
    const size_t bitOffset = outputIdx % elemBits;
    const T mask = 0x1U << bitOffset;
    outputMatrix[elemOffset] ^= mask;
}

template <typename T> inline bool matmulGetSingleBit(const T* inputMatrix, const size_t inputIdx) {
    constexpr size_t elemBits = sizeof(T) * 8;
    const size_t elemOffset = inputIdx / elemBits;
    const size_t bitOffset = inputIdx % elemBits;
    const T mask = 0x1U << bitOffset;
    return (bool)(inputMatrix[elemOffset] & mask);
}

template <>
inline void matmulPerformIOAddition<uint8_t, bool, true>(uint8_t* outputMatrix,
                                                         const size_t outputIdx,
                                                         const bool* inputMatrix,
                                                         const size_t inputIdx) {
    if (inputMatrix[inputIdx]) {
        matmulFlipSingleBit(outputMatrix, outputIdx);
    }
}

template <>
inline void matmulPerformIOAddition<bool, uint8_t, true>(bool* outputMatrix,
                                                         const size_t outputIdx,
                                                         const uint8_t* inputMatrix,
                                                         const size_t inputIdx) {
    outputMatrix[outputIdx] = (outputMatrix[outputIdx] != matmulGetSingleBit(inputMatrix, inputIdx));
}

template <>
inline void matmulPerformIOAddition<uint8_t, bool, true>(uint8_t* outputMatrix,
                                                         const size_t outputIdx,
                                                         const bool inputValue,
                                                         const size_t count) {
    release_assert(false); // Unsupported
}