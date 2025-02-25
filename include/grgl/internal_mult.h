
template <typename T> inline void vectorAdd(T* vect1, const T* vect2, size_t count);

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
                    vectorAdd<T>(&m_nodeValues[base], &m_nodeValues[cbase], m_numRows);
                }
            } else {
                for (const auto& parent : grg->getUpEdges(nodeId)) {
                    const size_t pbase = parent * m_numRows;
                    vectorAdd<T>(&m_nodeValues[base], &m_nodeValues[pbase], m_numRows);
                }
            }
        }
        return true;
    }

private:
    std::vector<T>& m_nodeValues;
    const size_t m_numRows;
};

// TODO: unit test
template <> inline void vectorAdd<double>(double* vect1, const double* vect2, const size_t count) {
    size_t j = 0;
#if USE_AVX
    constexpr size_t batchSize = 4;
    __m256d dst, src;
    const size_t leftover = count & (batchSize - 1);
    for (size_t i = 0; i + batchSize <= count; i += batchSize) {
        dst = __builtin_ia32_loadupd256(&vect1[i]);
        src = __builtin_ia32_loadupd256(&vect2[i]);
        dst = __builtin_ia32_addpd256(dst, src);
        __builtin_ia32_storeupd256(&vect1[i], dst);
    }
    j = count - leftover;
#endif
    for (; j < count; j++) {
        vect1[j] += vect2[j];
    }
}

template <> inline void vectorAdd<float>(float* vect1, const float* vect2, const size_t count) {
    size_t j = 0;
#if USE_AVX
    constexpr size_t batchSize = 8;
    __m256 dst, src;
    const size_t leftover = count & (batchSize - 1);
    for (size_t i = 0; i + batchSize <= count; i += batchSize) {
        dst = __builtin_ia32_loadups256(&vect1[i]);
        src = __builtin_ia32_loadups256(&vect2[i]);
        dst = __builtin_ia32_addps256(dst, src);
        __builtin_ia32_storeups256(&vect1[i], dst);
    }
    j = count - leftover;
#endif
    for (; j < count; j++) {
        vect1[j] += vect2[j];
    }
}

template <> inline void vectorAdd<int64_t>(int64_t* vect1, const int64_t* vect2, const size_t count) {
    size_t j = 0;
#if USE_AVX
    constexpr size_t batchSize = 2;
    __m128i dst, src;
    const size_t leftover = count & (batchSize - 1);
    for (size_t i = 0; i + batchSize <= count; i += batchSize) {
        dst = _mm_loadu_si128((__m128i*)&vect1[i]);
        src = _mm_loadu_si128((__m128i*)&vect2[i]);
        dst = _mm_add_epi64(dst, src); // Cycle 2: [A AB BC CD]
        _mm_storeu_si128((__m128i*)&vect1[i], dst);
    }
    j = count - leftover;
#endif
    for (; j < count; j++) {
        vect1[j] += vect2[j];
    }
}

template <> inline void vectorAdd<int32_t>(int32_t* vect1, const int32_t* vect2, const size_t count) {
    size_t j = 0;
#if USE_AVX
    constexpr size_t batchSize = 4;
    __m128i dst, src;
    const size_t leftover = count & (batchSize - 1);
    for (size_t i = 0; i + batchSize <= count; i += batchSize) {
        dst = _mm_loadu_si128((__m128i*)&vect1[i]);
        src = _mm_loadu_si128((__m128i*)&vect2[i]);
        dst = _mm_add_epi32(dst, src); // Cycle 2: [A AB BC CD]
        _mm_storeu_si128((__m128i*)&vect1[i], dst);
    }
    j = count - leftover;
#endif
    for (; j < count; j++) {
        vect1[j] += vect2[j];
    }
}