#ifndef GRGL_BF_HASH_H
#define GRGL_BF_HASH_H

#include <cstdint>
#include <vector>

namespace grgl {

/**
 * Like a BFHash, but way paired down. Only a single hash function, and only exposes the bit
 * vector.
 */
class BFHash {
public:
    using ElemType = uint32_t;
    static constexpr size_t ELEM_SIZE = (sizeof(ElemType) * 8);

    explicit BFHash(size_t numBits)
            : m_bitVector((numBits + (ELEM_SIZE-1)) / ELEM_SIZE) {
    }

    void addHash(const size_t hashIdx) {
        const size_t K = m_bitVector.size() * ELEM_SIZE;
        const size_t bitIndex = hashIdx % K;
        const size_t element = bitIndex / ELEM_SIZE;
        const ElemType mask = 0x1U << (bitIndex % ELEM_SIZE);
        m_bitVector[element] |= mask;
    }

    const std::vector<ElemType>& vector() const {
        return m_bitVector;
    }

    std::vector<ElemType> stealVector() {
        return std::move(m_bitVector);
    }
private:
    std::vector<ElemType> m_bitVector;
};

}


#endif /* GRGL_BF_HASH_H */
