#ifndef GRG_VARINT_H
#define GRG_VARINT_H

#include <limits>
#include <iostream>
#include <cstdint>

namespace grgl {

// https://sqlite.org/src4/doc/trunk/www/varint.wiki

inline uint8_t readByte(std::istream& inStream) {
    uint8_t byteValue = 0;
    inStream.read(reinterpret_cast<char*>(&byteValue), sizeof(byteValue));
    return byteValue;
}

inline uint64_t readVarInt(std::istream& inStream) {
    uint8_t controlByte = readByte(inStream);
    if (controlByte <= 240) {
        return static_cast<uint64_t>(controlByte);
    } else if (controlByte <= 248) {
        return 240 + (256 * (static_cast<uint64_t>(controlByte) - 241))
            + static_cast<uint64_t>(readByte(inStream));
    } else {
        switch (controlByte) {
            case 249:
                return 2288 + (256 * (static_cast<uint64_t>(readByte(inStream))))
                    + static_cast<uint64_t>(readByte(inStream));
            case 250:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U);
            case 251:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 24U);
            case 252:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 24U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 32U);
            case 253:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 24U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 32U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 40U);
            case 254:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 24U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 32U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 40U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 48U);
            default:
                return static_cast<uint64_t>(readByte(inStream))
                    + (static_cast<uint64_t>(readByte(inStream)) << 8U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 16U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 24U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 32U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 40U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 48U)
                    + (static_cast<uint64_t>(readByte(inStream)) << 56U);
        }
    }
}

inline void writeByte(uint8_t value, std::ostream& outStream) {
    outStream.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

inline void writeVarInt(uint64_t intValue, std::ostream& outStream) {
    const uint64_t case1 = 240;
    const uint64_t case2 = 2287;
    const uint64_t case3 = 67823;
    const uint64_t case6 = 1099511627775;
    const uint64_t case7 = 281474976710655;
    const uint64_t case8 = 72057594037927935;
    if (intValue <= case1) {
        writeByte(static_cast<uint8_t>(intValue), outStream);
    } else if (intValue <= case2) {
        writeByte(static_cast<uint8_t>((intValue-case1) / 256) + case1 + 1, outStream);
        writeByte(static_cast<uint8_t>((intValue-case1) % 256), outStream);
    } else if (intValue <= case3) {
        writeByte(static_cast<uint8_t>(249), outStream);
        writeByte(static_cast<uint8_t>((intValue - (case2+1)) / 256), outStream);
        writeByte(static_cast<uint8_t>((intValue - (case2+1)) % 256), outStream);
    } else if (intValue <= std::numeric_limits<uint16_t>::max()) {
        writeByte(static_cast<uint8_t>(250), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
    } else if (intValue <= std::numeric_limits<uint32_t>::max()) {
        writeByte(static_cast<uint8_t>(251), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 24), outStream);
    } else if (intValue <= case6) {
        writeByte(static_cast<uint8_t>(252), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 24), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 32), outStream);
    } else if (intValue <= case7) {
        writeByte(static_cast<uint8_t>(253), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 24), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 32), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 40), outStream);
    } else if (intValue <= case8) {
        writeByte(static_cast<uint8_t>(254), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 24), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 32), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 40), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 48), outStream);
    } else {
        writeByte(static_cast<uint8_t>(255), outStream);
        writeByte(static_cast<uint8_t>(intValue), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 8), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 16), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 24), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 32), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 40), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 48), outStream);
        writeByte(static_cast<uint8_t>(intValue >> 56), outStream);
    }
}

}

#endif /* GRG_VARINT_H */