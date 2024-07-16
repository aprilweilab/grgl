#ifndef GRGL_UTIL_H
#define GRGL_UTIL_H

#include "grgl/grgnode.h"
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define release_assert(condition) do { \
    if (!(condition)) { \
        std::cerr << "release_assert(" #condition ") failed at " << __FILE__ << ":" << __LINE__ << std::endl; \
        abort(); \
    } \
} while(0)

inline bool ends_with(std::string const &string1, std::string const &string2) {
    if (string1.length() < string2.length()) {
        return false;
    }
    return (string2 == string1.substr(string1.length() - string2.length()));
}

template <typename Out>
inline void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

template <typename T>
static double trailingMean(const std::vector<T>& hist, double tmThreshold) {
    release_assert(tmThreshold > 0.0 && tmThreshold <= 1.0);
    size_t totalValues = 0;
    for (size_t i = 0; i < hist.size(); i++) {
        totalValues += hist[i];
    }
    double mean = 0.0;
    const auto useValues = (size_t)((double)totalValues * tmThreshold);
    auto remainingValues = useValues;
    for (size_t i = 0; i < hist.size(); i++) {
        size_t val = 0;
        if (hist[i] < useValues) {
            remainingValues -= hist[i];
            val = hist[i];
        } else {
            val = remainingValues;
            remainingValues = 0;
        }
        mean += ((double)(i * val)) / ((double)useValues);
        if (remainingValues == 0) {
            break;
        }
    }
    return mean;
}

inline bool parseExactDouble(const std::string& exactValue, double& result) {
    char* endPtr = nullptr;
    double _result = std::strtod(exactValue.c_str(), &endPtr);
    if (endPtr != (exactValue.c_str() + exactValue.size())) {
        return false; 
    }
    result = _result;
    return true;
}

inline bool parseExactInt32(const std::string& exactValue, int32_t& result) {
    char* endPtr = nullptr;
    auto _result = static_cast<int32_t>(std::strtoul(exactValue.c_str(), &endPtr, 10));
    if (endPtr != (exactValue.c_str() + exactValue.size())) {
        return false; 
    }
    result = _result;
    return true;
}

inline bool parseExactUint32(const std::string& exactValue, uint32_t& result) {
    char* endPtr = nullptr;
    auto _result = static_cast<uint32_t>(std::strtoull(exactValue.c_str(), &endPtr, 10));
    if (endPtr != (exactValue.c_str() + exactValue.size())) {
        return false; 
    }
    result = _result;
    return true;
}

inline int32_t getEnvInt(const char* varName, int32_t defaultValue) {
    const char *value = getenv(varName);
    if (nullptr == value) {
        return defaultValue;
    }
    std::string valueStr(value);
    int32_t result = defaultValue;
    if (!parseExactInt32(valueStr, result)) {
        std::cerr << "Could not parse \"" << varName << "\" as integer" << std::endl;
    }
    return result;
}

inline grgl::NodeIDList loadNodeIDs(const std::string& filename) {
    grgl::NodeIDList result;
    std::ifstream infile(filename);
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) {
            continue;
        }
        uint32_t value = 0;
        if (!parseExactUint32(line, value)) {
            std::stringstream ssErr;
            ssErr << "Invalid unsigned integer in file" << std::endl;
            throw std::runtime_error(ssErr.str().c_str());
        }
        result.emplace_back(value);
    }
    return std::move(result);
}


#endif /* GRGL_UTIL_H */
