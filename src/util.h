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
#ifndef GRGL_UTIL_H
#define GRGL_UTIL_H

#include "grgl/grgnode.h"
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

#define release_assert(condition)                                                                                      \
    do {                                                                                                               \
        if (!(condition)) {                                                                                            \
            std::cerr << "release_assert(" #condition ") failed at " << __FILE__ << ":" << __LINE__ << std::endl;      \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)

inline bool ends_with(std::string const& string1, std::string const& string2) {
    if (string1.length() < string2.length()) {
        return false;
    }
    return (string2 == string1.substr(string1.length() - string2.length()));
}

template <typename Out> inline void split(const std::string& s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return std::move(elems);
}

template <typename T> static double trailingMean(const std::vector<T>& hist, double tmThreshold) {
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

inline bool parseExactUint64(const std::string& exactValue, uint64_t& result) {
    char* endPtr = nullptr;
    auto _result = static_cast<uint64_t>(std::strtoull(exactValue.c_str(), &endPtr, 10));
    if (endPtr != (exactValue.c_str() + exactValue.size())) {
        return false;
    }
    result = _result;
    return true;
}

inline int32_t getEnvInt(const char* varName, int32_t defaultValue) {
    const char* value = getenv(varName);
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

template <class T> inline void helper_addToStringMap(T& result, const std::string& v1, const std::string& v2) {
    result.emplace(v1, v2);
};

template <>
inline void helper_addToStringMap(std::vector<std::pair<std::string, std::string>>& result,
                                  const std::string& v1,
                                  const std::string& v2) {
    result.emplace_back(v1, v2);
};

/**
 * Helper to convert data from a tab-separate file into a string->string map.
 */
template <typename MapType = std::map<std::string, std::string>>
inline MapType loadMapFromTSV(const std::string& filename,
                              const std::string& lhsField,
                              const std::string& rhsField,
                              const char separator = '\t') {
    release_assert(lhsField != rhsField);
    std::ifstream infile(filename);
    if (!infile) {
        std::stringstream ssErr;
        ssErr << "Could not read file " << filename;
        throw grgl::BadInputFileFailure(ssErr.str().c_str());
    }
    std::string line;
    size_t numCols = 0;
    size_t lineNum = 0;
    size_t lhsIndex = std::numeric_limits<size_t>::max();
    size_t rhsIndex = std::numeric_limits<size_t>::max();
    MapType result;
    while (std::getline(infile, line)) {
        auto tokens = split(line, separator);
        if (lineNum == 0) {
            numCols = tokens.size();
            for (size_t i = 0; i < tokens.size(); i++) {
                if (tokens[i] == lhsField) {
                    lhsIndex = i;
                }
                if (tokens[i] == rhsField) {
                    rhsIndex = i;
                }
            }
            std::stringstream ssErr;
            if (lhsIndex >= numCols) {
                ssErr << "Could not find TSV header named \"" << lhsField << "\". ";
            }
            if (rhsIndex >= numCols) {
                ssErr << "Could not find TSV header named \"" << rhsField << "\". ";
            }
            auto errorString = ssErr.str();
            if (!errorString.empty()) {
                throw std::runtime_error(errorString);
            }
        } else {
            if (numCols != tokens.size()) {
                std::stringstream ssErr;
                ssErr << "Malformed TSV file: wrong number of columns at line " << lineNum;
                throw std::runtime_error(ssErr.str());
            }
            helper_addToStringMap(result, tokens.at(lhsIndex), tokens.at(rhsIndex));
        }
        lineNum++;
    }
    return std::move(result);
}

inline bool pathExists(const std::string& pathname) {
    struct stat statBuf;
    return (stat(pathname.c_str(), &statBuf) == 0);
}

inline void makeDir(const std::string& pathname) {
    release_assert(0 == mkdir(pathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH));
}

inline std::string removeExt(const std::string& pathname) {
    size_t pos = pathname.find_last_of('.');
    if (pos != std::string::npos) {
        return std::move(pathname.substr(0, pos));
    }
    return pathname;
}

inline std::string basename(const std::string& pathname) {
    size_t pos = pathname.find_last_of('/');
    if (pos != std::string::npos) {
        return std::move(pathname.substr(pos + 1));
    }
    return pathname;
}

#endif /* GRGL_UTIL_H */
