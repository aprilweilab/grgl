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

/**
 * Index input file formats that are not already indexable.
 * Creates a new file with the same basename as the input file.
 *
 * Usage:
 *  gindexer <infile name>
 */
#include "util.h"
#include <args.hxx>
#include <iostream>
#include <memory>
#include <sstream>

#if BGEN_ENABLED
extern "C" {
#include <bgen/bgen.h>
}
#endif

#include <sys/stat.h>

inline bool fileExists(const std::string& filename) {
    struct stat statBuffer;
    return stat(filename.c_str(), &statBuffer) == 0;
}

#if BGEN_ENABLED
static bool bgenIndex(const std::string& filename) {
    std::stringstream metafileSS;
    metafileSS << filename << ".meta";
    std::string metafilename = metafileSS.str();

    if (fileExists(metafilename)) {
        std::cerr << "BGEN index already exists at " << metafilename << std::endl;
        return false;
    }

    bgen_file* bgenFile = bgen_file_open(filename.c_str());
    if (nullptr == bgenFile) {
        std::cerr << "Failed to open " << filename << std::endl;
        return false;
    }

    bgen_metafile* metafile = bgen_metafile_create(bgenFile, metafilename.c_str(), 1, 0);
    release_assert(metafile != nullptr);
    bgen_metafile_close(metafile);
    bgen_file_close(bgenFile);
    return true;
}
#endif

int main(int argc, char* argv[]) {
    args::ArgumentParser parser("Convert between input file formats.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> infile(parser, "infile", "The input file (.vcf or .bgen)");
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (!infile) {
        std::cout << parser;
        return 0;
    }

#if BGEN_ENABLED
    if (ends_with(*infile, ".bgen")) {
        return bgenIndex(*infile) ? 0 : 2;
    }
#endif

    std::cerr << "Unknown filetype; the file extension is used to determine type." << std::endl;
    std::cerr << "Supported filetypes: *.bgen" << std::endl;
    return 1;
}
