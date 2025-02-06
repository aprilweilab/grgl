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
#include <args.hxx>
#include <fstream>
#include <iostream>
#include <string>

#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/serialize.h"
#include "grgl/visitor.h"
#include "util.h"

int main(int argc, char** argv) {
    args::ArgumentParser parser("Merge two GRGs.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Flag showStats(parser, "show-stats", "Show statistics about the GRG", {'s', "show-stats"});
    args::Flag noSimplify(parser, "no-simplify", "Simplify the resulting GRG", {'l', "no-simplify"});
    args::Flag noCombine(parser, "no-combine", "Do not combine nodes with same samplesets", {'c', "no-combine"});
    args::Positional<std::string> outfile(parser, "outfile", "The output file (resulting .grg)");
    args::PositionalList<std::string> inputs(parser, "inputs", "The input files (must be .grg)");
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
    if (!inputs || !outfile || inputs->empty()) {
        std::cout << parser;
        return 0;
    }

    std::string firstFile;
    std::list<std::string> otherFiles;
    for (const auto& filename : *inputs) {
        if (firstFile.empty()) {
            firstFile = filename;
        } else {
            otherFiles.push_back(filename);
        }
    }
    grgl::MutableGRGPtr grg1 = grgl::loadMutableGRG(firstFile);
    if (!grg1) {
        return 2;
    }

    grg1->merge(otherFiles, !noCombine);

    if (showStats) {
        dumpStats(grg1);
    }

    std::ofstream outStream(*outfile, std::ios::binary);
    auto counts = grgl::writeGrg(grg1, outStream, !noSimplify);
    std::cout << "Wrote simplified GRG with:" << std::endl;
    std::cout << "  Nodes: " << counts.first << std::endl;
    std::cout << "  Edges: " << counts.second << std::endl;
    std::cout << "Wrote GRG to " << *outfile << std::endl;
    return 0;
}
