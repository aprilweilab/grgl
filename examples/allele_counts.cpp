// C++ example usage of GRGL
//
// First, build GRGL and then run "cmake --install . --prefix /some/path/you/choose/" to install
// the libraries and headers someplace you can easily use them.
//
// Build this example:
//   g++ allele_counts.cpp -I/some/path/you/choose/include -L/some/path/you/choose/lib -lgrgl -lvbyte -o allele_counts
//
// Run this example:
//   ./allele_counts <grg file>
//
//
#include <iostream>
#include "grgl/grg.h"
#include "grgl/serialize.h"

using namespace grgl;

int main(int argc, char *argv[]) {
    GRGPtr grg = loadImmutableGRG(argv[1]); // Load the GRG

    // Sum the allele counts over all nodes in bottom-up topo order
    std::vector<double> nodeSampleCounts(grg->numNodes());
    for (NodeID node = 0; node < grg->numNodes(); node++) {
        NodeIDSizeT count = 0;
        if (grg->isSample(node)) {
            count = 1;
        }
        for (NodeID child : grg->getDownEdges(node)) {
            count += nodeSampleCounts.at(child);
        }
        nodeSampleCounts[node] = count;
    }

    // Iterate over all the mutation nodes and print their counts
    for (const auto& nodeAndMut : grg->getNodesAndMutations()) {
        const NodeID nodeId = nodeAndMut.first;
        const MutationId mutId = nodeAndMut.second;
        Mutation mut = grg->getMutationById(mutId);
        NodeIDSizeT count = 0;
        if (nodeId != INVALID_NODE_ID) {
            count = nodeSampleCounts.at(nodeId);
        }
        std::cout << "Mutation(" << mutId << ": (" << mut.getPosition() << ", " <<
            mut.getAllele() << ") has sample count " << count << ")" << std::endl;
    }


    return 0;
}
