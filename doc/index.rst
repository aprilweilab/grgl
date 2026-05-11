GRGL Overview
=============

GRGL is a C++ library with Python bindings for accessing Genotype Representation Graphs (GRGs). 
A GRG is a compact way to store reference-aligned genotype data for large
genetic datasets. These datasets are typically stored in tabular formats (VCF,
BCF, BGEN, etc.) and then compressed using off-the-shelf compression. In
contrast, a GRG contains Mutation nodes (representing variants) and Sample nodes
(representing haploid samples), where there is a path from a Mutation node to a
Sample node if-and-only-if that sample contains that mutation. These paths go
through internal nodes that represent common ancestry between multiple samples,
and this can result in significant compression (20-40x smaller than ``.vcf.gz``).
Calculations on the whole dataset can be performed very quickly on GRG, using
GRGL. See our paper `"Enabling efficient analysis of biobank-scale data with
genotype representation graphs" <https://www.nature.com/articles/s43588-024-00739-9>`_
for more details.

This documentation tries to pull together all the important GRG-related functionality
in one spot. If you want detailed API documentation, they can be found at:

* Core library: `pygrgl <https://grgl.readthedocs.io/en/stable/python_api.html>`_
* Statgen and popgen applications: `grapp <https://grapp.readthedocs.io/en/latest/>`_
* Phenotype simulation: `grg_pheno_sim <https://grg-pheno-sim.readthedocs.io/en/latest/>`_


.. toctree::
  :caption: Concepts and APIs
  :maxdepth: 1

  Installation <installation>
  Concepts <concepts>
  Constructing GRGs <construct>
  Converting Tree-Sequences <ts_convert>
  Graph Traversals <traversal>
  Calculating dot products <dot_products>
  python_api
  cpp_api
  Examples and applications <examples_and_applications>
  Command-line Recipes <cli_recipes>

.. toctree::
  :caption: Tutorials
  :maxdepth: 1

  tutorials/GWAS
  tutorials/GWASCovariates
  tutorials/PCA
  tutorials/LinearOperators
  tutorials/WorkingWithRealData
  tutorials/WorkingWithSimData
  tutorials/SimulatingPhenotypes
  tutorials/IGDToGRG
  tutorials/VCFToGRG


Simple Usage Examples
---------------------

Python: Load a GRG and compute allele frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The below example computes allele frequency by visiting every node in the graph. The graph
can be visited in topological order by just visiting the nodes in their natural NodeID order.
However, if you only want to visit a subset of the graph, you'll need to use one of the
visitor functions like :py:meth:`pygrgl.get_topo_order()`. We could just use
`grapp.util.simple.allele_counts <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.util.simple.allele_counts>`_
to calculate this, or use :py:meth:`pygrgl.matmul()` to compute it much more quickly, but
the below illustrates more concepts about how the graph can be used.

::

  import pygrgl
  import sys

  # Load the graph from disk
  g = pygrgl.load_immutable_grg(sys.argv[1])

  # Holds a mapping from node_id to the number of samples (transitively) beneath it.
  node_sample_counts = [0 for _ in range(g.num_nodes)]

  # The natural order of NodeIDs follows a topological order in the graph: all children
  # will be visited prior to their parents. Reversing this list produces the opposite
  # where all parents are visited prior to their children.
  for node_id in range(g.num_nodes):
      count = 1 if g.is_sample(node_id) else 0
      for child_id in g.get_down_edges(node_id):
          count += node_sample_counts[child_id]
      node_sample_counts[node_id] = count

  # Iterate over all mutation nodes.
  for node_id, mutation_id in g.get_node_mutation_pairs():
      mutation = g.get_mutation_by_id(mutation_id)
      if node_id != pygrgl.INVALID_NODE:
          count = node_sample_counts[node_id]
      else:
          count = 0
      print(f"Mutation({mutation_id}): "
            f"({mutation.position}, {mutation.allele}) has sample count {count}")


C++: Load a GRG and compute allele frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build C++ code that uses GRGL, you need to first build GRGL.

::

  mkdir build && cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release
  make -j4

Then, install the library into a folder where you can access the headers and static libraries.
(Alternatively, you could include GRGL as a CMake project in another CMake project.)

::

  cmake --install . --prefix /some/path/you/choose/

Finally, you can build the example below using a modern C++ compiler. You need to link against
libvbyte (which is built and included with the GRGL install above).

::

  g++ allele_counts.cpp -I/some/path/you/choose/include -L/some/path/you/choose/lib -lgrgl -lvbyte -o allele_counts

Here is the example code, which computes the same thing as the Python example above:

::

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

