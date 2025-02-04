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
and this can result in significant compression (10-15x smaller than .vcf.gz).
Calculations on the whole dataset can be performed very quickly on GRG, using
GRGL. See our paper `"Enabling efficient analysis of biobank-scale data with
genotype representation graphs" <https://www.nature.com/articles/s43588-024-00739-9>`_
for more details.

.. contents::
   :depth: 2

.. toctree::
  :hidden:
  :maxdepth: 2

  Installation <installation>
  Constructing GRGs <construct>
  Converting Tree-Sequences <ts_convert>
  Graph Traversals <traversal>
  Calculating dot produces <dot_products>
  python_api
  cpp_api
  Examples and applications <examples_and_applications>

Simple Usage Examples
---------------------

Python: Load a GRG and compute allele frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The below example computes allele frequency by visiting every node in the graph. The graph
can be visited in topological order by just visiting the nodes in their natural NodeID order.
However, if you only want to visit a subset of the graph, you'll need to use one of the
visitor functions like :py:meth:`pygrgl.get_topo_order()`.

::

  import pygrgl

  # Load the graph from disk
  g = pygrgl.load_immutable_grg("my.grg")

  # Holds a mapping from node_id to the number of samples (transitively) beneath it.
  node_sample_counts = {}

  # The natural order of NodeIDs is guaranteed to follow a topological order in the graph:
  # all children will be visited prior to their parents. Reversing this list produces the
  # opposite where all parents are visited prior to their children.
  for node_id in range(0, g.num_nodes):
      count = 1 if g.is_sample(node_id) else 0
      for child_id in g.get_down_edges(node_id):
          count += node_sample_counts[child_id]
      node_sample_counts[node_id] = count

  # Iterate over all mutation nodes.
  for node_id, mutation_id in g.get_node_mutation_pairs():
      mutation = g.get_mutation_by_id(mutation_id)
      count = node_sample_counts[node_id]
      print(f"Mutation({mutation_id}): "
            f"({mutation.position}, {mutation.allele}) has sample count {count}")


C++ TBD1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  #include <iostream>
  #include "grgl/grg.hpp"

  int main(int argc, char *argv[]) {
    return 0;
  }


Simple Usage Examples
---------------------
