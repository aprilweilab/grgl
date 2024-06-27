GRGL Documentation
==================

GRGL is a C++ library with Python bindings for accessing Genotype Representation Graphs (GRGs). 

.. contents::
   :depth: 2

.. toctree::
  :maxdepth: 2

  python_api
  cpp_api

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

