
.. _traversal:

Graph Traversals
----------------

**TODO: Expand this into a more detailed description of graph traversals**

See :py:meth:`pygrgl.get_bfs_order`, :py:meth:`pygrgl.get_dfs_order`, and
:py:meth:`pygrgl.get_topo_order`. These methods take as input a list of *seed nodes* to
start traversing from, and emit a node ordering for all nodes reachable from those seeds.

When the :py:attr:`pygrgl.GRG.nodes_are_ordered` property is ``True``, the NodeIDs themselves provide
a topological order of the graph. That is, iterating from ``0...(grg.num_nodes-1)`` gives you
a bottom-up topological order of the entire graph, and iterating from ``(grg.num_nodes-1)...0``
gives you the top-down topological order. This property should be ``True`` for any GRG that has
been read from disk (until it is modified).

In the Python API, often using :py:meth:`pygrgl.dot_product` is significantly faster than performing
a traversal (even if you do not need to traverse the entire graph). See the `next section <dot_products.html>`_
for details on the dot product and when it is applicable.