GRG Concepts
============

GRG
---

A GRG is a graph datastructure and a file format that captures the relationship between
a set of Mutations (variants) and samples. One way to think of this is as encoding the haploid
genotype matrix, which is a matrix where there is a row per haplotype sample and a column per
mutation/variant. Each value is either 0 (the sample does not have the variant allele) or 1
(the sample does have the variant allele). There can be multiple variants at a given genetic
position (site), and if a sample has 0 for all of them then that sample must have the reference
allele at that site. When this matrix is represented densely, it takes ``O(NM)`` memory, which
is huge for large genetic datasets. There are sparse matrix methods, which instead of representing
all ``N`` rows and ``M`` columns, represents pairs of indices ``i, j`` where the presence of
such a pair indicates a ``1`` and the absence of the pair ``i, j`` indicates that the i'th
sample has a ``0`` for the j'th variant.

A GRG represents this same matrix, but even more compactly than the sparse matrix representation.
It does this by capturing shared hierarchical subsets of samples that have the same mutations (variants).
A GRG is a multi-tree (the paths above any node form a tree, and the paths below any node form a tree).
There are two types of special nodes: mutation nodes, and sample nodes. A mutation node has one or more
mutations attached to it. A sample node represents a haplotype sample. A sample contains a mutation (has a ``1``
for that mutation in the genotype matrix) if and only if there is a path from the sample to a node
containing the mutation.

A GRG can be created from phased data (`IGD <https://github.com/aprilweilab/picovcf?tab=readme-ov-file#indexable-genotype-data-igd>`_,
`VCF <https://samtools.github.io/>`_, or `BGEN <https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html>`_)
or by converting an Ancestral Recombination Graph (in `tskit format <https://tskit.dev/tskit/docs/stable/introduction.html>`_) to a GRG.

Unphased data is also supported, but the resulting graph is much less optimal than for phased data. 

The "down edges" in a GRG go from mutations to samples, and are what is stored on disk. When requested,
the "up edges" are created during GRG load from disk, and allow bottom-up traversal of the graph in
arbitrary ways.

Once you have a GRG, computations are done by performing a `traversal <traversal.html>`_ of the graph
in either the upwards or downwards directions.

Samples
-------

Samples can have arbitrary ploidy in GRG. Given a ploidy of ``P``, and ``N`` individuals, there will
be ``P * N`` haplotypes stored in the GRG.

Sample Node
~~~~~~~~~~~

A node which has an ID between ``0...((N*P) - 1)``, and is a leaf in the GRG downward direction. The
haploid nodes for an individual are consecutive, so for example the first individual will have ``P``
nodes with IDs ``0, 1, ..., P-1``.

Individuals
~~~~~~~~~~~

Individuals can have identifiers associated with them, which helps keep track of properties on those
individuals when the dataset is tranformed or filtered. Individuals are not explicitly represented
in a GRG, but given a sample node with ID ``i`` you can find the associated individual "ID" by
dividing by the ploidy (``i / P``).

Populations
~~~~~~~~~~~

A GRG can have ``K`` populations associated with it. Each sample node then have a population ID
between ``0...(K-1)`` which can be used to lookup the population description.

Mutations
---------

Mutations are just a triple of ``(position, ref, alt)`` where ``position`` is a base-pair position (integer),
``ref`` is an arbitrary string for the reference allele, and ``alt`` is an arbitrary string for the alternate
allele. Mutations have IDs in the range ``0...(M-1)``. You can have multiple mutation IDs for "the same"
mutation (i.e., the ``(position, ref, alt)`` is the same).

Mutation Node
~~~~~~~~~~~~~

Any node in the graph can have a Mutation "attached" to it, at which point we call it a mutation node. There can
be more than one mutation attached to any given node. The mapping between mutation ID and node is 1-to-1, so
if you need to attach a mutation to multiple nodes you need to make a "copy" of the mutation (so it gets another
ID) and then associate the new ID with the additional node.

In practice, ``grg construct`` associates mutations with a single node. ``grg convert``, when converting an Ancestral
Recombination Graph (ARG) to a GRG *may* associate the same mutation with multiple nodes, depending on how the ARG
was created.

Mutation Time
~~~~~~~~~~~~~

Mutations can optionally have a time associated with them, indicating how old the mutation is. Units are unspecified,
and depend on the use-case.

Missing Data
------------

Missing data is represented as a node-to-samples relationship, just like Mutations are. However, there is no Mutation
object directly associated with the missing data. Instead, each Mutation is mapped to two (optional) ``NodeID``s: the
mutation's node, and the missingness node associated with the site (genetic position) of the Mutation.

Obtaining the Missingness Node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can obtain the missingness node for a Mutation the same way you obtain it's "regular" node: by iterating over all
Mutations. In C++, this is via :cpp:func:`GRG::getNodesAndMutations` (ordered by Mutation node) or 
:cpp:func:`GRG::getMutationsToNodeOrdered` (ordered by Mutation position/allele). In Python, this is via
:py:meth:`get_node_mutation_miss` and :py:meth:`get_mutation_node_miss`.

Given a node, you can also lookup the Mutations and their missingness nodes via :cpp:func:`GRG::getMutationsForNode` (C++)
or :py:meth:`get_muts_and_miss_for_node`.

Missing Data in Matrix Multiplication
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRG matrix multiplication handles missing data by an input/output vector (of length ``num_mutations``) for the value of
each missingness node associated with each Mutation. For multi-allelic sites, multiple Mutations will share a missingness
node, in which case the values provided for each such Mutation will be added together when applied to the missingness node.

Unphased Data
-------------

Unphased data is supported by GRG, though GRG construction is not yet optimized for unphased data. An unphased GRG is
identical to a phased GRG: all data is stored at the haploid level, and all samples represent haplotype samples. The
only difference is that :cpp:func:`GRG::isPhased` (:py:attr:`GRG.is_phased`) returns false. It is up to the client code
to ensure that all calculations are performed at the *individual* level (e.g., with the ``by_individual`` flag to
:py:meth:`matmul`).

Variations of GRGs
------------------

The typical GRG that most users/tools will interact with has the following properties:

1. It is immutable (:py:meth:`pygrgl.load_immutable_grg`), meaning the edges and nodes cannot be changed.
2. The nodes are numbered according to bottom-up topological order. That means a node with ID ``i`` is gauranteed to have all of it's children visited if you visit all nodes with IDs less than ``i``. This can be checked with :py:attr:`pygrgl.GRG.nodes_are_ordered`.
3. The :py:class:`pygrgl.Mutation`s are numbered in ascending order of ``(position, ref_allele, allele)``. This can be checked with :py:attr:`pygrgl.GRG.mutations_are_ordered`.
4. The Sample nodes are numbered ``0...(num_samples - 1)``, and are leaves of the graph. This can be checked with :py:attr:`pygrgl.GRG.samples_are_ordered`.

The above are true for any immutable GRG that has been loaded. However, you can (optionally) break the above properties
by performing modifications to the graph. For example:

* An immutable GRG can still change it's Mutations (just not the nodes/edges of the graph). So if you call :py:meth:`pygrgl.GRG.add_mutation` or :py:meth:`pygrgl.GRG.remove_mutation` then :py:attr:`pygrgl.GRG.mutations_are_ordered` will become ``False``.
* Changes to mutable GRGs can result in violating the above properties.
  * :py:meth:`pygrgl.MutableGRG.connect` can affect the node topological ordering (see :py:meth:`pygrgl.MutableGRG.connect` for details). If :py:attr:`pygrgl.GRG.nodes_are_ordered` is ``False``, then you need to use :py:attr:`pygrgl.GRG.get_topo_order` (with ``seed_list=None``) to do the traversal.
  * :py:meth:`pygrgl.MutableGRG.set_samples` will cause :py:attr:`pygrgl.GRG.samples_are_ordered` to become ``False``. You can no longer rely on the numbering of the samples, but this is not usually impactful, as :py:meth:`pygrgl.GRG.get_sample_nodes` and :py:meth:`pygrgl.GRG.is_sample` abstracts this away anyways.

If you have modified a GRG, you can always get all of these properties back by writing it to disk and reloading it
as an immutable GRG.

Negative nodes
~~~~~~~~~~~~~~

"Negative" nodes are just regular nodes in the GRG that have been created with ``grg.make_node(negative=True)``. This tells the GRG
that the node should be considered from the **bottom** of the graph for it's topological order. Regular nodes are always checked
against the **top** of the graph. For example:

* We create a regular node with ID ``A``. We attach ``A`` as a parent to an existing node ``B``. This maintains the topological order, because the NodeIDs increase in number, and that matched the topological order (``A`` is a parent of a node with a smaller ID)
* We create a regular node with ID ``A``. We attach ``A`` as a **child** to an existing node ``B``. This breaks topological order! ``B`` is a smaller NodeID than ``A``, but is its parent.
* We create a negative node with ID ``A``. We attach ``A`` as a child to an existing node ``B``. This maintains the topological order, as negative nodes are considered "less than" all regular nodes, for topological order.

The actual NodeID returned by :py:meth:`pygrgl.GRG.make_node` is still a positive value. The only thing that makes a node negative is how
the GRG tracks it internally. Also, when you connect a negative node to another node you must put a negative sign in front of it, like this:

::

    new_node = grg.make_node(negative=True)
    grg.connect(old_node, -new_node)         # Passed the negation! This keeps the topological order
    ...
    grg.get_up_edges(new_node)               # All other GRG APIs use the positive node ID


Advantages of immutable GRGs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* An immutable GRG uses 2-3x less RAM than a mutable GRG. 
* An immutable GRG has a faster matrix multiplication than a mutable GRG.
* :py:attr:`pygrgl.GRG.nodes_are_ordered` and :py:attr:`pygrgl.GRG.samples_are_ordered` will *always* be true in an immutable GRG.