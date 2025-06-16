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