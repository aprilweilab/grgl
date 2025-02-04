
.. _dot_products:

Computing dot products with GRG
===============================

A GRG is a graph that represents the same data as the genotype matrix (as represented by a standard
tabular format like VCF). Since a GRG has internal nodes that maintain hierarchy and capture sharing
between Mutations and samples, it can be significantly more compact than even the sparse representation
of the genotype matrix. Furthermore, values can be stored these internal nodes, enabling dynamic-
programming like computations.

Let :math:`N_H` be the number of haploid samples, :math:`N` be the number of individual samples, and
:math:`M` be the number of mutations. There are two different genotype matrices that we can consider:

- The phased genotype matrix (we denote as :math:`X`) which has entries of either :math:`0` or :math:`1`, and dimensions :math:`N_H \times M`.
- The unphased diploid genotype matrix (we denote as :math:`X_u`) which has entries of either :math:`0`, :math:`1`, or :math:`2`, and dimensions :math:`N \times M`.

GRGs can only be constructed from phased data. However, we can still efficiently perform some operations on the
unphased genotype matrix :math:`X_u`.

For some applications, it might be useful to treat :math:`M` as the number of polymorphic sites as opposed to the
number of mutations. When the data is bi-allelic, these are equivalent. When the data is multi-allelic, users can
sum the results for the mutations that share sites.

Phased genotype matrix operations
---------------------------------

Dot product
~~~~~~~~~~~

The GRG directly captures the information in the phased genotype matrix :math:`X`. Traversals of the GRG
that *only sum* values calculated transitively for their children (upward traversal)
are equivalent to the matrix product :math:`\overrightarrow{v} \cdot X`. Similarly, traversals that *only sum*
values calculated transitively for parents (downward traversal) are equivalent to the matrix product
:math:`\overrightarrow{v} \cdot X^{T}`. :py:attr:`pygrgl.GRG.num_samples` is equal to :math:`N_H`, and
:py:attr:`pygrgl.GRG.num_mutations` is equal to :math:`M`.

**TODO: Expand this more.**

See :py:meth:`pygrgl.dot_product` for more details on usage.

Allele Frequency
^^^^^^^^^^^^^^^^

Allele frequency calculates for each allele the value :math:`\frac{a}{N_H}`, where :math:`a` is the number of
samples containing the particular allele. This can be computed as a bottom-up graph traversal where every sample
node is given the initial value :math:`\frac{1}{N_H}` and the sums are propagated upwards until they reach mutation nodes.
The value saved at the mutation node is then the allele frequency.

To calculate allele frequency using :py:meth:`pygrgl.dot_product` we construct an input vector :math:`\overrightarrow{v}` of length
:math:`N_H` with :math:`\frac{1}{N_H}` for every value. The resulting vector :math:`\overrightarrow{y}` (of length :math:`M`)
is the allele frequency for each mutation.

::

	import numpy as np
	import pygrgl

	grg = pygrg.load_immutable_grg("test.grg")
	input = np.ones(grg.num_samples) / grg.num_samples
	result = pygrgl.dot_product(grg, input, pygrgl.TraversalDirection.UP)


Hamming Distance
^^^^^^^^^^^^^^^^

**TODO**

Unphased genotype matrix operations
-----------------------------------

See :py:meth:`pygrgl.GRG.get_num_individual_coals`, which gets the number of individuals that
coalesced *at* the given node. See section "GWAS Computation" in the Methods of our paper
(`"Enabling efficient analysis of biobank-scale data with genotype representation graphs" <https://www.nature.com/articles/s43588-024-00739-9>`_)
for more details.

**TODO: Fill in this section. The above is a very vague pointer to how we handle this.**
