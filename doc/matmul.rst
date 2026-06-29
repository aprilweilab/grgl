
.. _matmul:

Matrix multiplication with GRG
==============================

GRG is a graph representation of the genotype matrix. GRG's internal nodes maintain a hierarchy 
between Mutations and samples, making it significantly more compact than even the sparse representation
of the genotype matrix. Non-GRG methods typically have to perform matrix multiplication block-wise,
because the dense/sparse genotype matrix cannot fit in RAM, whereas GRG does not usually require
block-wise computation.

Let :math:`N_H` be the number of haploid samples, :math:`N` be the number of individual samples, and
:math:`M` be the number of mutations. There are two different genotype matrices that we can consider:

- The haploid genotype matrix (we denote as :math:`X`) which has entries of either :math:`0` or :math:`1`, and dimensions :math:`N_H \times M`.
- The diploid genotype matrix (we denote as :math:`X_u`) which has entries of either :math:`0`, :math:`1`, or :math:`2`, and dimensions :math:`N \times M`.

GRGs can be constructed from phased data or unphased data. Haploid matrix operations can only be used with
phased data, but the diploid matrix operations can be used with either phased or unphased data.

For some applications, it might be useful to treat :math:`M` as the number of polymorphic sites as opposed to the
number of mutations. When the data is bi-allelic, these are equivalent. When the data is multi-allelic, users can
sum the results for the mutations that share sites.

See `the grapp preprint <https://doi.org/10.64898/2026.04.10.717786>`_ for a detailed description of how the
matrix multiplication works as a graph traversal. Practically speaking, given a matrix :math:`A` that is
of dimension :math:`K \times N` or a matrix :math:`B` that is :math:`K \times M`, you can use :py:meth:`pygrgl.matmul`
to get the matrix product :math:`A \times X_u` or :math:`B \times X_u^T` (or :math:`B \times X^T`).
The parameter `by_individual` controls whether the multiplication is against :math:`X_u` (`by_individual=True`)
or :math:`X` (`by_individual=False`).

See :py:meth:`pygrgl.matmul` for more details on usage.

Allele Frequency
----------------

Allele frequency calculates for each allele the value :math:`\frac{a}{N_H}`, where :math:`a` is the number of
samples containing the particular allele. This can be computed as a bottom-up graph traversal where every sample
node is given the initial value :math:`\frac{1}{N_H}` and the sums are propagated upwards until they reach mutation nodes.
The value saved at the mutation node is then the allele frequency.

To calculate allele frequency using :py:meth:`pygrgl.matmul` we construct an input vector :math:`\overrightarrow{v}` of length
:math:`N_H` with :math:`\frac{1}{N_H}` for every value. The resulting vector :math:`\overrightarrow{y}` (of length :math:`M`)
is the allele frequency for each mutation. Because ``matmul`` expects matrix input/output (not vector), we use row vectors.

::

	import numpy as np
	import pygrgl

	grg = pygrg.load_immutable_grg("test.grg")
	input = np.ones( (1, grg.num_samples) ) / grg.num_samples
	result = pygrgl.matmul(grg, input, pygrgl.TraversalDirection.UP)[0]

The above causes ``pygrgl`` to perform floating point addition at every node in the graph. Instead, we can just get the count
of every allele and then divide by ``grg.num_samples`` afterwards, which will be much faster (if we use the correct ``numpy.dtype``).

::

	import numpy as np
	import pygrgl

	grg = pygrg.load_immutable_grg("test.grg")
	input = np.ones( (1, grg.num_samples), dtype=np.uint32 )
	result = pygrgl.matmul(grg, input, pygrgl.TraversalDirection.UP)[0] / grg.num_samples

For very large datasets, using integers can make a non-trvial performance difference.


Hamming Distance
----------------

To compute the Hamming distance between two samples (i.e., how many variants differ between them), we can also use ``matmul``.
The operation that we can do efficiently is a 1-vs-all query, which can be used for finding the nearest neighbors to a sample,
for example. Given a GRG containing the sample set :math:`S`, consider a query sample (haplotype) :math:`q \in S`. Then to find
the Hamming distance between :math:`q` and all of :math:`S` we perform two ``matmul`` operations:

1. Create a vector :math:`v` which is all 0, except that :math:`v_q = 1` (the position associated with sample :math:`q`). Perform an upward ``matmul`` with :math:`v` to get result :math:`h`, which is :math:`q`'s haplotype.
2. Perform a ``matmul`` in the downward direction with input :math:`1 - (2 \times h)`. Given :math:`r`, the result of the matrix multiplication, the Hamming distance is then :math:`r + \sum h`.

::

	import numpy as np
	import pygrgl

	grg = pygrg.load_immutable_grg("test.grg")
	q = 0  # Use the 0th haplotype as the query.

	# Get q's haplotype
	v = np.zero( (1, grg.num_samples), dtype=np.int32 )
	v[q] = 1
	h = pygrgl.matmul(grg, v, pygrgl.TraversalDirection.UP)

	# Get all N_H hamming distances (note: hamming_q[q] = 0)
	r = pygrgl.matmul(grg, np.ones( (1, grgl.num_mutations) ) - (2 * h), pygrgl.TraversalDirection.DOWN)[0]
	hamming_q = r + numpy.sum(h)

Explanation
~~~~~~~~~~~

The Hamming distance between two binary sequences can be described in terms of sets of elements, where the presence of that element is a :math:`1` and absence is a :math:`0`.
Then :math:`H(A, B) = |A| + |B| - 2|A \cap B|`, i.e., the total number of :math:`1` s minus twice the shared :math:`1` s. Given :math:`q`, our query, and all other samples :math:`b`,
we get :math:`H(q, b) = |q| + |b| - 2|q \cap b|` where :math:`|b| - 2|q \cap b|` is computed by the second ``matmul`` above. The dot product between the query haplotype and the
genotype matrix, :math:`hX` produces the sum of the intersection between :math:`q` and each other haplotype. However, we perform :math:`(1 - 2h)X = 1X - 2hX` which gives us the
sum of each haplotype (i.e., :math:`|b|` for each :math:`b`) minus twice the intersection. Then we just need to add :math:`|q|` (``numpy.sum(h)``) to get the hamming distance.

Transformation of the genotype matrix
-------------------------------------

The `grapp <https://github.com/aprilweilab/grapp>`_ library provides instantiations of
`scipy.sparse.linalg.LinearOperator <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html>`_
that perform matrix multiplication against transformations of the genotype matrix. Some of the more useful
ones are:

1. `grapp.linalg.SciPyStdXOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXOperator>`_ - the standardized genotype matrix.
2. `grapp.linalg.SciPyStdXTXOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXTXOperator>`_ - the :math:`M \times M` covariance matrix based on the standardized genotype matrix.
3. `grapp.linalg.SciPyStdXXTOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXXTOperator>`_ - the :math:`N \times N` genetic related matrix (GRM), based on the standardized genotype matrix.

There are unstandardized versions of these operators as well.

Whole genome (all chromosomes)
------------------------------

The "multi" operators in `grapp <https://github.com/aprilweilab/grapp>`_ behave the same way as the other operators, except they operate on a list of
GRGs instead of a single one. Typically, a GRG is created for each chromosome, and then whole genome operations use a
`"multi" operator <https://grapp.readthedocs.io/en/latest/grapp.html#linear-operators-for-multiple-grgs>`_.
