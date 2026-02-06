LinearOperators Tutorial
========================

Unlike some of the other tutorials, this one is focused on the Python
API only, and a generic mechanism to integrate with iterative linear
algebra methods from SciPy (or other libraries). The PCA and GWAS
functions actually make use of the
`LinearOperators <https://grapp.readthedocs.io/en/latest/grapp.html#non-standardized-linear-operators>`__
under the hood.

**What you’ll need:**

-  Python dependencies “grapp”, “pandas”, “seaborn”, “scipy”:
   ``pip install grapp pandas seaborn scipy``
-  Command line tool “wget”: ``sudo apt install wget`` (or your system’s
   equivalent)

Get Dataset
-----------

We’re going to use a really small existing (simulated) dataset so that
we can compare GRG to standard dense matrix algebra. As soon as datasets
get even a little bit big, the dense matrix can take up too much RAM for
demonstration purposes.

.. code:: bash

    %%bash
    
    if [[ ! -e linop.example.igd ]]; then
      # Download a small example dataset
      wget https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz -O linop.example.vcf.gz
    
      # Convert to IGD; this isn't necessary, but most of the time you will want to do this
      igdtools linop.example.vcf.gz -o linop.example.igd
    fi
    
    # Just show some stats about the dataset
    igdtools -s linop.example.igd


.. parsed-literal::

    --2026-02-06 12:39:52--  https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving github.com (github.com)... 140.82.112.4
    Connecting to github.com (github.com)|140.82.112.4|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz [following]
    --2026-02-06 12:39:52--  https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.110.133, 185.199.109.133, 185.199.108.133, ...
    Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.110.133|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 494022 (482K) [application/octet-stream]
    Saving to: ‘linop.example.vcf.gz’
    
         0K .......... .......... .......... .......... .......... 10% 2.07M 0s
        50K .......... .......... .......... .......... .......... 20% 2.52M 0s
       100K .......... .......... .......... .......... .......... 31% 4.32M 0s
       150K .......... .......... .......... .......... .......... 41% 5.12M 0s
       200K .......... .......... .......... .......... .......... 51% 2.94M 0s
       250K .......... .......... .......... .......... .......... 62% 3.31M 0s
       300K .......... .......... .......... .......... .......... 72%  119M 0s
       350K .......... .......... .......... .......... .......... 82% 4.37M 0s
       400K .......... .......... .......... .......... .......... 93% 3.36M 0s
       450K .......... .......... .......... ..                   100% 6.98M=0.1s
    
    2026-02-06 12:39:53 (3.75 MB/s) - ‘linop.example.vcf.gz’ saved [494022/494022]
    


.. parsed-literal::

    Wrote 10893 total variants
    Of which 6228 were written sparsely
    Stats for linop.example.igd
    ... in range 0 - 18446744073709551615
      Variants in range: 10893
      Average samples/var: 50.2075
      Stddev samples/var: 86.1985
      Average var/sample: 1367.28
      Stddev var/sample: 25.9211
      Variants with missing data: 0
      Total missing alleles: 0
      Total unique sites: 10885


Now construct the GRG.

.. code:: bash

    %%bash
    
    if [[ ! -e linop.example.grg ]]; then
      # -j controls how many threads to use.
      grg construct -j 1 linop.example.igd -o linop.example.grg
    fi


.. parsed-literal::

    Processing input file in 85 parts.
    Auto-calculating number of trees per part.
    Converting segments of input data to graphs
    100%|██████████| 85/85 [00:00<00:00, 134.66it/s]
    Merging...


.. parsed-literal::

    === GRG Statistics ===
    Nodes: 15481
    Edges: 93351
    Samples: 400
    Mutations: 10893
    Ploidy: 2
    Phased: true
    Populations: 0
    Range of mutations: 55829 - 9999127
    Specified range: 0 - 10894
    ======================
    Wrote simplified GRG with:
      Nodes: 15481
      Edges: 93351
    Wrote GRG to linop.example.grg


The first thing we’ll do is define a function that *very inneficiently*
converts a GRG to a dense matrix. This is only valuable for illustration
and testing purposes, as these matrices will get huge while a GRG will
remain small enough to fit in RAM still. For reasons of clarity, we’re
only going to generate the diploid genotype matrix (0, 1, 2 values) and
we’re going to pretend there is no such thing as missing data. For a
fuller implementation that handles these things, see
`testing_utils.py <https://github.com/aprilweilab/grapp/blob/main/test/testing_utils.py>`__.

.. code:: ipython3

    import pygrgl
    import numpy
    
    def grgToDiploidX(grg: pygrgl.GRG):
        """
        Create a diploid genotype matrix from a GRG.
        """
        result = numpy.zeros((grg.num_individuals, grg.num_mutations))
        samps_below = [list() for _ in range(grg.num_nodes)]
        for node_id in range(grg.num_nodes):
            sb = []
            if grg.is_sample(node_id):
                sb.append(node_id)
            for child_id in grg.get_down_edges(node_id):
                sb.extend(samps_below[child_id])
            samps_below[node_id] = sb
    
            muts = grg.get_mutations_for_node(node_id)
            if muts:
                for sample_id in sb:
                    indiv = sample_id // grg.ploidy
                    for mut_id in muts:
                        result[indiv][mut_id] += 1
        return result

.. code:: ipython3

    # Load the GRG dataset
    grg = pygrgl.load_immutable_grg("linop.example.grg", load_up_edges=False)
    
    # Convert GRG to dense matrix
    X = grgToDiploidX(grg)
    
    assert X.shape == (grg.num_individuals, grg.num_mutations)

The simplest thing to start with is just a random matrix multiplication
against the dense matrix and the
`grapp.linalg.ops_scipy.SciPyXOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyXOperator>`__.
This operator, and the others in ``grapp.linalg.ops_scipy`` conform to
the ``scipy`` `LinearOperator
interface <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html#scipy.sparse.linalg.LinearOperator>`__,
which means that ``_matmat()``, ``_rmatmat()``, ``_matvec()``, and
``_rmatvec()`` implement multiplications between the GRG and a matrix or
vector.

.. code:: ipython3

    from grapp.linalg.ops_scipy import SciPyXOperator
    
    # Perform a multiplication against X (NxM matrix, where N is number of individuals and M is number of mutations).
    random_mat = numpy.random.standard_normal((7, grg.num_individuals))
    
    # Dense matrix multiplication
    numpy_result = random_mat @ X
    
    # GRG matrix multiplication
    grg_X = SciPyXOperator(grg, pygrgl.TraversalDirection.UP)
    grg_result = random_mat @ grg_X
    
    numpy.allclose(numpy_result, grg_result)




.. parsed-literal::

    True



That very simple result is equivalent to using the
`pygrgl.matmul <https://grgl.readthedocs.io/en/stable/python_api.html#pygrgl.matmul>`__
function with ``by_individual=True``, it just has the nice feature of
supporting the ``numpy`` matrix multiplication operator ``@``, and also
can be used with ``scipy``\ ’s linear algebra methods.

You can also do matrix-vector products, or products with the transpose,
using ``numpy`` syntax.

.. code:: ipython3

    vect_result = random_mat[2] @ grg_X  # Vector x Matrix
    print(f"Vector result has shape {vect_result.shape}")
    
    B = numpy.random.standard_normal((3, grg.num_mutations))
    transpose_result = B @ grg_X.T  # B x X^T
    print(f"Transpose result has shape {transpose_result.shape}")


.. parsed-literal::

    Vector result has shape (10893,)
    Transpose result has shape (3, 200)


The real power of GRG LinearOperators comes from the transformations
they can make to the genotype matrix. Here we list some out: \*
`grapp.linalg.ops_scipy.SciPyStdXOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXOperator>`__:
Operations against the standardized genotype matrix. \*
`grapp.linalg.ops_scipy.SciPyStdXTXOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXTXOperator>`__:
Operations against the matrix :math:`X^T \times X` for standardized
:math:`X`, which is the correlation matrix of the genotype matrix. \*
`grapp.linalg.ops_scipy.SciPyStdXXTOperator <https://grapp.readthedocs.io/en/latest/grapp.html#grapp.linalg.ops_scipy.SciPyStdXXTOperator>`__:
Operations against the matrix :math:`X \times X^T` for standardized
:math:`X`, which is the genetic relatedness matrix (GRM).

Below we illustrate how these can interact with ``scipy``, by getting
the first :math:`k` eigenvalues from the GRM using
`scipy.sparse.linalg.eigs <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html#scipy.sparse.linalg.eigs>`__.

.. code:: ipython3

    from grapp.linalg.ops_scipy import SciPyStdXXTOperator
    from scipy.sparse.linalg import eigs
    from grapp.util import allele_frequencies
    
    # We need the allele frequencies for standardizing the genotype matrix
    freqs = allele_frequencies(grg)
    
    # Let's standardize the dense matrix X first, by column
    Xstd = numpy.zeros(X.shape)
    sigma = numpy.sqrt(2 * freqs * (1 - freqs))
    Xstd = numpy.divide(X - 2*freqs, sigma, out=Xstd, where=(sigma != 0))
    
    # First 5 eigen values of the dense GRM
    np_eigvals, _ = eigs(Xstd @ Xstd.T, k=5)
    
    # Now we can implicitly get the standardized GRM from the GRG by using a LinearOperator
    grg_GRM = SciPyStdXXTOperator(grg, freqs)
    grg_eigvals, _ = eigs(grg_GRM, k=5)
    
    assert numpy.allclose(np_eigvals, grg_eigvals)
    print(f"Dense result: {np_eigvals}")
    print(f"GRG result: {np_eigvals}")


.. parsed-literal::

    Dense result: [31851.93496578+0.j 28802.06816356+0.j 25479.32804366+0.j
     23571.51144699+0.j 21472.83667251+0.j]
    GRG result: [31851.93496578+0.j 28802.06816356+0.j 25479.32804366+0.j
     23571.51144699+0.j 21472.83667251+0.j]


Related Topics
--------------

-  `PCA <PCA.html>`__ is an example functionality that uses
   LinearOperators under the hood.
-  Documentation links:

   -  `grapp.linalg <https://grapp.readthedocs.io/en/latest/grapp.html#linear-algebra>`__:
      Python APIs for Linear Algebra
