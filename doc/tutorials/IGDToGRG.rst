Create IGDs and Converting IGD to GRG
=====================================

This tutorial covers two topics: 1. Creating an IGD file from a
``.vcf.gz``. 2. Converting an IGD file to a GRG

IGD files store referenced-aligned genotype data in a very simple sparse
matrix format, without compression. An IGD file is often smaller than
the corresponding ``.vcf.gz`` file, and is almost always faster (usually
``15x`` or more) to access. See `the
paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC11838554/>`__ if you’re
interested in comparison experiments.

Many datasets are stored in ``.vcf.gz`` format by “default.” If these
datasets are large, they are usually stored using
`BGZIP <https://www.htslib.org/doc/bgzip.html>`__ so that they can be
indexed for semi-random access. The two different kinds of index files
for ``BGZIP`` are `tabix <https://www.htslib.org/doc/tabix.html>`__ or
`bcftools <https://samtools.github.io/bcftools/bcftools.html>`__. We
only support ``tabix``-style indexes.

**What you’ll need:**

-  Python dependencies “pygrgl” and “igdtools”:
   ``pip install pygrgl igdtools``
-  Command line tool “wget” and “tabix”: ``sudo apt install wget tabix``
   (or your distribution’s equivalent)

Get Dataset
~~~~~~~~~~~

For our example, we’ll just download a very small simulated dataset that
is stored as ``.vcf.gz``.

.. code:: bash

    %%bash
    
    # Download a small example dataset
    wget https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz -O igd_convert.example.vcf.gz


.. parsed-literal::

    --2026-05-11 09:49:04--  https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving github.com (github.com)... 140.82.112.3
    Connecting to github.com (github.com)|140.82.112.3|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz [following]
    --2026-05-11 09:49:04--  https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.108.133, 185.199.111.133, 185.199.110.133, ...
    Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.108.133|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 494022 (482K) [application/octet-stream]
    Saving to: ‘igd_convert.example.vcf.gz’
    
         0K .......... .......... .......... .......... .......... 10% 2.25M 0s
        50K .......... .......... .......... .......... .......... 20% 6.40M 0s
       100K .......... .......... .......... .......... .......... 31% 5.05M 0s
       150K .......... .......... .......... .......... .......... 41% 7.50M 0s
       200K .......... .......... .......... .......... .......... 51% 6.98M 0s
       250K .......... .......... .......... .......... .......... 62% 6.04M 0s
       300K .......... .......... .......... .......... .......... 72% 12.1M 0s
       350K .......... .......... .......... .......... .......... 82% 11.2M 0s
       400K .......... .......... .......... .......... .......... 93% 21.0M 0s
       450K .......... .......... .......... ..                   100% 21.6M=0.07s
    
    2026-05-11 09:49:04 (6.47 MB/s) - ‘igd_convert.example.vcf.gz’ saved [494022/494022]
    


Convert to IGD
--------------

If we want to convert from ``.vcf.gz`` to IGD with a single thread, then
we do not need a ``tabix`` index. However, if we want to use multiple
threads (i.e., ``-j 2`` below) then a tabix index will provide a
performance improvement. So here we first ``tabix`` index the dataset,
and then convert to IGD with ``igdtools``.

.. code:: bash

    %%bash
    
    # Index the file. Note: usually your dataset will come with an index already.
    tabix igd_convert.example.vcf.gz
    
    # -j controls how many threads to use.
    igdtools -j 2 igd_convert.example.vcf.gz -o igd_convert.example.igd


.. parsed-literal::

    Wrote 5446 total variants
    Of which 3170 were written sparsely
    Wrote 5447 total variants
    Of which 3058 were written sparsely


Convert to GRG
--------------

If we have an IGD file (either because we converted it, or that is how
our dataset came) then we can construct a GRG easily:

.. code:: bash

    %%bash
    
    # -j controls how many threads to use.
    grg construct -j 1 igd_convert.example.igd -o igd_convert.example.grg


.. parsed-literal::

    Processing input file in 85 parts.
    Auto-calculating number of trees per part.
    Converting segments of input data to graphs
    100%|██████████| 85/85 [00:00<00:00, 203.33it/s]
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
    Wrote GRG to igd_convert.example.grg


What if I have metadata?
------------------------

IGD and GRG have the general philosophy that a lot of metadata is kept
separate from the genotype data. The metadata that the two formats
contain natively are: \* IGD: Individual identifiers (one per sampled
individual, not haplotype) - see
`pyigd.IGDReader.get_individual_ids <https://pyigd.readthedocs.io/en/latest/pyigd.html#pyigd.IGDReader.get_individual_ids>`__.
\* IGD: Variant identifiers (one per variant) - see
`pyigd.IGDReader.get_variant_ids <https://pyigd.readthedocs.io/en/latest/pyigd.html#pyigd.IGDReader.get_variant_ids>`__.
\* GRG: Individual identifiers (same as IGD) - see
`pygrgl.GRG.get_individual_id <https://grgl.readthedocs.io/en/stable/python_api.html#pygrgl.GRG.get_individual_id>`__.

Beyond that, it is suggested that you keep metadata in a simple format.
For example, ``igdtools`` supports `exporting metadata to .txt
files <https://picovcf.readthedocs.io/en/latest/igdtools.html#convert-vcf-gz-to-igd-and-export-metadata>`__
in a format that is loadable by
`numpy.loadtxt <https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html>`__.
This metadata can be accessed by the index of the variant (mutation),
:math:`i`, and you can also keep a mapping from variant identifier
:math:`v_i` to :math:`i` so that you can easily lookup other metadata by
the variant id (``igdtools`` exports such a mapping).

Related Topics
--------------

-  The `igdtools
   documentation <https://picovcf.readthedocs.io/en/latest/igdtools.html>`__
-  An overview of the `IGD file
   format <https://picovcf.readthedocs.io/en/latest/igd_overview.html>`__
