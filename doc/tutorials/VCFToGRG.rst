Converting .vcf.gz to GRG
=========================

Many datasets are stored in ``.vcf.gz`` format by “default.” If these
datasets are large, they are usually stored using
`BGZIP <https://www.htslib.org/doc/bgzip.html>`__ so that they can be
indexed for semi-random access. The two different kinds of index files
for ``BGZIP`` are `tabix <https://www.htslib.org/doc/tabix.html>`__ or
`bcftools <https://samtools.github.io/bcftools/bcftools.html>`__.

In this tutorial we’ll show the (very simple) process of converting
``.vcf.gz`` data to ``GRG`` format, which is much smaller (usually at
least ``25x`` smaller) and faster (many orders of magnitude) for
performing computations.

**What you’ll need:** \* Python dependencies “pygrgl”:
``pip install pygrgl`` \* Command line tools “wget” and “tabix”:
``sudo apt install wget tabix`` (or your distribution’s equivalent)

Get Dataset
~~~~~~~~~~~

For our example, we’ll just download a very small simulated dataset that
is stored as ``.vcf.gz``.

.. code:: bash

    %%bash
    
    # Download a small example dataset
    wget https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz -O vcf_convert.example.vcf.gz


.. parsed-literal::

    --2026-02-06 12:31:22--  https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving github.com (github.com)... 140.82.114.4
    Connecting to github.com (github.com)|140.82.114.4|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz [following]
    --2026-02-06 12:31:23--  https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.111.133, 185.199.110.133, 185.199.108.133, ...
    Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.111.133|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 494022 (482K) [application/octet-stream]
    Saving to: ‘vcf_convert.example.vcf.gz’
    
         0K .......... .......... .......... .......... .......... 10% 1.87M 0s
        50K .......... .......... .......... .......... .......... 20% 4.30M 0s
       100K .......... .......... .......... .......... .......... 31% 3.73M 0s
       150K .......... .......... .......... .......... .......... 41% 3.93M 0s
       200K .......... .......... .......... .......... .......... 51% 4.32M 0s
       250K .......... .......... .......... .......... .......... 62% 4.08M 0s
       300K .......... .......... .......... .......... .......... 72% 4.49M 0s
       350K .......... .......... .......... .......... .......... 82% 4.11M 0s
       400K .......... .......... .......... .......... .......... 93% 4.59M 0s
       450K .......... .......... .......... ..                   100% 4.37M=0.1s
    
    2026-02-06 12:31:23 (3.71 MB/s) - ‘vcf_convert.example.vcf.gz’ saved [494022/494022]
    


Convert to GRG
--------------

Lets first attempt to convert to GRG without having an index for the
file.

.. code:: bash

    %%bash
    
    # -j controls how many threads to use.
    grg construct -j 1 vcf_convert.example.vcf.gz -o vcf_convert.example.grg || true


.. parsed-literal::

    Will not count variants in VCF files (too slow)
    Could not count number of variants in vcf_convert.example.vcf.gz. Using the default of 100 (use --parts to override).
    Processing input file in 100 parts.
    Auto-calculating number of trees per part.
    Converting segments of input data to graphs
      0%|          | 0/100 [00:00<?, ?it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    terminate called after throwing an instance of 'grgl::ApiMisuseFailure'
      what():  WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      0%|          | 0/100 [00:00<?, ?it/s]
    multiprocessing.pool.RemoteTraceback: 
    """
    Traceback (most recent call last):
      File "/usr/lib/python3.10/multiprocessing/pool.py", line 125, in worker
        result = (True, func(*args, **kwds))
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 300, in star_build_grg
        return build_grg(*args)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 266, in build_grg
        shape_grg = build_shape(range_triple, args, auto_args, input_file, output_file)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 253, in build_shape
        tb_time = time_call(command, stdout=sys.stdout)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/common.py", line 48, in time_call
        subprocess.check_call(command, **kwargs)
      File "/usr/lib/python3.10/subprocess.py", line 369, in check_call
        raise CalledProcessError(retcode, cmd)
    subprocess.CalledProcessError: Command '['/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/../../grgl', 'vcf_convert.example.vcf.gz', '--trees', 'optimal', '--lf-no-tree', '10', '--reduce', '5', '-r', '0.0:0.01', '-o', 'vcf_convert.example.grg.part0.grg']' died with <Signals.SIGABRT: 6>.
    """
    
    The above exception was the direct cause of the following exception:
    
    Traceback (most recent call last):
      File "/home/ddehaas/Py3Env/bin/grg", line 8, in <module>
        sys.exit(main())
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/cli.py", line 70, in main
        construct.from_tabular(args)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 421, in from_tabular
        list(
      File "/home/ddehaas/Py3Env/lib/python3.10/site-packages/tqdm/std.py", line 1181, in __iter__
        for obj in iterable:
      File "/usr/lib/python3.10/multiprocessing/pool.py", line 873, in next
        raise value
    subprocess.CalledProcessError: Command '['/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/../../grgl', 'vcf_convert.example.vcf.gz', '--trees', 'optimal', '--lf-no-tree', '10', '--reduce', '5', '-r', '0.0:0.01', '-o', 'vcf_convert.example.grg.part0.grg']' died with <Signals.SIGABRT: 6>.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    terminate called after throwing an instance of 'grgl::ApiMisuseFailure'
      what():  WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.


The ``grg`` tool did not like that. Why?
``WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.``

For large datasets, trying to convert a ``.vcf.gz`` to GRG without an
index will be very slow. This warning (and failure) is there to prevent
you from accidentally doing this. However, we know that our dataset is
really small so we don’t care – we can use the ``--force`` flag to force
GRG construction.

.. code:: bash

    %%bash
    
    # -j controls how many threads to use.
    grg construct --force -j 1 vcf_convert.example.vcf.gz -o vcf_convert.example.grg


.. parsed-literal::

    Will not count variants in VCF files (too slow)
    Could not count number of variants in vcf_convert.example.vcf.gz. Using the default of 100 (use --parts to override).
    Processing input file in 100 parts.
    Auto-calculating number of trees per part.
    Converting segments of input data to graphs
      0%|          | 0/100 [00:00<?, ?it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      1%|          | 1/100 [00:00<00:11,  8.59it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      3%|▎         | 3/100 [00:00<00:10,  9.21it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      4%|▍         | 4/100 [00:00<00:11,  8.41it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      5%|▌         | 5/100 [00:00<00:11,  8.09it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      6%|▌         | 6/100 [00:00<00:12,  7.48it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      8%|▊         | 8/100 [00:00<00:11,  8.01it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     10%|█         | 10/100 [00:01<00:10,  8.26it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     12%|█▏        | 12/100 [00:01<00:08, 10.44it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     14%|█▍        | 14/100 [00:01<00:06, 12.38it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     16%|█▌        | 16/100 [00:01<00:06, 13.96it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     18%|█▊        | 18/100 [00:01<00:05, 14.55it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     20%|██        | 20/100 [00:01<00:05, 15.77it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     22%|██▏       | 22/100 [00:01<00:04, 16.75it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     24%|██▍       | 24/100 [00:01<00:04, 16.86it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     26%|██▌       | 26/100 [00:02<00:04, 17.09it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     28%|██▊       | 28/100 [00:02<00:04, 16.53it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     30%|███       | 30/100 [00:02<00:04, 17.21it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     32%|███▏      | 32/100 [00:02<00:03, 17.71it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     34%|███▍      | 34/100 [00:02<00:03, 18.03it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     36%|███▌      | 36/100 [00:02<00:03, 17.70it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     38%|███▊      | 38/100 [00:02<00:03, 17.99it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     40%|████      | 40/100 [00:02<00:03, 18.33it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     42%|████▏     | 42/100 [00:02<00:03, 18.60it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     44%|████▍     | 44/100 [00:03<00:03, 18.28it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     46%|████▌     | 46/100 [00:03<00:02, 18.39it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     48%|████▊     | 48/100 [00:03<00:02, 18.11it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     50%|█████     | 50/100 [00:03<00:02, 18.22it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     52%|█████▏    | 52/100 [00:03<00:02, 18.52it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     54%|█████▍    | 54/100 [00:03<00:02, 18.66it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     56%|█████▌    | 56/100 [00:03<00:02, 18.71it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     58%|█████▊    | 58/100 [00:03<00:02, 18.85it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     60%|██████    | 60/100 [00:03<00:02, 19.00it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     62%|██████▏   | 62/100 [00:04<00:01, 19.06it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     64%|██████▍   | 64/100 [00:04<00:01, 18.94it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     66%|██████▌   | 66/100 [00:04<00:01, 18.87it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     68%|██████▊   | 68/100 [00:04<00:01, 18.93it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     70%|███████   | 70/100 [00:04<00:01, 19.03it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     72%|███████▏  | 72/100 [00:04<00:01, 19.07it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     74%|███████▍  | 74/100 [00:04<00:01, 19.14it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     76%|███████▌  | 76/100 [00:04<00:01, 18.74it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     78%|███████▊  | 78/100 [00:04<00:01, 17.23it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     80%|████████  | 80/100 [00:05<00:01, 16.92it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     82%|████████▏ | 82/100 [00:05<00:01, 16.11it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     84%|████████▍ | 84/100 [00:05<00:01, 15.69it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     86%|████████▌ | 86/100 [00:05<00:00, 15.05it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     88%|████████▊ | 88/100 [00:05<00:00, 14.98it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     90%|█████████ | 90/100 [00:05<00:00, 15.08it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     92%|█████████▏| 92/100 [00:05<00:00, 15.19it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     94%|█████████▍| 94/100 [00:05<00:00, 15.04it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     96%|█████████▌| 96/100 [00:06<00:00, 15.07it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     98%|█████████▊| 98/100 [00:06<00:00, 14.80it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    100%|██████████| 100/100 [00:06<00:00, 15.57it/s]
    Merging...


.. parsed-literal::

    === GRG Statistics ===
    Nodes: 12945
    Edges: 117253
    Samples: 400
    Mutations: 10893
    Ploidy: 2
    Phased: true
    Populations: 0
    Range of mutations: 55829 - 9999127
    Specified range: 0 - 100000001
    ======================
    Wrote simplified GRG with:
      Nodes: 12945
      Edges: 117253
    Wrote GRG to vcf_convert.example.grg


That gave us a valid GRG, by essentially “brute force” accessing the VCF
file without an index. This works fine on a small file, but not a large
one. So instead, lets try indexing this file with ``tabix`` and trying
against. **NOTE**: ``grg`` only supports ``tabix``-style indexes, and
not ``bcftools``-style indexes on ``.vcf.gz`` files.

.. code:: bash

    %%bash
    tabix vcf_convert.example.vcf.gz

Now we can construct a GRG without having to use ``--force``, and we
don’t get any warnings.

.. code:: bash

    %%bash
    
    # -j controls how many threads to use.
    grg construct -j 1 vcf_convert.example.vcf.gz -o vcf_convert.example.grg


.. parsed-literal::

    Will not count variants in VCF files (too slow)
    Could not count number of variants in vcf_convert.example.vcf.gz. Using the default of 100 (use --parts to override).
    Processing input file in 100 parts.
    Auto-calculating number of trees per part.
    Converting segments of input data to graphs
    100%|██████████| 100/100 [00:03<00:00, 30.38it/s]
    Merging...


.. parsed-literal::

    === GRG Statistics ===
    Nodes: 12945
    Edges: 117253
    Samples: 400
    Mutations: 10893
    Ploidy: 2
    Phased: true
    Populations: 0
    Range of mutations: 55829 - 9999127
    Specified range: 0 - 100000001
    ======================
    Wrote simplified GRG with:
      Nodes: 12945
      Edges: 117253
    Wrote GRG to vcf_convert.example.grg


Related Topics
--------------

-  Often, it is better to convert a ``vcf.gz`` to
   `IGD <https://picovcf.readthedocs.io/en/latest/igd_overview.html>`__
   first and then convert to GRG. IGD files can be substantially faster
   to access than ``.vcf.gz`` files. See `Converting IGD to
   GRG <IGDToGRG.html>`__.
