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

**What you’ll need:**

-  Python dependencies “pygrgl”: ``pip install pygrgl``
-  Command line tools “wget” and “tabix”:
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

    --2026-05-11 09:46:08--  https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving github.com (github.com)... 140.82.113.3
    Connecting to github.com (github.com)|140.82.113.3|:443... connected.
    HTTP request sent, awaiting response... 302 Found
    Location: https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz [following]
    --2026-05-11 09:46:08--  https://raw.githubusercontent.com/aprilweilab/grg_pheno_sim/refs/heads/main/demos/data/test-200-samples.vcf.gz
    Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 185.199.111.133, 185.199.109.133, 185.199.110.133, ...
    Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|185.199.111.133|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 494022 (482K) [application/octet-stream]
    Saving to: ‘vcf_convert.example.vcf.gz’
    
         0K .......... .......... .......... .......... .......... 10% 2.43M 0s
        50K .......... .......... .......... .......... .......... 20% 7.11M 0s
       100K .......... .......... .......... .......... .......... 31% 3.88M 0s
       150K .......... .......... .......... .......... .......... 41% 7.37M 0s
       200K .......... .......... .......... .......... .......... 51% 26.7M 0s
       250K .......... .......... .......... .......... .......... 62% 10.3M 0s
       300K .......... .......... .......... .......... .......... 72% 5.27M 0s
       350K .......... .......... .......... .......... .......... 82% 14.6M 0s
       400K .......... .......... .......... .......... .......... 93% 10.1M 0s
       450K .......... .......... .......... ..                   100% 7.97M=0.07s
    
    2026-05-11 09:46:09 (6.35 MB/s) - ‘vcf_convert.example.vcf.gz’ saved [494022/494022]
    


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
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 301, in star_build_grg
        return build_grg(*args)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 267, in build_grg
        shape_grg = build_shape(range_triple, args, auto_args, input_file, output_file)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 254, in build_shape
        tb_time = time_call(command, stdout=sys.stdout)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/common.py", line 48, in time_call
        subprocess.check_call(command, **kwargs)
      File "/usr/lib/python3.10/subprocess.py", line 369, in check_call
        raise CalledProcessError(retcode, cmd)
    subprocess.CalledProcessError: Command '['/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/../../grgl', 'vcf_convert.example.vcf.gz', '--trees', 'optimal', '--lf-no-tree', '10', '--reduce', '5', '-r', '0.0:0.01', '-o', 'vcf_convert.example.grg.part0.grg']' died with <Signals.SIGABRT: 6>.
    """
    
    The above exception was the direct cause of the following exception:
    
    Traceback (most recent call last):
      File "/home/ddehaas/Py3Env/bin/grg", line 6, in <module>
        sys.exit(main())
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/cli.py", line 70, in main
        construct.from_tabular(args)
      File "/home/ddehaas/GrgProject/public/grgl/pygrgl/clicmd/construct.py", line 422, in from_tabular
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
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      2%|▏         | 2/100 [00:00<00:07, 13.33it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      4%|▍         | 4/100 [00:00<00:09,  9.66it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      6%|▌         | 6/100 [00:00<00:11,  7.85it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      7%|▋         | 7/100 [00:00<00:11,  8.13it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      8%|▊         | 8/100 [00:00<00:12,  7.64it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
      9%|▉         | 9/100 [00:01<00:11,  8.09it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     10%|█         | 10/100 [00:01<00:12,  7.23it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     12%|█▏        | 12/100 [00:01<00:08,  9.82it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     14%|█▍        | 14/100 [00:01<00:07, 11.77it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     16%|█▌        | 16/100 [00:01<00:06, 13.41it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     18%|█▊        | 18/100 [00:01<00:05, 14.54it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     20%|██        | 20/100 [00:01<00:05, 15.41it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     22%|██▏       | 22/100 [00:01<00:04, 15.93it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     24%|██▍       | 24/100 [00:02<00:04, 15.59it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     26%|██▌       | 26/100 [00:02<00:04, 16.06it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     28%|██▊       | 28/100 [00:02<00:04, 16.36it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     30%|███       | 30/100 [00:02<00:04, 16.54it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     32%|███▏      | 32/100 [00:02<00:04, 16.83it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     34%|███▍      | 34/100 [00:02<00:04, 16.15it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     36%|███▌      | 36/100 [00:02<00:03, 16.42it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     38%|███▊      | 38/100 [00:02<00:03, 16.24it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     40%|████      | 40/100 [00:03<00:03, 16.32it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     42%|████▏     | 42/100 [00:03<00:03, 16.39it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     44%|████▍     | 44/100 [00:03<00:03, 16.48it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     46%|████▌     | 46/100 [00:03<00:03, 16.44it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     48%|████▊     | 48/100 [00:03<00:03, 16.14it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     50%|█████     | 50/100 [00:03<00:03, 16.19it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     52%|█████▏    | 52/100 [00:03<00:03, 15.94it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     54%|█████▍    | 54/100 [00:03<00:02, 16.04it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     56%|█████▌    | 56/100 [00:04<00:02, 16.14it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     58%|█████▊    | 58/100 [00:04<00:02, 16.16it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     60%|██████    | 60/100 [00:04<00:02, 15.78it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     62%|██████▏   | 62/100 [00:04<00:02, 16.20it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     64%|██████▍   | 64/100 [00:04<00:02, 16.56it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     66%|██████▌   | 66/100 [00:04<00:02, 16.91it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     68%|██████▊   | 68/100 [00:04<00:01, 17.43it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     70%|███████   | 70/100 [00:04<00:01, 17.04it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     72%|███████▏  | 72/100 [00:04<00:01, 16.72it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     74%|███████▍  | 74/100 [00:05<00:01, 16.84it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     76%|███████▌  | 76/100 [00:05<00:01, 16.46it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     78%|███████▊  | 78/100 [00:05<00:01, 16.35it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     80%|████████  | 80/100 [00:05<00:01, 17.08it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     82%|████████▏ | 82/100 [00:05<00:01, 17.30it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     84%|████████▍ | 84/100 [00:05<00:00, 17.32it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     86%|████████▌ | 86/100 [00:05<00:00, 17.48it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     88%|████████▊ | 88/100 [00:05<00:00, 17.09it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     90%|█████████ | 90/100 [00:06<00:00, 17.15it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     92%|█████████▏| 92/100 [00:06<00:00, 16.25it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     94%|█████████▍| 94/100 [00:06<00:00, 15.67it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     96%|█████████▌| 96/100 [00:06<00:00, 16.01it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
     98%|█████████▊| 98/100 [00:06<00:00, 16.36it/s]WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    WARNING: Conversion from VCF without a tabix index is very slow, and not recommended.
    100%|██████████| 100/100 [00:06<00:00, 15.01it/s]
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
    100%|██████████| 100/100 [00:03<00:00, 28.33it/s]
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
