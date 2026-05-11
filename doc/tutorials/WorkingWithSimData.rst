Working with Simulated Data
===========================

Simulated data from
`msprime <https://tskit.dev/msprime/docs/stable/intro.html>`__ or
`SLiM <https://messerlab.org/slim/>`__ can be saved in the `tskit
TreeSequence <https://tskit.dev/tskit/docs/stable/introduction.html>`__
format, which can be directly converted to GRG. It is worth noting that
GRG’s converted from ``tskit.TreeSequence`` do not have the same
performance characteristics as GRG’s converted from the corresponding
tabular (e.g., ``.vcf.gz`` or ``.igd``) data (a consideration when doing
benchmarking).

We’ll use the ``msprime`` engine from
`stdpopsim <https://popsim-consortium.github.io/stdpopsim-docs/stable/index.html>`__
for this tutorial.

**What you’ll need:**

-  Python dependencies “stdpopsim”, “msprime”, “grapp”:
   ``pip install stdpopsim msprime grapp``

Simulate Data
-------------

First lets simulate a small example dataset.

.. code:: ipython3

    import stdpopsim
    
    # We will simulate 200 diploid individuals in total, equally spread across populations.
    individuals = 200
    i_per_pop = individuals // 2
    
    # Setup the options for stdpopsim
    species = stdpopsim.get_species("HomSap")
    samples = {"AFR": i_per_pop, "EUR": i_per_pop}
    contig = species.get_contig("chr22", genetic_map="HapMapII_GRCh38")
    engine = stdpopsim.get_engine("msprime")
    model = species.get_demographic_model("OutOfAfrica_2T12")
    
    # Simulate the data, which returns an ARG "ts" and then we also save that ARG to disk.
    ts = engine.simulate(model, contig, samples)
    ts_filename = "stdpop.ooa2.chr22.trees"
    ts.dump(ts_filename)


.. parsed-literal::

    /home/ddehaas/Py3Env/lib/python3.10/site-packages/stdpopsim/engines.py:111: UserWarning: The demographic model has mutation rate 2.36e-08, but this simulation used the contig's mutation rate 1.29e-08. Diversity levels may be different than expected for this species. For details see documentation at https://popsim-consortium.github.io/stdpopsim-docs/stable/tutorial.html
      warnings.warn(


.. code:: ipython3

    # Count the recurrent and back mutations in the TreeSequence
    recurrent = 0
    back = 0
    already_seen = set()
    for s in ts.sites():
        for m in s.mutations:
            if m.parent != -1:
                back += 1
            key = (s.position, m.derived_state)
            if key in already_seen:
                recurrent += 1
            already_seen.add(key)

.. code:: ipython3

    print("## TreeSequence Info ##")
    print(f"Mutations: {ts.num_mutations}")
    print(f"Haplotypes: {ts.num_samples}")
    print(f"Nodes: {ts.num_nodes}")
    print(f"Edges: {ts.num_edges}")
    print(f"Recurrent Mutations: {recurrent}")
    print(f"Back Mutations: {back}")


.. parsed-literal::

    ## TreeSequence Info ##
    Mutations: 175356
    Haplotypes: 400
    Nodes: 151705
    Edges: 1021127
    Recurrent Mutations: 128
    Back Mutations: 51


Convert to GRG
--------------

Converting directly to GRG can be much faster than exporting the
``tskit.TreeSequence`` to VCF and then building a GRG. We illustrate
both the command-line and API methods for this conversion.

.. code:: bash

    %%bash
    # You can also use pygrgl.grg_from_trees(), but `grg convert` produces a simplified graph.
    grg convert stdpop.ooa2.chr22.trees stdpop.ooa2.chr22.trees.grg


.. parsed-literal::

    Warning: no individual (diploid) coalescence information will be calculated unless you specify --ts-coals. This is fine, but downstream applications will need to use approximations for Mutation variance.


.. code:: ipython3

    import pygrgl
    
    # Load the GRG we converted above
    grg = pygrgl.load_immutable_grg(ts_filename + ".grg")
    
    # Convert the TS to GRG using the API - just for comparison
    api_grg = pygrgl.grg_from_trees(ts_filename)

.. code:: ipython3

    print("## GRG Info ##")
    print(f"Mutations: {grg.num_mutations}")
    print(f"Haplotypes: {grg.num_samples}")
    print(f"Nodes: {grg.num_nodes}")
    print(f"Edges: {grg.num_edges}")


.. parsed-literal::

    ## GRG Info ##
    Mutations: 175208
    Haplotypes: 400
    Nodes: 197416
    Edges: 971022


As you can see, the GRG is not just an exact copy of the TreeSequence:
\* There are fewer Mutations: this is due to recurrent and back
mutations that we counted above. GRG tries to keep only a single node
per unique variant (position and alleles): both ``grg construct`` and
``grg convert`` (API: ``pygrgl.grg_from_trees``) create a single node
(and Mutation object) for each unique variant. Some back mutations are
no effect and are dropped by GRG: e.g., when two mutations for the same
site occur on the same ARG branch (i.e. they share a node). \* For this
ARG, the TreeSequence has fewer nodes. As the datasets get large, this
becomes increasingly unlikely. A GRG is a multitree that is created by
duplicating nodes in the coalescent trees - i.e., a node in
``tskit.TreeSequence`` can have different sets of samples beneath it in
different coalescent trees, but in GRG each node has a single sample
set. However, GRG also *removes* nodes from the ARG when that node does
not provide compression or information, hence in large datasets there
are typically fewer nodes.

Now if we look at the GRG that was converted via the API method, we see
there are substantially more nodes and edges. This is because GRG
simplifies the graph *during serialization to the disk*. So when you
convert via API you are getting the unsimplified graph - it has not yet
been written to disk.

.. code:: ipython3

    print("## GRG Info ##")
    print(f"Mutations: {api_grg.num_mutations}")
    print(f"Haplotypes: {api_grg.num_samples}")
    print(f"Nodes: {api_grg.num_nodes}")
    print(f"Edges: {api_grg.num_edges}")


.. parsed-literal::

    ## GRG Info ##
    Mutations: 175208
    Haplotypes: 400
    Nodes: 730624
    Edges: 1460403


Converting to IGD
-----------------

Now lets say we want a non-graph representation of the data. We could
convert the ARG to VCF with the ``tskit vcf`` command, but (a) VCF is
large and slow to work with and (b) converting to VCF can also be really
slow! Instead, now that we have a GRG we can convert it to
`IGD <https://github.com/aprilweilab/picovcf#indexable-genotype-data-igd>`__
and then we can use
`pyigd <https://pyigd.readthedocs.io/en/latest/igd_overview.html>`__ to
examine the data.

.. code:: ipython3

    from grapp.util.igd import export_igd
    
    # Convert to IGD using 4 threads
    export_igd(grg, ts_filename + ".igd", jobs=4, verbose=True)


.. parsed-literal::

    Using temporary directory /tmp/tmpsbd8vq95.
    Splitting GRG into 8 parts..
    Converting GRG parts to IGD files...
    Merging 8 parts into single IGD stdpop.ooa2.chr22.trees.igd...


.. code:: ipython3

    import pyigd
    
    with open(ts_filename + ".igd", "rb") as f:
        igd = pyigd.IGDReader(f)
    
        print("## IGD Info ##")
        print(f"Variants: {igd.num_variants}")
        print(f"Haplotypes: {igd.num_samples}")


.. parsed-literal::

    ## IGD Info ##
    Variants: 175208
    Haplotypes: 400


Related Topics
--------------

-  See `IGDToGRG <IGDToGRG.html>`__ for information on constructing GRGs
   from IGD data
-  Documentation links:

   -  `grapp.util <https://grapp.readthedocs.io/en/latest/grapp.html#filtering-export-etc>`__:
      Operations on GRG related to filtering, export, etc.

