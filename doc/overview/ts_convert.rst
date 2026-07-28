
.. _ts_convert:

Converting tree-sequence (TS) to GRG
------------------------------------

A tskit tree sequence (TS) is a particular format for Ancestral Recombination Graphs (ARGs).
ARGs, and tskit TS specifically, can be converted to GRGs. The sample-to-mutation (variant) mapping
is retained by this conversion, as well as the times associated with mutations, but there is some
information lost upon conversion to GRG:

- Tree sequences (sometimes implicitly) encode the location of recombinations along the genome.
  Converting to GRG will lose some of this information.
- Tree sequences have a node for each inferred/known coalescence between two lineages within the
  dataset. Conversion to GRG may drop such nodes if they do not contain any information
  about the mutation-to-sample relationship. It may combine many such nodes into a
  single node.

For real datasets, GRGs are usually much faster to infer than a tree sequence, see `GRG construction <construct.html>`_.
When simulated datasets get very large, calculations with GRG can be significantly faster than with
tree sequences. Conversion from TS to GRG is very fast (usually on the order of seconds or minutes).

The ``grg convert`` command will create a GRG from a TS. Example:

::

	grg convert /path/to/my.trees my.grg


The following properties of a tree sequence *are* maintained when converted to a GRG:

- The time associated with a Mutation is stored on the GRG :py:class:`pygrgl.Mutation` object.
- The sample IDs (NodeID associated with each sample) are identical between the TS and the GRG.
- All Mutations from the TS are copied into the GRG. Mutations in the TS can have no samples
  reachable beneath them (due to recurrent mutations "blocking" them), and these are copied into
  the GRG as a :py:class:`pygrgl.Mutation` with no associated node in the graph.

The following properties of a tree sequence are *not* maintained when converted to a GRG:

- The samples reached by a GRG *reflect the samples containing the variant*, and this can different
  from the TS topology in the presence of recurrent and back mutations. In the TS the topology reflects
  the coalescent tree, and in the GRG the topology reflects the actual sample sets for those mutations.
  For example, there could be 100 samples beneath a mutation ``A`` in a TS local tree, with a back
  mutation ``B`` beneath it that makes the "true" sample set only have 10 samples in it (``B`` acts
  like a set subtraction from ``A``). In the GRG, ``A`` will only have 10 samples beneath it.
- tskit has an ID number associated with each Mutation in the TS. GRG also has an ID number
  associated with each Mutation. These IDs *will not match* between the two formats. When comparing
  Mutations between the two formats, use the Mutation's position and allele values. For most 
  simulated data, using only the alternate allele and position are sufficient, but for some real
  datasets there are multiple different reference alleles at a site.
- Some nodes are dropped when converting to GRG, unless the ``--no-simplify`` option is provided (*NOT RECOMMENDED*).

Mapping between TS and GRG mutations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It can be useful to know which mutations are the same between a tskit TreeSequence and a GRG. The IDs
do not necessarily match after conversion! I.e., the tskit mutation ID is not the same as the GRG
MutationID. Below is illustrative code for how to perform the mapping, which should be based on the
position and allele values:

::

  # Map TS mutations to GRG mutations. This dictionary contains a mapping from (position, allele) to
  # the tskit mutation ID.
  pos_allele_to_tsid = {}
  for i, tree in enumerate(ts.trees()):
    for mut in tree.mutations():
      site = ts.site(mut.site)
      key = (site.position, mut.derived_state)
      # Map from position,allele to tskit ID
      pos_allele_to_tsid[key] = mut.id

  # Now lookup each GRG mutation to find the corresponding tskit ID
  for mut_id in range(grg.num_mutations):
      mut = grg.get_mutation_by_id(mut_id)
      tsid = pos_allele_to_tsid.get( (mut.position, mut.allele) )
      print(f"GRG ID = {mut_id}, tskit ID = {tsid}")


``MutationID`` to ``NodeID`` is a 1-to-1 relationship, and in *constructed* GRGs we only ever create
a single ``Mutation`` for every unique ``(position, allele, ref_allele)`` triple. However, when converting
from *tskit* you can have multiple ``Mutation``s for the same ``(position, allele, ref_allele)`` that represents
either the actual ancestry (two separate mutations, mostly in simulated ARGs) or is just an
artifact of the ARG inference process (e.g., *tsinfer*). In these cases, you need to adjust
any relevant calculations (you may *want* to keep separate ``Mutations``).


A note about Mutation times
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a TS, Mutations are mapped to *edges*, whereas in a GRG Mutations are mapped to *nodes*.
When converting a TS to GRG we map the nodes, so the edge information is potentially collapsed. This
means that Mutations which were on different edges in the TS may end up on the same node in the GRG.
In the presence of recombination you can also end up with a GRG child node containing a Mutation that is
older than a Mutation on the parent node (though rarely). The GRG conversion flag ``--use-node-times``
can avoid this Mutations-out-of-order problem by associating the coalesence time from the node below
the Mutation in the TS with the :py:class:`pygrgl.Mutation` in the GRG.


Node coalescence counts
~~~~~~~~~~~~~~~~~~~~~~~

`GWAS <examples_and_applications.html>`_ and zygosity information (see ``grg process zygosity --help``)
both require the GRG to have coalescence information at each node. This is just a count of the number of
diploid individuals that coalesced at any given node (both of their haploid samples are reachable from
the node, but not reachable together from any of the node's children). This information is automatically
calculated when you construct a GRG via ``grg construct``. In order to get this information in TS-converted
GRGs, you need to pass ``--ts-coals`` flags to the ``grg convert`` command.

Conversion from tree-sequence to GRG is very fast, on the order of seconds or minutes for very large datasets.
Setting flag ``--ts-coals`` can slow it down significantly, especially for datasets with a lot of samples.

You can still do most operations on a GRG without having the coalescence counts; even diploid-based GWAS
can be performed by using `grapp's <https://github.com/aprilweilab/grapp>`_ ``--binomial`` flag.

See also :py:meth:`pygrgl.GRG.get_num_individual_coals` and :py:meth:`pygrgl.GRG.set_num_individual_coals`.
