
.. _ts_convert:

Converting tree-sequence (TS) to GRG
------------------------------------

A tskit tree sequence (TS) is a particular format for Ancestral Recombination Graphs (ARGs).
ARGs, and tskit TS specifically, can be converted to GRGs. This is a lossy
conversion:

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

- tskit has an ID number associated with each Mutation in the TS. GRG also has an ID number
  associated with each Mutation. These IDs *will not match* between the two formats. When comparing
  Mutations between the two formats, use the Mutation's position and allele values. For most 
  simulated data, using only the alternate allele and position are sufficient, but for some real
  datasets there are multiple different reference alleles at a site.
- Some nodes are dropped when converting to GRG, unless the ``--no-simplify`` option is provided.

A note about Mutation times
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a TS, Mutations are mapped to *edges*, whereas in a GRG Mutations are mapped to *nodes*.
When converting a TS to GRG we map the nodes, so the edge information is potentially collapsed. This
means that Mutations which were on different edges in the TS may end up on the same node in the GRG.
In the presence of recombination you can also end up with a GRG child node containing a Mutation that is
older than a Mutation on the parent node (though rarely). The GRG conversion flag ``--use-node-times``
can avoid this Mutations-out-of-order problem by associating the coalesence time from the node below
the Mutation in the TS with the :py:class:`pygrgl.Mutation` in the GRG.
