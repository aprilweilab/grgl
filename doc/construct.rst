.. _construct:

Constructing GRGs from raw data
-------------------------------

GRGs can be constructed from phased VCF, VCF.GZ, or IGD files. Phased BGEN is also supported,
but only on Linux and it requires installing GRGL from source (see the 
`installation instructions <installation.html>`_).

IGD is the most efficient input format to use. When constructing GRGs from any non-trivial
dataset, use IGD instead of VCF(.GZ).

Convert VCF to IGD
~~~~~~~~~~~~~~~~~~

A VCF (compressed or uncompressed) file can be converted to IGD using the ``grg convert``
command:

::

	grg convert path/to/foo.vcf foo.igd

If further manipulation of the IGD file is needed, you can use `igdtools
<https://github.com/aprilweilab/picovcf>`_ or `pyigd <https://github.com/aprilweilab/pyigd>`_.


Construct GRG from IGD 
~~~~~~~~~~~~~~~~~~~~~~

There are three key parameters for GRG construction:
  1. ``-p``: How many segments to split the input genome (typically a single
  chromosome) into. If ``-p`` is too small, GRG construction will likely take
  a long time. If ``-p`` is too large, the resulting GRG file might be large.
  2. ``-t``: How many trees *within* each segment to create. This affects the shape
  of the graph that Mutations are mapped to. The largest ``-t`` value we have
  ever used is ``32``, and most datasets will want a value less than ``10``.
  3. ``-j``: How many cores/threads to use during construction. This should be no
  greater than the ``-p`` value.

Guidance:
- If your data is similar to the 1,000 Genomes dataset (high diversity but
few samples) you likely want a small ``-p`` value and a large ``-t`` value: something
like ``-p 50 -t 16``.
- If your data is simulated, you probably want a larger ``-p`` and smaller ``-t`` value.
We typically use something like ``-p 100 -t 5`` for data simulated with a non-trivial
(e.g., out-of-Africa) demographic model and just ``-p 100`` or ``-p 50`` for panmictic
data.
- Large, real datasets typically requires a large ``-p`` and small ``-t``. Making ``-p``
larger can improve performance when constructing GRGs, so even something as large
as ``-p 400 -t 4`` can work well on biobank-scale dataset.
- You can always construct a GRG for a small region of the genome to try to experimentally
determine the best ``-p`` and ``-t`` values. See ``grg construct --range``. If you have A
chromosome of length ``L`` (base-pairs) and theorize that ``-p P`` is a good value, then
the region size to test is ``s = L/P``. You can use ``grg construct --range`` to create
GRGs of size ``s``, ``s - delta`` and ``s + delta`` to see which is smallest and what
the time-to-construct is. 

Example for constructing a GRG from 1,000 Genomes chromosome 22:

::

	grg construct -p 50 -j 20 -t 16 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.v3.igd


Another flag that can affect the resulting GRG size is ``--maf-flip``. By default, a GRG
*implicitly* encodes the reference allele for each site and explicitly encodes
each alternate allele as a :py:class:`pygrgl.Mutation` attached to a graph node. Each of these
Mutation objects will list the reference allele. In contrast, a GRG built with 
``--maf-flip`` will *change* the reference to be the major allele (the one with the
highest frequency). This means the information about which allele was originally the
reference is lost, but the graph is faster to construct and smaller in size (because the
alleles that are explicitly represented as Mutations have fewer samples associated with
them).

See ``grg construct --help`` for more options that affect GRG construction.


Construct GRG via Python API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A GRG can be constructed using arbitrary input data, by making use of the Python API.
See the methods on :py:class:`pygrgl.MutableGRG` which can be used to create nodes and connect
them, as well as merge two GRGs.
