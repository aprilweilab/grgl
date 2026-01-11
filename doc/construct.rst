.. _construct:

Constructing GRGs from raw data
-------------------------------

GRGs are generally constructed from phased VCF, VCF.GZ, or IGD files. Unphased inputs of the same types
can also be used. While GRG is not yet optimized for unphased data, it still typically stores the data
in less space than most other techniques (though constructing the GRG will be slower than if the data
is phased).

IGD is the most efficient input format to use. When constructing GRGs from any non-trivial
dataset, avoid the use of unindexed VCF files.

Convert VCF to IGD
~~~~~~~~~~~~~~~~~~

A VCF (compressed or uncompressed) file can be converted to IGD using
`igdtools <https://picovcf.readthedocs.io/en/latest/igdtools.html>`_ (``pip install igdtools``).
If the VCF is a BGZF file with a `Tabix <https://www.htslib.org/doc/tabix.html>`_ index,
then you can use the ``-j <threads>`` flag to significantly speedup conversion.

Example1:

::

	igdtools path/to/foo.vcf -o foo.igd

Example2:

::

	igdtools -j 20 path/to/indexed.vcf.gz -o foo.igd


If further manipulation of the IGD file is needed, you can use `picovcf
<https://github.com/aprilweilab/picovcf>`_ (C++) or `pyigd <https://github.com/aprilweilab/pyigd>`_ (Python).


Construct GRG from IGD 
~~~~~~~~~~~~~~~~~~~~~~

The key parameter for GRG construction is ``-j``, how many cores/threads to use during construction.
See ``grg construct --help`` for other parameters, though most users will not need them.

Example for constructing a GRG from 1,000 Genomes chromosome 22:

::

	grg construct -j 6 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.v3.igd

On a typical machine (e.g., my laptop) the above takes about 30 seconds. To give an example of a
more large-scale dataset, constructing chromosome 22 of the phased UK Biobank dataset (200,000 individuals)
takes about 12 minutes (using ``-j 70``).

Construct GRG from indexed VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the VCF is BGZF compressed with a `Tabix <https://www.htslib.org/doc/tabix.html>`_ index,
then GRG construction will be reasonably fast. For datasets with a small-to-medium number of
samples (e.g., 10,000 or less) construction from indexed VCF is not much slower than from IGD.
However, for really large sample sizes the GRG construction can be many times slower than from IGD.

Example for constructing a GRG from 1,000 Genomes chromosome 22:

::

	grg construct -j 6 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


Construct GRG from unindexed VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

THIS IS NOT RECOMMENDED! This functionality primarily exists for testing purposes, as building a GRG
from an unindexed VCF will be absurdly slow for a large dataset. A GRG is constructed by partitioning
a dataset along the chromosome/genome and performing many small graph constructions, which are then 
merged into a single graph -- this cannot be done efficiently on unindexed VCF.

For this case, you must specify the ``--force`` option to override warnings about the performance:

::

	grg construct -j 10 --force my_test.vcf


Construct GRG via Python API
----------------------------

A GRG can be constructed using arbitrary input data, by making use of the Python API.
See the methods on :py:class:`pygrgl.MutableGRG` which can be used to create nodes and connect
them, as well as merge two or more GRGs.

Extra Options
-------------

Run ``grg construct -h`` and ``grg convert -h`` to see a full list of GRG construction options. We highlight a few here.

Population Labels
~~~~~~~~~~~~~~~~~

GRG stores population labels against sample nodes. These can be retrieved from a GRG in Python using
:py:meth:`pygrgl.GRG.get_population_id` (passing the sample node ID). The population ID can then be
looked up in list of population descriptions (via :py:meth:`pygrgl.GRG.get_populations`).

When converting from a tskit tree sequence, populations are automatically attached to the resulting
GRG. When constructing a GRG from genotype data, you have to specify the ``--population-ids`` flag.
Population data is specified via a tab-separated input text file. Here is an example, ``pop.tsv``:

::

	SAMPLE_NAME POP_NAME
	samp1    CEU
	samp99   CHB
	samp4    YRI

Assume the whitespace is tab characters. Then the GRG would be constructed with:

::

	grg construct input.vcf.gz --population-ids pop.tsv:SAMPLE_NAME:POP_NAME

The syntax of the ``--population-ids`` flag is ``filename:field1:field2`` where ``field1`` is the
column name of the sample identifier column. ``field2`` is the column name for the population identifier.