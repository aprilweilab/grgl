.. _examples:

Examples and applications using GRG
===================================

GWAS
----

The ``grg`` command implements GWAS using linear association affects for diploid data. The user provides
phenotype information for each diploid individual, which is a file containing a line per individual. The
order of the individuals must match the order of the individuals in the GRG (which *will* match the order
of the tabular file that the GRG was constructed from).

There are three space-separate columns for each line. Currently, the first two columns are ignored, and only
the third column (phenotype value) is used.

The output is emitted to stdout, so usually it is wise to redirect stdout to a file.

Example usage:

::

	grg process gwas test.grg --phenotype test.pheno > gwas_results.txt


Phenotype Simulation
--------------------

Given a GRG, you can simulate phenotypes for the individual samples contained within it. This is done using the
external package `grg_pheno_sim <https://github.com/aprilweilab/grg_pheno_sim/>`_, which can be installed via:

::

	pip install grg_pheno_sim

Usage examples can be found in the `example jupyter notebooks <https://github.com/aprilweilab/grg_pheno_sim/tree/main/demos>`_.


**TODO: Provide a simple end-to-end example that performs both phenotype simulation and GWAS.**

Splitting GRGs
--------------

See command ``grg split`` and Python API :py:meth:`save_subset`.

**TODO: Expand description.**