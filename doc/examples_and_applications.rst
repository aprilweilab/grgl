.. _examples:

Examples and applications using GRG
===================================

GWAS
----

The ``grapp`` tool supports GWAS using linear association affects for diploid data, using a GRG. install
via

::

  pip install grapp


See

::

	grapp --help

for details.


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

See command ``grg split --help`` and Python API :py:meth:`save_subset`.

Splitting equally
~~~~~~~~~~~~~~~~~

You can use ``grg split input.grg -s <size in base-pair>`` to split a GRG into equal graphs, each of which cover
a base-pair range as specified by the ``-s`` flag. If you specify the ``--rec-map <hapmap-style file>`` option,
the size is assumed to be in centimorgans (cM).

If you have a list of ranges that you want to split a GRG into, you can put them in a text file with two columns
(space separated) and a header line. Here is an example file:

::

  start end
  0     1000000
  1000000 2000000
  5000000 7000000

Then save that file (e.g., as ``ranges.txt``) and pass it to the split command like: ``grg split input.grg -f ranges.txt``.
This will create 3 GRG files, each spanning one of the ranges from the text file. These ranges are inclusive on "start"
and exlusive on "end", so to get full coverage of the genome you want each consecutive range to have a "start" that is
identical to the previous range's "end".

By default, for a file ``input.grg`` the split command will produce a directory ``input.grg.split/`` and place all the
resulting split GRGs in that directory. You can change the directory using the ``-o`` flag, like:
``grg split input.grg -s 1000000 -o my_split_grgs``.

The output directory must not already exist, or splitting will fail and ask you to remove the directory or
specify a different one.