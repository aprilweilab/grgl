Command-line Recipes
====================

A bunch of command-line examples provided without much explanation.
See the Tutorials or other documentation for full explanations.

.vcf.gz to GRG
~~~~~~~~~~~~~~

.. code:: bash

    # Make sure my_input.vcf.gz is tabix indexed!
    grg construct -j 6 my_input.vcf.gz -o my_input.grg


.vcf.gz to IGD
~~~~~~~~~~~~~~

.. code:: bash

    # Make sure my_input.vcf.gz is tabix indexed!
    igdtools -j 6 my_input.vcf.gz -o my_input.igd


IGD to GRG
~~~~~~~~~~

.. code:: bash

    # Will always be faster than the .vcf.gz conversion
    grg construct -j 6 my_input.igd -o my_input.grg

View info about GRG
~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grg process stats my_input.grg

.. code:: bash

    grapp show -i my_input.grg


Split GRG into pieces
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # Each resulting GRG will be a 5MB of the genome, with all the samples
    grg split -s 5000000 my_input.grg


Perform GWAS
~~~~~~~~~~~~

.. code:: bash

    grapp assoc -p my_phenotype.phen my_input.grg -o my_input.assoc.tsv


Perform PCA
~~~~~~~~~~~

.. code:: bash

    # Get the top 20 PCs
    grapp pca -d 20 my_input.grg


Filter by list of individual IDs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # list_of_individuals.txt has one ID per line
    grapp filter -S list_of_individuals.txt my_input.grg my_input.filtered.grg

Filter bi-allelic SNPs with frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequency :math:`\ge 0.01`, bi-allelic SNPs only. Command flags are similar to ``bcftools``.

.. code:: bash

    # list_of_individuals.txt has one ID per line
    grapp filter -q 0.01 -v snps -m 2 -M 2 my_input.grg my_input.filtered.grg


Show the individual IDs
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grapp show -S my_input.grg

Show variants with allele frequencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grapp show -c my_input.grg

Show variants with HWE p-values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grapp show -j 4 -H my_input.grg


Simulate phenotypes with h^2=0.4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grapp pheno -e 0.4 my_input.grg


GRG to IGD
~~~~~~~~~~

.. code:: bash

    # This can be slow! Use more threads (-j) if possible
    grapp export -j 4 my_input.grg --igd exported.igd
