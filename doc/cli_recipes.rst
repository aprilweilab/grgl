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


Split GRG into pieces
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # Each resulting GRG will be a 5MB of the genome, with all the samples
    grg split -s 5000000 my_input.grg


Perform GWAS
~~~~~~~~~~~~

.. code:: bash

    grg assoc -p my_phenotype.phen my_input.grg -o my_input.assoc.tsv


Perform PCA
~~~~~~~~~~~

.. code:: bash

    # Get the top 20 PCs
    grg pca -d 20 my_input.grg


Filter by list of individual IDs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # list_of_individuals.txt has one ID per line
    grapp filter -S list_of_individuals.txt my_input.grg my_input.filtered.grg


Show the individual IDs
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    grapp show -S my_input.grg

GRG to IGD
~~~~~~~~~~

.. code:: bash

    # This can be slow! Use more threads (-j) if possible
    grapp export -j 4 my_input.grg --igd exported.igd

