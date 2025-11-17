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
 package `grg_pheno_sim <https://github.com/aprilweilab/grg_pheno_sim/>`_, which can be installed via::

	pip install grg_pheno_sim

Basic Example
~~~~~~~~~~~~~

Below is a minimal Python example using the `sim_phenotypes` function:

.. code-block:: python

    import pygrgl
    from grg_pheno_sim.phenotype import sim_phenotypes

    # 1. Load an immutable GRG
    grg = pygrgl.load_immutable_grg("test.grg")

    # 2. Set desired heritability (h²)
    heritability = 0.33

    # 3. Simulate phenotypes
    phenotypes = sim_phenotypes(grg, heritability=heritability)

The returned DataFrame contains one row per individual and has the following columns:

- ``causal_mutation_id`` – the identifier for each causal mutation (for multivariate simulations)  
- ``individual_id`` – the external individual identifier  
- ``genetic_value`` – the sum of all causal-mutation effects for the individual  
- ``environmental_noise`` – the simulated environmental component  
- ``phenotype`` – the total phenotype value  

Extended Example with additional parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sim_phenotypes`` function supports a range of optional parameters for customizing 
the simulation process. This example shows how to specify the effect-size model, 
control the number of causal variants, set reproducibility parameters, normalize phenotypes, 
and save intermediate outputs.

.. code-block:: python

    import pygrgl
    from grg_pheno_sim.phenotype import sim_phenotypes
    from grg_pheno_sim.model import grg_causal_mutation_model

    # 1. Load an immutable GRG
    grg_1 = pygrgl.load_immutable_grg("test.grg")

    # 2. Define effect-size model
    model_type = "normal"
    mean = 0
    var = 1
    model = grg_causal_mutation_model(model_type, mean=mean, var=var)

    # 3. Simulation parameters
    num_causal = 1000
    random_seed = 1
    normalize_phenotype = True                     # normalize final phenotypes
    normalize_genetic_values_before_noise = True   # normalize genetic values before adding noise
    heritability = 0.33
    effect_output_required = True                  # save effect sizes to .par file
    effect_path = 'univariate_sample_effect_sizes.par'
    standardized_output = True
    output_path = 'normal_pheno_normalized.phen'   # output phenotype file
    header = True                                  # include header row in output

    # 4. Run simulation
    phenotypes = sim_phenotypes(
        grg_1,
        model,
        num_causal,
        random_seed,
        normalize_phenotype=normalize_phenotype,
        normalize_genetic_values_before_noise=normalize_genetic_values_before_noise,
        heritability=heritability,
        save_effect_output=effect_output_required,
        effect_path=effect_path,
        standardized_output=standardized_output,
        path=output_path,
        header=header
    )

This produces a phenotype DataFrame (and optional saved files) based on the chosen model 
and parameters.

Inspecting output files
~~~~~~~~~~~~~~~~~~~~~~~

When ``save_effect_output=True`` is specified an effect sizes .par file is produced and when ``standardized_output=True`` is specified a phenotype file :

1. **Effect sizes file** – contains effect sizes for each mutation in the GRG.  
   Example from ``univariate_sample_effect_sizes.par``:

   .. list-table::
      :header-rows: 1
      :widths: 15 15 15 15 15 15

      * - mutation_id
        - AlternateAllele
        - Position
        - RefAllele
        - Frequency
        - Effect
      * - 20
        - A
        - 73883
        - C
        - 0.000
        - -1.810258
      * - 28
        - G
        - 81128
        - A
        - 0.335
        - 1.151768
      * - 62
        - G
        - 117367
        - C
        - 0.840
        - 1.681257
      * - 76
        - A
        - 134028
        - C
        - 0.010
        - 2.346698
      * - 119
        - A
        - 180670
        - T
        - 0.840
        - -0.286668

   - ``mutation_id`` – unique ID for each mutation in the GRG  
   - ``AlternateAllele`` / ``RefAllele`` – allele information  
   - ``Position`` – genomic position  
   - ``Frequency`` – allele frequency in the dataset  
   - ``Effect`` – effect size used in phenotype simulation  

2. **Phenotype file** – contains the simulated phenotypes for each individual.  
   Example from ``normal_pheno_normalized.phen``:

   .. list-table::
      :header-rows: 1
      :widths: 20 20

      * - person_id
        - phenotype
      * - 0
        - 0.715620
      * - 1
        - -0.706513
      * - 2
        - 0.416523
      * - 3
        - 0.657675
      * - 4
        - -0.747519

   - ``person_id`` – numeric identifier for each individual  
   - ``phenotype`` – simulated phenotype value  




Binary phenotype (case–control) example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supports binary phenotype simulation.  
This step follows the standard pipeline to simulate a continuous phenotype,  
then employs a ``population_prevalence`` parameter to establish a Gaussian threshold  
to convert continuous values into binary outcomes.  
Lastly, the system overwrites continuous phenotypes with their binary counterparts.

.. code-block:: python

    import pygrgl
    from grg_pheno_sim.phenotype import sim_binary_phenotypes

    # 1. Load an immutable GRG
    grg_1 = pygrgl.load_immutable_grg("test.grg")

    # 2. Set model parameters
    heritability = 0.33                  # liability-scale h^2
    population_prevalence = 0.10         # 1 in 10 individuals are cases on average

    # 3. Simulate binary phenotypes (0 = control, 1 = case)
    phenotypes_binary = sim_binary_phenotypes(
        grg_1,
        population_prevalence,
        heritability=heritability,
        # optional:
        # num_causal=1000,
        # random_seed=42,
        # standardized=True,
    )

The output is the same DataFrame as the previous example except ``phenotype`` is now a binary outcome (``0`` control, ``1`` case) determined by the prevalence threshold.


Standardized genotype example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sim_phenotypes`` function also supports internally standardizing the genotype matrix 
before computing genetic values. This is done by setting ``standardized=True``, which 
centers and scales each mutation's genotype values to have mean zero and variance one.

.. code-block:: python

    import pygrgl
    from grg_pheno_sim.phenotype import sim_phenotypes

    # 1. Load an immutable GRG
    grg_1 = pygrgl.load_immutable_grg("test.grg")

    # 2. Set desired heritability (h²)
    heritability = 0.33

    # 3. Simulate phenotypes with standardized genotypes
    phenotypes_standardized_genes = sim_phenotypes(
        grg_1,
        heritability=heritability,
        standardized=True
    )

This produces the same output format as the basic example, but with genotypes standardized before effect application.


Custom effect sizes example
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sim_phenotypes_custom`` function allows you to provide your own effect sizes instead 
of sampling from a distribution-based model. This can be useful for testing, reproducing known 
phenotypes, or running simulations with fixed effect patterns.

Supported effect size input types:

- **List** – a list of effect sizes in mutation order  
- **Dictionary** – keys are mutation IDs, values are effect sizes  
- **Pandas DataFrame** – must contain at least ``mutation_id`` and ``effect_size`` columns  
  (for multivariate simulations, also include ``causal_mutation_id``)  

.. code-block:: python

    import numpy as np
    import pandas as pd
    import pygrgl
    import random
    from grg_pheno_sim.phenotype import sim_phenotypes_custom

    # 1. Load an immutable GRG
    grg_1 = pygrgl.load_immutable_grg("test.grg")

    # 2. Prepare effect sizes
    n = grg_1.num_mutations
    specific_effects = [1.0 for _ in range(n)]  # list input, fixed effects
    effect_sizes = np.random.randn(n)           # NumPy array
    mutation_dict = {i: effect_sizes[i] for i in range(n)}  # dictionary input
    input_df = pd.DataFrame(list(mutation_dict.items()), columns=['mutation_id', 'effect_size'])  # DataFrame input

    # 3. Simulation parameters
    heritability = 0.33
    normalize_genetic_values_before_noise = True
    standardized_output = True
    output_path = 'custom_pheno.phen'

    # 4. Simulate phenotypes using custom effect sizes
    phenotypes_list = sim_phenotypes_custom(
        grg_1,
        specific_effects,
        normalize_genetic_values_before_noise=normalize_genetic_values_before_noise,
        heritability=heritability,
        standardized_output=standardized_output,
        path=output_path
    )

The returned DataFrame has the same structure as in the basic example, with phenotypes generated using 
your provided effect sizes rather than randomly sampled ones.


Quantile normalization example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates quantile normalization for a single causal mutation's phenotypic DataFrame.

.. code-block:: python

    import numpy as np
    import pygrgl
    import matplotlib.pyplot as plt
    from grg_pheno_sim.effect_size import sim_grg_causal_mutation, additive_effect_sizes, samples_to_individuals
    from grg_pheno_sim.model import grg_causal_mutation_model
    from grg_pheno_sim.noise_sim import sim_env_noise
    from grg_pheno_sim.normalization import quantile_normalize

    # 1. Load an immutable GRG
    grg_1 = pygrgl.load_immutable_grg("test.grg")

    # 2. Define model for causal effects
    mean_1 = 0.0
    var_1 = 1.0
    model_normal = grg_causal_mutation_model("normal", mean=mean_1, var=var_1)

    # 3. Simulate genetic values
    trait_df_normal = sim_grg_causal_mutation(grg_1, num_causal=1000, model=model_normal, random_seed=1)
    sample_nodes_df = additive_effect_sizes(grg_1, trait_df_normal)
    individual_genetic_value_df = samples_to_individuals(sample_nodes_df)

    # 4. Add environmental noise without normalizing genetic values
    phenotypes = sim_env_noise(individual_genetic_value_df, h2=0.5)
    phenotype_df = phenotypes.phenotype_df

    # 5. Quantile normalize to the normal distribution
    quantile_normalize_phenotype_df = quantile_normalize(phenotype_df)

More usage examples can be found in the `example jupyter notebooks <https://github.com/aprilweilab/grg_pheno_sim/tree/main/demos>`_.


**TODO: Provide a simple end-to-end example that performs both phenotype simulation and GWAS.**

Splitting equally
-----------------

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