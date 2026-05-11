GWAS Tutorial
=============

This tutorial describes how to perform a GWAS with simple linear
regression, where each SNP is treated independently when determining its
effect on the phenotype. Lets consider our input data to be a genotype
matrix :math:`X` which has :math:`N` rows (number of individuals) and
:math:`M` columns (number of variants, usually SNPs). We’re using
`Genotype Representation
Graphs <https://grgl.readthedocs.io/en/stable/concepts.html#grg>`__ to
represent our genotype matrix. If you have ``.vcf.gz`` data, see
`Converting .vcf.gz to GRG <VCFToGRG.html>`__. If you have (or want) a
more efficicient representation than ``.vcf.gz`` see `Converting IGD to
GRG <IGDToGRG.html>`__.

For this tutorial, we’ll just use a (really) small simulated dataset of
200 individuals that we can quickly download. Replace it with your own
dataset as desired.

The other piece of information we need for performing a GWAS is the
phenotype, :math:`Y`. We are just going to simulate a phenotype for this
tutorial, but replace it with your own phenotype as desired.

**What you’ll need:**

-  Python dependencies “grapp”, “igdtools”, “seaborn”:
   ``pip install grapp igdtools seaborn``
-  Command line tool “wget”: ``sudo apt install wget`` (or your system’s
   equivalent)

Get Dataset
~~~~~~~~~~~

.. code:: bash

    %%bash
    
    if [[ ! -e gwas.example.igd ]]; then
      # Download a small example dataset
      wget https://github.com/aprilweilab/grg_pheno_sim/raw/refs/heads/main/demos/data/test-200-samples.vcf.gz -O gwas.example.vcf.gz
    
      # Convert to IGD; this isn't necessary, but most of the time you will want to do this
      igdtools gwas.example.vcf.gz -o gwas.example.igd
    fi
    
    # Just show some stats about the dataset
    igdtools -s gwas.example.igd


.. parsed-literal::

    Stats for gwas.example.igd
    ... in range 0 - 18446744073709551615
      Variants in range: 10893
      Average samples/var: 50.2075
      Stddev samples/var: 86.1985
      Average var/sample: 1367.28
      Stddev var/sample: 25.9211
      Variants with missing data: 0
      Total missing alleles: 0
      Total unique sites: 10885


Convert to GRG
--------------

.. code:: bash

    %%bash
    
    if [[ ! -e gwas.example.grg ]]; then
      # -j controls how many threads to use.
      grg construct -j 1 gwas.example.igd -o gwas.example.grg
    fi

Simulate Phenotype
------------------

Since this phenotype is just for demonstration purposes, we can use the
default settings for phenotype simulation. See `Simulating
Phenotypes <SimulatingPhenotypes.html>`__ for a more detailed tutorial
that describes the range of options. Here, we just want a phenotype that
has been derived from the genotype: that is, our simulation will choose
some “causal mutations” that are responsible for some proportion of the
phenotype (this is called the heritability), and compute for each
individual a phenotype value from the actual mutations the individual
has.

By default, phenotype simulation draws effect sizes for causal SNPs from
a standard normal distribution ~\ :math:`N(0, 1)`.

.. code:: bash

    %%bash
    
    grapp pheno gwas.example.grg


.. parsed-literal::

    The initial effect sizes are 
           mutation_id  effect_size  causal_mutation_id
    0                0    -1.657757                   0
    1                1    -0.537562                   0
    2                2     1.146189                   0
    3                3    -0.697705                   0
    4                4    -0.157002                   0
    ...            ...          ...                 ...
    10888        10888    -0.217740                   0
    10889        10889     0.936594                   0
    10890        10890     0.976964                   0
    10891        10891     0.274724                   0
    10892        10892     0.711005                   0
    
    [10893 rows x 3 columns]
    The genetic values of the individuals are 
         individual_id  genetic_value  causal_mutation_id
    0                0     -62.758200                   0
    1                1       8.271010                   0
    2                2     -93.726026                   0
    3                3     -58.657459                   0
    4                4     -47.272752                   0
    ..             ...            ...                 ...
    195            195     -22.379358                   0
    196            196      34.579149                   0
    197            197     -40.854451                   0
    198            198     -21.653168                   0
    199            199     -69.561989                   0
    
    [200 rows x 3 columns]
    
    Wrote phenotypes to gwas.example.grg.phen


Perform association between genotype and phenotype
--------------------------------------------------

Now we want to learn the association between the genotype :math:`X` and
the phenotype :math:`Y`. We are using simple linear regression to do
this. The result is a value :math:`\beta_i` for each mutation :math:`i`:
a larger absolute value indicates a larger effect of that mutation
(variant) on the phenotype. We are finding *correlations* between
:math:`X` and :math:`Y`. A variant :math:`i` that is highly correlated
with a phenotype :math:`Y` *might* have a causal relationship to
:math:`Y`.

GWAS via command line
~~~~~~~~~~~~~~~~~~~~~

First we’ll run the shell commands that perform a GWAS.

.. code:: bash

    %%bash
    
    # Using our phenotype file, emit a tab-separated (tsv) pandas dataframe containing the results of our GWAS between X (test-200-samples.grg) and Y (test-200-samples.grg.phen)
    grapp assoc -p gwas.example.grg.phen -o gwas.example.gwas.tsv gwas.example.grg


.. parsed-literal::

    Using the column 2 (the last one) for phenotype.


.. parsed-literal::

    Wrote results to gwas.example.gwas.tsv


Now we can examine the results by loading the dataframe into pandas.

.. code:: ipython3

    import pandas
    
    gwas_df = pandas.read_csv("gwas.example.gwas.tsv", delimiter="\t")
    gwas_df




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>POS</th>
          <th>ALT</th>
          <th>REF</th>
          <th>COUNT</th>
          <th>BETA</th>
          <th>B0</th>
          <th>SE</th>
          <th>R2</th>
          <th>T</th>
          <th>P</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>55829</td>
          <td>G</td>
          <td>A</td>
          <td>4</td>
          <td>-0.512129</td>
          <td>0.010243</td>
          <td>0.505040</td>
          <td>0.005166</td>
          <td>-1.014035</td>
          <td>0.311804</td>
        </tr>
        <tr>
          <th>1</th>
          <td>56812</td>
          <td>T</td>
          <td>G</td>
          <td>3</td>
          <td>0.794622</td>
          <td>-0.011919</td>
          <td>0.580456</td>
          <td>0.009376</td>
          <td>1.368961</td>
          <td>0.172562</td>
        </tr>
        <tr>
          <th>2</th>
          <td>57349</td>
          <td>G</td>
          <td>T</td>
          <td>1</td>
          <td>-1.466317</td>
          <td>0.007332</td>
          <td>0.999621</td>
          <td>0.010750</td>
          <td>-1.466873</td>
          <td>0.143997</td>
        </tr>
        <tr>
          <th>3</th>
          <td>58785</td>
          <td>T</td>
          <td>C</td>
          <td>10</td>
          <td>-0.358335</td>
          <td>0.017917</td>
          <td>0.324263</td>
          <td>0.006130</td>
          <td>-1.105077</td>
          <td>0.270467</td>
        </tr>
        <tr>
          <th>4</th>
          <td>59367</td>
          <td>A</td>
          <td>G</td>
          <td>2</td>
          <td>-0.118437</td>
          <td>0.001184</td>
          <td>0.712412</td>
          <td>0.000140</td>
          <td>-0.166248</td>
          <td>0.868131</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>10888</th>
          <td>9997601</td>
          <td>G</td>
          <td>C</td>
          <td>3</td>
          <td>0.113106</td>
          <td>-0.001697</td>
          <td>0.583141</td>
          <td>0.000190</td>
          <td>0.193959</td>
          <td>0.846407</td>
        </tr>
        <tr>
          <th>10889</th>
          <td>9998038</td>
          <td>A</td>
          <td>C</td>
          <td>21</td>
          <td>-0.072860</td>
          <td>0.007650</td>
          <td>0.209914</td>
          <td>0.000608</td>
          <td>-0.347095</td>
          <td>0.728889</td>
        </tr>
        <tr>
          <th>10890</th>
          <td>9998412</td>
          <td>G</td>
          <td>T</td>
          <td>42</td>
          <td>0.068393</td>
          <td>-0.014362</td>
          <td>0.164342</td>
          <td>0.000874</td>
          <td>0.416160</td>
          <td>0.677743</td>
        </tr>
        <tr>
          <th>10891</th>
          <td>9999031</td>
          <td>C</td>
          <td>G</td>
          <td>295</td>
          <td>-0.015402</td>
          <td>0.022718</td>
          <td>0.116634</td>
          <td>0.000088</td>
          <td>-0.132055</td>
          <td>0.895075</td>
        </tr>
        <tr>
          <th>10892</th>
          <td>9999126</td>
          <td>T</td>
          <td>A</td>
          <td>2</td>
          <td>1.195692</td>
          <td>-0.011957</td>
          <td>0.707376</td>
          <td>0.014225</td>
          <td>1.690320</td>
          <td>0.092540</td>
        </tr>
      </tbody>
    </table>
    <p>10893 rows × 10 columns</p>
    </div>



``BETA`` is the effect size for the variant at base-pair position
``POS`` with alternate allele ``ALT``. We can plot the histogram of our
inferred ``BETA`` values and see that it does indeed recover a normal
distribution centered at :math:`0`, which is what we simulated with our
phenotype.

.. code:: ipython3

    import seaborn
    seaborn.histplot(data=gwas_df, x="BETA")




.. parsed-literal::

    <Axes: xlabel='BETA', ylabel='Count'>




.. image:: GWAS_files/GWAS_12_1.png


GWAS via Python APIs
~~~~~~~~~~~~~~~~~~~~

Now we do the same GWAS, but using Python code instead of running
``grapp assoc ...`` on the command line. This gives us more options. In
this example, we’ll standardize the genotype matrix :math:`X` prior to
performing the linear regression, which will make all mutations
(variants) have the same variance. This has the effect of putting
variants with different allele frequencies on the same scale.

.. code:: ipython3

    from grapp.assoc import linear_assoc_no_covar
    import pygrgl
    
    GRG_FILE = "gwas.example.grg"
    PHEN_FILE = "gwas.example.grg.phen"
    
    # Load the GRG into memory
    grg = pygrgl.load_immutable_grg(GRG_FILE)
    
    # Load the phenotype into memory
    Y = pandas.read_csv(PHEN_FILE, delimiter="\t")
    
    # Perform the GWAS
    gwas_df = linear_assoc_no_covar(grg, Y["phenotypes"].to_numpy(), standardize=True)
    gwas_df




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>POS</th>
          <th>ALT</th>
          <th>REF</th>
          <th>COUNT</th>
          <th>BETA</th>
          <th>B0</th>
          <th>SE</th>
          <th>R2</th>
          <th>T</th>
          <th>P</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>55829</td>
          <td>G</td>
          <td>A</td>
          <td>4</td>
          <td>-0.071335</td>
          <td>0.001427</td>
          <td>0.070707</td>
          <td>0.005116</td>
          <td>-1.008875</td>
          <td>0.314266</td>
        </tr>
        <tr>
          <th>1</th>
          <td>56812</td>
          <td>T</td>
          <td>G</td>
          <td>3</td>
          <td>0.096223</td>
          <td>-0.001443</td>
          <td>0.070558</td>
          <td>0.009307</td>
          <td>1.363732</td>
          <td>0.174200</td>
        </tr>
        <tr>
          <th>2</th>
          <td>57349</td>
          <td>G</td>
          <td>T</td>
          <td>1</td>
          <td>-0.103295</td>
          <td>0.000516</td>
          <td>0.070508</td>
          <td>0.010724</td>
          <td>-1.465014</td>
          <td>0.144503</td>
        </tr>
        <tr>
          <th>3</th>
          <td>58785</td>
          <td>T</td>
          <td>C</td>
          <td>10</td>
          <td>-0.077090</td>
          <td>0.003854</td>
          <td>0.070676</td>
          <td>0.005988</td>
          <td>-1.090740</td>
          <td>0.276713</td>
        </tr>
        <tr>
          <th>4</th>
          <td>59367</td>
          <td>A</td>
          <td>G</td>
          <td>2</td>
          <td>-0.011755</td>
          <td>0.000118</td>
          <td>0.070884</td>
          <td>0.000139</td>
          <td>-0.165830</td>
          <td>0.868460</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>10888</th>
          <td>9997601</td>
          <td>G</td>
          <td>C</td>
          <td>3</td>
          <td>0.013696</td>
          <td>-0.000205</td>
          <td>0.070882</td>
          <td>0.000189</td>
          <td>0.193225</td>
          <td>0.846981</td>
        </tr>
        <tr>
          <th>10889</th>
          <td>9998038</td>
          <td>A</td>
          <td>C</td>
          <td>21</td>
          <td>-0.026328</td>
          <td>0.002764</td>
          <td>0.070864</td>
          <td>0.000704</td>
          <td>-0.371527</td>
          <td>0.710643</td>
        </tr>
        <tr>
          <th>10890</th>
          <td>9998412</td>
          <td>G</td>
          <td>T</td>
          <td>42</td>
          <td>0.029327</td>
          <td>-0.006159</td>
          <td>0.070857</td>
          <td>0.000903</td>
          <td>0.413891</td>
          <td>0.679402</td>
        </tr>
        <tr>
          <th>10891</th>
          <td>9999031</td>
          <td>C</td>
          <td>G</td>
          <td>295</td>
          <td>-0.009143</td>
          <td>0.013486</td>
          <td>0.070880</td>
          <td>0.000267</td>
          <td>-0.128993</td>
          <td>0.897494</td>
        </tr>
        <tr>
          <th>10892</th>
          <td>9999126</td>
          <td>T</td>
          <td>A</td>
          <td>2</td>
          <td>0.118671</td>
          <td>-0.001187</td>
          <td>0.070386</td>
          <td>0.014155</td>
          <td>1.686008</td>
          <td>0.093369</td>
        </tr>
      </tbody>
    </table>
    <p>10893 rows × 10 columns</p>
    </div>



We can see that the output format (dataframe) is the same as the command
line version. We should expect some differences in the distribution of
betas because we standardized :math:`X` this time, though the shape
should still be a normal distribution centered at :math:`0`.

.. code:: ipython3

    seaborn.histplot(data=gwas_df, x="BETA")




.. parsed-literal::

    <Axes: xlabel='BETA', ylabel='Count'>




.. image:: GWAS_files/GWAS_16_1.png


You can see that the range of beta values is much smaller (from
:math:`-0.3` to :math:`+0.3`).

Finally, we can perform one more GWAS where we use an approximate
distribution for the variance of each mutation in :math:`X`. Each value
in the :math:`X` matrix can take on a value of :math:`0`, :math:`1`, or
:math:`2`, based on how many copies each of the (diploid) individuals
has for that allele. By default, the linear regression is performed with
the exact sum of squared errors, based on these diploid values in
:math:`X`. Theoretically, in the absence of many complicating effects
(selection, genotyping error, inbreeding, etc.) the variance of each
mutation should follow a binomial distribution with mean
:math:`2 \times f_i` (where :math:`f_i` is the allele frequency for
mutation :math:`i`). You can choose to perform the regression using the
binomial distribution, instead of the exact values, via
``dist="binomial"``, as shown below.

.. code:: ipython3

    gwas_df = linear_assoc_no_covar(grg, Y["phenotypes"].to_numpy(), dist="binomial")
    seaborn.histplot(data=gwas_df, x="BETA")




.. parsed-literal::

    <Axes: xlabel='BETA', ylabel='Count'>




.. image:: GWAS_files/GWAS_18_1.png


You can see that this (non-standardized) distribution follows our
original (exact) result pretty closely, even with a small sample size of
200 individuals. Partially this is because our genotype matrix :math:`X`
was generated using a neutral simulation; in real datasets, the binomial
distribution may cause more beta values to differ from the exact value
than is shown here.

Related Topics
--------------

-  In practice, GWAS is often performed while handling covariates (such
   as age, sex, principal components, etc). See `GWAS with
   Covariates <GWASCovariates.html>`__.
-  See `Simulating Phenotypes <SimulatingPhenotypes.html>`__ for a more
   details on generating synthetic phenotypes.
-  Documentation links:

   -  `grapp.assoc <https://grapp.readthedocs.io/en/latest/grapp.html#module-grapp.assoc>`__:
      Python APIs for GWAS
