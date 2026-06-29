import io
import numpy
import os
import pygrgl
import tskit
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
from testing_utils import allele_frequencies

# Tree defined looks like this (see tskit.TreeSequence.draw_text()):
#
# 100.00|       14        |        14       |
#       |    +---+---+    |     +---+---+   |
# 95.00 |   12       |    |    12       |   |
#       |  +-+-+     |    |   +-+--+    |   |
# 70.00 |  8   |     |    |   8    |    |   |
#       | +++  |     |    | +-+-+  |    |   |
# 61.00 | | |  9     |    | | | |  9    |   |
#       | | | +++    |    | | | | +++   |   |
# 60.00 | | | | |   13    | | | | | |  13   |
#       | | | | |  +-+-+  | | | | | |  ++-+ |
# 50.00 | | | | | 10   |  | | | | | | 10  | |
#       | | | | | +++  |  | | | | | | +++ | |
# 43.00 | | | | | | | 11  | | | | | | | |11 |
#       | | | | | | | +++ | | | | | | | | | |
# 0.00  | 0 1 2 3 4 5 6 7 | 0 1 7 2 3 4 5 6 |
#       0                 5                10

TEST_NODES = """\
is_sample   individual   time
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
0           -1           70.0
0           -1           61.0
0           -1           50.0
0           -1           43.0
0           -1           95.0
0           -1           60.0
0           -1           100.0
"""

TEST_EDGES = """\
left   right   parent  child
0.0    10.0    8       0
0.0    10.0    8       1
0.0    10.0    9       2
0.0    10.0    9       3
0.0    10.0    10      4
0.0    10.0    10      5
0.0    10.0    11      6
0.0     5.0    11      7
5.0    10.0    8       7
0.0    10.0    12      8
0.0    10.0    12      9
0.0    10.0    13      10
0.0    10.0    13      11
0.0    10.0    14      13
0.0    10.0    14      12
"""

TEST_SITES = """\
position      ancestral_state
1.0           A
2.0           A
4.0           A
8.0           T
"""

TEST_MUTS = """\
site   node    derived_state    time    parent
0      13      T                62.0     -1
1      12      G                97.0     -1
2      13      C                73.0     -1
3      8       A                77.0     -1
"""

TEST_TS = tskit.load_text(
    nodes=io.StringIO(TEST_NODES),
    edges=io.StringIO(TEST_EDGES),
    sites=io.StringIO(TEST_SITES),
    mutations=io.StringIO(TEST_MUTS),
    strict=False,
)


# Back/multiply-mapped mutation test:
#
# 100.00┊      14         ┊
#       ┊    ┏━━┻━━━┓     ┊
# 95.00 ┊   12      ┃     ┊
#       ┊  ┏━┻━┓    ┃     ┊
# 70.00 ┊  8   ┃    ┃     ┊
#       ┊ ┏┻┓  ┃    ┃     ┊
# 61.00 ┊ ┃ ┃  9    ┃     ┊
#       ┊ ┃ ┃ ┏┻┓   ┃     ┊
# 60.00 ┊ ┃ ┃ ┃ ┃  13     ┊
#       ┊ ┃ ┃ ┃ ┃  ┏┻━┓   ┊
# 50.00 ┊ ┃ ┃ ┃ ┃ 10  ┃   ┊
#       ┊ ┃ ┃ ┃ ┃ ┏┻┓ ┃   ┊
# 43.00 ┊ ┃ ┃ ┃ ┃ ┃ ┃11   ┊
#       ┊ ┃ ┃ ┃ ┃ ┃ ┃ ┃   ┊
# 0.00  ┊ 0 1 2 3 4 5 6 7 ┊
#       0                10

NESTED_NODES = """\
is_sample   individual   time
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
1           -1           0.0
0           -1           70.0
0           -1           61.0
0           -1           50.0
0           -1           43.0
0           -1           95.0
0           -1           60.0
0           -1           100.0
"""

NESTED_EDGES = """\
left   right   parent  child
0.0    10.0    8       0
0.0    10.0    8       1
0.0    10.0    9       2
0.0    10.0    9       3
0.0    10.0    10      4
0.0    10.0    10      5
0.0    10.0    11      6
0.0    10.0    12      8
0.0    10.0    12      9
0.0    10.0    13      10
0.0    10.0    13      11
0.0    10.0    14      13
0.0    10.0    14      12
"""

NESTED_SITES = """\
position      ancestral_state
1.0           A
8.0           T
"""

# Notice the two back muts with parents
NESTED_MUTS = """\
site   node    derived_state    time    parent
0      12      T                96.0     -1
0      8       A                80.0     0
0      9       C                70.0     0
1      8       G                81.0     -1
1      10      G                55.0     -1
1      9       G                69.0     -1
"""

NESTED_TS = tskit.load_text(
    nodes=io.StringIO(NESTED_NODES),
    edges=io.StringIO(NESTED_EDGES),
    sites=io.StringIO(NESTED_SITES),
    mutations=io.StringIO(NESTED_MUTS),
    strict=False,
)


class TestTS2GRG(unittest.TestCase):
    def test_simple_ts(self):
        ts_file = "test.simple_ts.trees"
        TEST_TS.dump(ts_file)
        grg = pygrgl.grg_from_trees(ts_file)
        self.assertEqual(grg.num_samples, 8)
        self.assertEqual(grg.num_mutations, 4)
        self.assertEqual(grg.ploidy, 2)
        self.assertTrue(grg.is_phased)
        af = allele_frequencies(grg)
        numpy.testing.assert_allclose(af, [0.5, 0.5, 0.5, 0.375])

    def test_back_mut(self):
        ts_file = "test.back_mut.trees"
        NESTED_TS.dump(ts_file)
        self.assertEqual(NESTED_TS.num_mutations, 6)

        # This tree sequence is really tricky. The back mutations look like this:
        #
        #             | A (ancestral)
        #             *
        #             | T (parent)
        #             *
        # A (child) /   \  C (child)
        #          *    *
        # [2 samples]  [2 samples]
        #
        # The A->T->A means that the left tree cancels out and has no mutations, from the GRG p.o.v.,
        # while the right tree has 2 samples with mutation (A->C). So while the original TS has 5 mutations,
        # 2 of them become no-ops (A->T and T->A) and should not show up in the GRG.
        #
        # The other site (site 1) has two recurrent mutations (not nested), which should be combined into
        # a single GRG mutation.

        grg = pygrgl.grg_from_trees(ts_file)
        grg.sort_mutations()
        self.assertTrue(grg.mutations_are_unique())
        self.assertEqual(grg.num_samples, 8)
        self.assertEqual(grg.num_mutations, 3)
        mut1 = grg.get_mutation_by_id(0)
        self.assertEqual((mut1.position, mut1.allele), (1.0, "C"))
        mut2 = grg.get_mutation_by_id(1)
        self.assertEqual((mut2.position, mut2.allele), (1.0, "T"))
        mut3 = grg.get_mutation_by_id(2)
        self.assertEqual((mut3.position, mut3.allele), (8.0, "G"))
        self.assertEqual(grg.ploidy, 2)
        self.assertTrue(grg.is_phased)
        af = allele_frequencies(grg)
        numpy.testing.assert_allclose(af, [0.25, 0.0, 0.75])
