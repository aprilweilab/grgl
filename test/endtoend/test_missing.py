import numpy
import os
import pygrgl
import sys
import unittest

JOBS = 4
CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(THIS_DIR)
from testing_utils import construct_grg, grg2X, allele_frequencies

EXPECT_DIR = os.path.join(THIS_DIR, "expect")
INPUT_DIR = os.path.join(THIS_DIR, "input")


class TestGrgMissingData(unittest.TestCase):
    # Properties of the input data.
    MISSING_INDIVS = 21
    MISSING_SAMPLES = 25

    @classmethod
    def setUpClass(cls):
        cls.grg_filename = construct_grg("test-200-samples.miss.igd")
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename, load_up_edges=False)

    # Expensive, so only use sparingly.
    def get_samples_below(self, node):
        result = []
        for c in self.grg.get_down_edges(node):
            if self.grg.is_sample(c):
                result.append(c)
            else:
                result.extend(self.get_samples_below(c))
        return result

    # Test that the mutations have corresponding missiness nodes
    def test_muts(self):
        total_miss_idv = 0
        total_miss = 0
        for mut, node, miss in self.grg.get_mutation_node_miss():
            mut = self.grg.get_mutation_by_id(mut)
            if node != pygrgl.INVALID_NODE and miss != pygrgl.INVALID_NODE:
                rs = self.get_samples_below(node)
                assert len(set(rs)) == len(rs)
                ms = self.get_samples_below(miss)
                assert len(set(ms)) == len(ms)
                total_miss += len(ms)
                total_miss_idv += len(set([i // 2 for i in ms]))
        self.assertEqual(total_miss, self.MISSING_SAMPLES)
        self.assertEqual(total_miss_idv, self.MISSING_INDIVS)

    # Test matmul() support for missing data.
    def test_mult(self):
        # X is the explicit genotype matrix, with allele frequency used for missing items. So the
        # only non-0,1,2 values should be missing items.
        X = grg2X(self.grg, diploid=True)
        self.assertEqual(
            len(numpy.where((X > 0) & (X != 1) & (X != 2))[0]), self.MISSING_INDIVS
        )

        #### UP direction (AX) ####
        K = 7
        rv = numpy.random.standard_normal((K, self.grg.num_individuals))
        # Multiply genotype matrix with a random matrix.
        numpy_result = rv @ X
        # Multiply GRG with the same random matrix.
        freqs = allele_frequencies(self.grg)
        M = numpy.zeros((K, self.grg.num_mutations))
        grg_result = pygrgl.matmul(
            self.grg, rv, pygrgl.TraversalDirection.UP, by_individual=True, miss=M
        )
        # We can adjust the GRG result with the M matrix to get an equivalent result to having
        # the explicit matrix with missing values set to the mean.
        numpy.testing.assert_allclose(numpy_result, grg_result + (M * freqs))

        #### DOWN direction (AX^T) ####
        rv = numpy.random.standard_normal((K, self.grg.num_mutations))
        # Multiply genotype matrix with a random matrix.
        numpy_result = rv @ X.T
        # Multiply GRG with the same random matrix.
        M = numpy.array([freqs]) * rv
        grg_result = pygrgl.matmul(
            self.grg, rv, pygrgl.TraversalDirection.DOWN, by_individual=True, miss=M
        )
        # We primed the missingness matrix with M=f_i*input, to weight those results, so our GRG result
        # should be identical to the explicit matrix result.
        numpy.testing.assert_allclose(numpy_result, grg_result)

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)


if __name__ == "__main__":
    unittest.main()
