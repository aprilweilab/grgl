import glob
import numpy
import os
import pygrgl
import subprocess
import sys
import unittest

CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(THIS_DIR)
from testing_utils import construct_grg

EXPECT_DIR = os.path.join(THIS_DIR, "expect")
INPUT_DIR = os.path.join(THIS_DIR, "input")


class TestGrgModify(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grg_filename = construct_grg("test-200-samples.vcf.gz")
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename)

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)

    def test_downsample(self):
        sub_grg_filename = "test.downsample.grg"
        pygrgl.save_subset(
            self.grg,
            sub_grg_filename,
            pygrgl.TraversalDirection.UP,
            list(range(0, 200, 10)),
        )
        sub_grg = pygrgl.load_immutable_grg(sub_grg_filename)
        self.assertEqual(sub_grg.num_mutations, self.grg.num_mutations)
        self.assertEqual(sub_grg.num_samples, 20)

        input_full = numpy.zeros((1, self.grg.num_samples), dtype=numpy.int32)
        for i in range(0, 200, 10):
            input_full[0][i] = 1
        acount_full = pygrgl.matmul(self.grg, input_full, pygrgl.TraversalDirection.UP)[
            0
        ]

        input_sub = numpy.ones((1, sub_grg.num_samples), dtype=numpy.int32)
        acount_sub = pygrgl.matmul(sub_grg, input_sub, pygrgl.TraversalDirection.UP)[0]

        self.assertEqual(acount_full.shape, acount_sub.shape)
        numpy.testing.assert_array_equal(acount_full, acount_sub)


if __name__ == "__main__":
    unittest.main()
