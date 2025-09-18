import glob
import numpy as np
import os
import pygrgl
import subprocess
import sys
import unittest
from typing import List, Dict, Tuple, Optional

JOBS = 4
CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(THIS_DIR)
from testing_utils import construct_grg

EXPECT_DIR = os.path.join(THIS_DIR, "expect")
INPUT_DIR = os.path.join(THIS_DIR, "input")


class TestGrgMissingData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grg_filename = construct_grg("test-200-samples.miss.igd")
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename)

    # Test that the mutations have corresponding missiness nodes
    def test_muts(self):
        pass

    # Test that the allele count computation handles missingness properly
    def test_allele_counts(self):
        pass

    # Test that ignoring the missing data results in them all becoming REF
    def test_drop_missing(self):
        pass

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)


if __name__ == "__main__":
    unittest.main()
