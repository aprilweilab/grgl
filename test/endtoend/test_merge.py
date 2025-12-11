import glob
import numpy as np
import os
import pygrgl
import subprocess
import sys
import unittest
from typing import List, Dict, Tuple

JOBS = 4
CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(THIS_DIR)
from testing_utils import construct_grg


class TestGrgMerge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.left_filename = construct_grg("merge_left.vcf", parts=1)
        cls.right_filename = construct_grg("merge_right.vcf", parts=1)
        cls.left_grg = pygrgl.load_mutable_grg(cls.left_filename)
        cls.right_grg = pygrgl.load_mutable_grg(cls.right_filename)

    def test_no_error(self):
        """
        The two VCFs we are merging overlap, in that the start/end position is the same. This
        is the only allowed overlap (if the positions overlap _past_ each other, that is bad)
        """
        self.left_grg.merge([self.right_filename])

    @classmethod
    def tearDownClass(cls):
        cls.left_grg = None
        cls.right_grg = None
        if CLEANUP:
            os.remove(cls.left_filename)
            os.remove(cls.right_filename)


if __name__ == "__main__":
    unittest.main()
