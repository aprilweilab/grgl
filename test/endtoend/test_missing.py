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
        self.assertEqual(total_miss, 25)
        self.assertEqual(total_miss_idv, 21)

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)


if __name__ == "__main__":
    unittest.main()
