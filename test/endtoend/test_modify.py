import numpy
import os
import pygrgl
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
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename, load_up_edges=True)

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)

    def test_downsample(self):
        sub_grg_filename = "test.downsample.grg"
        saved = pygrgl.save_subset(
            self.grg,
            sub_grg_filename,
            pygrgl.TraversalDirection.UP,
            list(range(0, 200, 10)),
        )
        self.assertTrue(saved)
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

    def test_downsample_fail(self):
        sub_grg_filename = "test.downsample.fail.grg"
        saved = pygrgl.save_subset(
            self.grg,
            sub_grg_filename,
            pygrgl.TraversalDirection.DOWN,
            [],
        )
        self.assertFalse(saved)
        self.assertFalse(os.path.isfile(sub_grg_filename))

    def test_reduce(self):
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=True)
        iterations = grg.reduce_until(iterations=5, min_dropped=10, verbose=True)
        self.assertLess(grg.num_edges, self.grg.num_edges)
        self.assertGreater(grg.num_nodes, self.grg.num_nodes)
        self.assertLessEqual(iterations, 5)
        self.assertGreaterEqual(iterations, 1)

    def test_mutations(self):
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=False)
        orig_count = grg.num_mutations
        # Pick a (deterministic) random mutation
        mut_id = int(grg.num_mutations / 2.5)
        for m, node_id in grg.get_mutation_node_pairs():
            if m == mut_id:
                break
        old_mut = grg.get_mutation_by_id(mut_id)
        self.assertNotEqual(node_id, pygrgl.INVALID_NODE)
        self.assertTrue(grg.mutations_are_ordered)
        grg.remove_mutation(mut_id, node_id)
        self.assertFalse(grg.mutations_are_ordered)
        # Pick a different random node and add two mutations to it.
        self.assertTrue(grg.num_nodes > 1001)
        new_mut_id1 = grg.add_mutation(pygrgl.Mutation(0, "FAKE ALT", "00"), 500)
        new_mut_id2 = grg.add_mutation(pygrgl.Mutation(500_000_000, "FAKE ALT2", "002"), 1001)
        self.assertNotEqual(new_mut_id1, new_mut_id2)
        self.assertFalse(grg.mutations_are_ordered)
        # The node->mut map should be updated properly as well
        self.assertEqual(grg.get_mutations_for_node(500)[-1], new_mut_id1)
        self.assertEqual(grg.get_mutations_for_node(1001)[-1], new_mut_id2)

        # Now sort them and verify the new ones are in the right order
        grg.sort_mutations()
        self.assertTrue(grg.mutations_are_ordered)
        self.assertEqual(grg.num_mutations, orig_count + 1)
        new_mut1 = grg.get_mutation_by_id(0)
        self.assertEqual(new_mut1.allele, "FAKE ALT")
        self.assertEqual(new_mut1.position, 0)
        new_mut2 = grg.get_mutation_by_id(grg.num_mutations - 1)
        self.assertEqual(new_mut2.allele, "FAKE ALT2")
        self.assertEqual(new_mut2.position, 500_000_000)
        # Make sure the old one is gone
        for i in range(grg.num_mutations):
            mut = grg.get_mutation_by_id(i)
            self.assertNotEqual(
                (old_mut.position, old_mut.allele),
                (mut.position, mut.allele),
            )


if __name__ == "__main__":
    unittest.main()
