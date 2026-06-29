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
        new_mut_id2 = grg.add_mutation(
            pygrgl.Mutation(500_000_000, "FAKE ALT2", "002"), 1001
        )
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

    def test_set_samples(self):
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=True)
        self.assertEqual(grg.num_samples, 400)
        self.assertTrue(grg.samples_are_ordered)
        # Just remove a couple samples, don't reorder them.
        grg.set_samples(list(range(50, grg.num_samples)))
        self.assertEqual(grg.num_samples, 350)
        self.assertFalse(grg.samples_are_ordered)

        # Do matmul and get the value for sample 10 (which is "old" sample 60)
        inmat = numpy.ones((1, grg.num_mutations))
        result = pygrgl.matmul(grg, inmat, pygrgl.TraversalDirection.DOWN)[0]
        old_value = result[10]

        # Now reorder the samples by swapping sample 10 and 11 (60 and 61)
        new_samples = grg.get_sample_nodes()
        tmp = new_samples[10]
        new_samples[10] = new_samples[11]
        new_samples[11] = tmp
        grg.set_samples(new_samples)
        self.assertEqual(grg.num_samples, 350)
        self.assertFalse(grg.samples_are_ordered)

        # Check that the value followed the node, not the sample "index"!
        result = pygrgl.matmul(grg, inmat, pygrgl.TraversalDirection.DOWN)[0]
        wrong_value = result[10]
        new_value = result[11]
        self.assertEqual(old_value, new_value)
        self.assertNotEqual(old_value, wrong_value)

    def test_negative_nodes(self):
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=True)
        self.assertTrue(grg.nodes_are_ordered)

        regular_new = grg.make_node()
        self.assertTrue(grg.nodes_are_ordered)
        grg.connect(regular_new, grg.num_nodes - 2)  # Connect as parent: ok
        self.assertTrue(grg.nodes_are_ordered)

        negative_new = grg.make_node(negative=True)
        self.assertFalse(grg.nodes_are_ordered)
        grg.connect(
            0, -negative_new
        )  # Connect negative as child: breaks counting order
        self.assertFalse(grg.nodes_are_ordered)

        topo_list = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.UP, [])
        self.assertEqual(len(topo_list), grg.num_nodes)
        self.assertEqual(topo_list[0], negative_new)

        # Now TRULY break the topological order, which will force the get_topo_order() to run a DFS
        # and produce a new order
        negative_new = grg.make_node(negative=True)
        grg.connect(-negative_new, 10)  # Break the order entirely

        topo_list = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.UP, [])
        # We have some arbitrary order that does not match the node ID order
        self.assertNotEqual(list(sorted(topo_list)), topo_list)

    def test_regr_negnodes_matmul(self):
        """
        Regression test for the following scenario: load a MutableGRG, add some positive and
        negative nodes, interleaved (alternating between them). Verify that the topological
        order produces the correct order (negative nodes first), other newly added nodes last.
        Verify that DOWN matmul() produces the expected result. Verify that matmul() after
        set_samples() produces the expected result.
        """
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=True)

        # Check our order invariant before any changes.
        topo_up = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.UP, [])
        topo_down = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.DOWN, [])
        self.assertEqual(topo_up, list(reversed(topo_down)))

        # Get our original matmul result
        rand_input = numpy.random.standard_normal((1, grg.num_mutations))
        orig_product = pygrgl.matmul(grg, rand_input, pygrgl.TraversalDirection.DOWN)[0]

        # Just add two rows of new nodes above the previous roots and below the previous samples.
        orig_node_ct = grg.num_nodes
        orig_roots = grg.get_root_nodes()
        orig_samples = grg.get_sample_nodes()

        above1 = []
        for r in orig_roots:
            a = grg.make_node()
            grg.connect(a, r)
            above1.append(a)

        below1 = []
        for s in orig_samples:
            b = grg.make_node(negative=True)
            grg.connect(s, -b)
            below1.append(b)

        above2 = []
        for r in above1:
            a = grg.make_node()
            grg.connect(a, r)
            above2.append(a)

        below2 = []
        for s in below1:
            b = grg.make_node(negative=True)
            grg.connect(s, -b)
            below2.append(b)

        self.assertFalse(grg.nodes_are_ordered)

        new_below = len(orig_samples) * 2
        new_above = len(orig_roots) * 2

        # Compare to our original matmul result. We have not modified the Mutations or Samples,
        # so this should be identical!
        product2 = pygrgl.matmul(grg, rand_input, pygrgl.TraversalDirection.DOWN)[0]
        numpy.testing.assert_allclose(orig_product, product2)

        topo_up = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.UP, [])

        # Correct counts.
        self.assertEqual(len(topo_up), grg.num_nodes)
        self.assertEqual(len(topo_up), orig_node_ct + new_above + new_below)
        # Every node appears only once.
        self.assertEqual(len(topo_up), len(set(topo_up)))
        # The new "below" nodes should be at the beginning
        self.assertEqual(
            topo_up[:new_below], list(reversed(below2)) + list(reversed(below1))
        )
        # The new "above" nodes should be at the end
        self.assertEqual(topo_up[-new_above:], above1 + above2)

        topo_down = pygrgl.get_topo_order(grg, pygrgl.TraversalDirection.DOWN, [])
        # Correct counts.
        self.assertEqual(len(topo_down), grg.num_nodes)
        self.assertEqual(len(topo_down), orig_node_ct + new_above + new_below)
        # Every node appears only once.
        self.assertEqual(len(topo_down), len(set(topo_down)))
        # The new "above" nodes should be at the beginning
        self.assertEqual(
            topo_down[:new_above], list(reversed(above2)) + list(reversed(above1))
        )
        # The new "below" nodes should be at the end
        self.assertEqual(topo_down[-new_below:], below1 + below2)

        # The two directions should be equivalent, just opposite orders.
        self.assertEqual(topo_up, list(reversed(topo_down)))

        # Save our new samples, and verify that this still doesn't change the matmul result.
        grg.set_samples(below2)
        self.assertEqual(grg.get_sample_nodes(), below2)
        product3 = pygrgl.matmul(grg, rand_input, pygrgl.TraversalDirection.DOWN)[0]
        numpy.testing.assert_allclose(orig_product, product3)

        # save_subset() should reject us, but save_grg() should work fine.
        with self.assertRaises(RuntimeError):
            pygrgl.save_subset(
                grg, "not_created.grg", pygrgl.TraversalDirection.UP, below2[:10]
            )
        pygrgl.save_grg(grg, "test.throwaway.grg")

    @unittest.skip(
        "Requires pip installing pygrgl with GRGL_CHECK_NEGATIVE=1 in the environment"
    )
    def test_check_negative(self):
        grg = pygrgl.load_mutable_grg(self.grg_filename, load_up_edges=True)

        # Scenario 1: negative node, passed to connect() without negative sign
        neg_i = grg.make_node(negative=True)
        with self.assertRaises(RuntimeError):
            grg.connect(0, neg_i)
        with self.assertRaises(RuntimeError):
            grg.connect(neg_i, 10)

        # Scenario 2: positive node, passed to connect() with negative sign
        pos_j = grg.make_node()
        with self.assertRaises(RuntimeError):
            grg.connect(0, -pos_j)
        with self.assertRaises(RuntimeError):
            grg.connect(-pos_j, 100)

        # Correct signs for both cases.
        grg.connect(0, -neg_i)
        grg.connect(0, pos_j)


if __name__ == "__main__":
    unittest.main()
