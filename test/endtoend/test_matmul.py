import unittest
import subprocess
import os
import numpy as np
import pygrgl
import glob
from typing import List, Dict, Tuple, Optional

JOBS = 4
CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

INPUT_DIR = os.path.join(THIS_DIR, "input")


# Absurdly innefficient
def grg2X(grg: pygrgl.GRG, individual: bool = False):
    num_samples = grg.num_samples if not individual else grg.num_individuals
    result = np.zeros((num_samples, grg.num_mutations))
    muts_above = {}
    for node_id in reversed(range(grg.num_nodes)):
        muts = grg.get_mutations_for_node(node_id)
        ma = []
        if muts:
            ma.extend(muts)
        for parent_id in grg.get_up_edges(node_id):
            ma.extend(muts_above[parent_id])
        muts_above[node_id] = ma
        if grg.is_sample(node_id):
            sample_index = node_id if not individual else (node_id // grg.ploidy)
            for mut_id in muts_above[node_id]:
                result[sample_index][mut_id] += 1
    return result


def construct_grg(input_file: str, output_file: Optional[str] = None) -> str:
    cmd = [
        "grg",
        "construct",
        "-p",
        "10",
        "-t",
        "2",
        "-j",
        str(JOBS),
        os.path.join(INPUT_DIR, input_file),
    ]
    if output_file is not None:
        cmd.extend(["-o", output_file])
    else:
        output_file = os.path.basename(input_file) + ".final.grg"
    subprocess.check_call(cmd)
    return output_file


class TestMatrixMultiplication(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grg_filename = construct_grg("test-200-samples.vcf.gz", "test.matmul.grg")
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename)
        np.random.seed(42)

    def direction_helper(
        self, grg: pygrgl.GRG, size: int, direction: pygrgl.TraversalDirection
    ):
        K = 100

        # Construct input vectors.
        list_of_vects = []
        for _ in range(K):
            list_of_vects.append(np.random.rand(size))

        # Result 1: K independent dot_product calls to GRG
        dot_prod_result = []
        for i in range(K):
            dot_prod_result.append(pygrgl.dot_product(grg, list_of_vects[i], direction))
        dot_prod_result = np.array(dot_prod_result)

        # Result 2: Numpy matrix multiplication on the genotype matrix
        genotype_matrix = grg2X(grg)
        if direction == pygrgl.TraversalDirection.UP:
            np_result = np.matmul(np.array(list_of_vects), genotype_matrix)
        else:
            np_result = np.matmul(
                np.array(list_of_vects), np.transpose(genotype_matrix)
            )

        # Result 3: GRG matrix multiplication
        matmul_result = pygrgl.matmul(grg, np.array(list_of_vects), direction)

        self.assertTrue(np.allclose(dot_prod_result, np_result))
        self.assertTrue(np.allclose(dot_prod_result, matmul_result))

    def test_different_methods(self):
        self.direction_helper(
            self.grg, self.grg.num_samples, pygrgl.TraversalDirection.UP
        )
        self.direction_helper(
            self.grg, self.grg.num_mutations, pygrgl.TraversalDirection.DOWN
        )

    def test_diploid(self):
        # Test that diploid operations are the same between GRG and numpy
        K = 10

        # Test the (KxN)(NxM) multiplication.
        genotype_matrix = grg2X(self.grg, individual=True)
        rand_matrix = np.random.rand(K, self.grg.num_individuals)
        np_result = np.matmul(rand_matrix, genotype_matrix)
        grg_result = pygrgl.matmul(
            self.grg, rand_matrix, pygrgl.TraversalDirection.UP, by_individual=True
        )
        self.assertTrue(np.allclose(grg_result, np_result))

        # Test the (KxM)(MxN) multiplication.
        genotype_matrix = np.transpose(grg2X(self.grg, individual=True))
        rand_matrix = np.random.rand(K, self.grg.num_mutations)
        np_result = np.matmul(rand_matrix, genotype_matrix)
        grg_result = pygrgl.matmul(
            self.grg, rand_matrix, pygrgl.TraversalDirection.DOWN, by_individual=True
        )
        self.assertTrue(np.allclose(grg_result, np_result))

    def test_xtx_init(self):
        # Test that (X^t * X) can be computed by setting init properly.
        node_XX_count = [0 for _ in range(self.grg.num_nodes)]
        assert self.grg.ploidy == 2
        for node_id in range(self.grg.num_nodes):
            curr_coals = self.grg.get_num_individual_coals(node_id)
            assert curr_coals != pygrgl.COAL_COUNT_NOT_SET
            assert curr_coals <= self.grg.num_samples / 2
            coal_modifier = 2 * curr_coals
            if self.grg.is_sample(node_id):
                node_XX_count[node_id] = 1
            else:
                count = sum(
                    [
                        node_XX_count[child_id]
                        for child_id in self.grg.get_down_edges(node_id)
                    ]
                )
                node_XX_count[node_id] = count + coal_modifier

        matmul_count = pygrgl.matmul(
            self.grg,
            np.ones((2, self.grg.num_samples), dtype=np.int32),
            pygrgl.TraversalDirection.UP,
            init="xtx",
        )
        for mut_id, node_id in self.grg.get_mutation_node_pairs():
            if node_id == pygrgl.INVALID_NODE:
                continue
            self.assertEqual(node_XX_count[node_id], matmul_count[0][mut_id])
            self.assertEqual(node_XX_count[node_id], matmul_count[1][mut_id])

    def test_vector_init(self):
        dtype = np.float64
        K = 10
        init = np.ones(K, dtype=dtype) * 2

        without_init = pygrgl.matmul(
            self.grg,
            np.ones((K, self.grg.num_samples), dtype=dtype),
            pygrgl.TraversalDirection.UP,
        )
        with_init = pygrgl.matmul(
            self.grg,
            np.ones((K, self.grg.num_samples), dtype=dtype),
            pygrgl.TraversalDirection.UP,
            init=init,
        )

        # The sample -> node relationship is one-to-one, so we know this should hold
        self.assertEqual(without_init.shape, with_init.shape)
        self.assertTrue(np.all(with_init >= 2 * without_init))

    def test_matrix_init(self):
        dtype = np.int64
        K = 10
        init = np.zeros((K, self.grg.num_nodes), dtype=dtype)

        without_init = pygrgl.matmul(
            self.grg,
            np.ones((K, self.grg.num_samples), dtype=dtype),
            pygrgl.TraversalDirection.UP,
        )
        with_init = pygrgl.matmul(
            self.grg,
            np.ones((K, self.grg.num_samples), dtype=dtype),
            pygrgl.TraversalDirection.UP,
            init=init,
        )

        # The sample -> node relationship is one-to-one, so we know this should hold
        self.assertEqual(without_init.shape, with_init.shape)
        self.assertTrue(np.array_equal(without_init, with_init))

        # Pick a random node, get all muts above.
        tweak_id = np.random.randint(self.grg.num_samples, self.grg.num_nodes - 100)
        collected_mut_ids = []

        def collect(start_id):
            muts = self.grg.get_mutations_for_node(start_id)
            collected_mut_ids.extend(muts)
            for parent_id in self.grg.get_up_edges(start_id):
                collect(parent_id)

        collect(tweak_id)
        # Change the init value _only_ for that node, so everything above it will be impacted.
        init[4, tweak_id] = 100
        # print(f"Muts at tweak_id = {self.grg.get_mutations_for_node(tweak_id)}")
        # print(f"Above: {collected_mut_ids}")

        with_init = pygrgl.matmul(
            self.grg,
            np.ones((K, self.grg.num_samples), dtype=dtype),
            pygrgl.TraversalDirection.UP,
            init=init,
        )

        # The sample -> node relationship is one-to-one, so we know this should hold
        self.assertEqual(without_init.shape, with_init.shape)
        for i in range(with_init.shape[0]):
            for j in range(with_init.shape[1]):
                if i != 4:
                    # All other rows should be unchanged.
                    self.assertEqual(with_init[i, j], without_init[i, j])
                else:
                    # Only mutations above our tweaked node should be impacted.
                    if j not in collected_mut_ids:
                        self.assertEqual(with_init[i, j], without_init[i, j])
                    else:
                        self.assertGreater(with_init[i, j], without_init[i, j])

    def test_init_failures(self):
        # Mismatching dimensions on init
        with self.assertRaises(RuntimeError):
            pygrgl.matmul(
                self.grg,
                np.ones((1, self.grg.num_samples)),
                pygrgl.TraversalDirection.UP,
                init=np.ones((1, 200000)),
            )
        # Mismatching dimensions on init
        with self.assertRaises(RuntimeError):
            pygrgl.matmul(
                self.grg,
                np.ones((1, self.grg.num_samples)),
                pygrgl.TraversalDirection.UP,
                init=np.ones((2, self.grg.num_samples)),
            )
        # Mismatching dtype on init
        with self.assertRaises(RuntimeError):
            pygrgl.matmul(
                self.grg,
                np.ones((1, self.grg.num_samples), dtype=np.int32),
                pygrgl.TraversalDirection.UP,
                init=np.ones((1, self.grg.num_samples)),
            )
        # Invalid string
        with self.assertRaises(RuntimeError):
            pygrgl.matmul(
                self.grg,
                np.ones((1, self.grg.num_samples), dtype=np.int32),
                pygrgl.TraversalDirection.UP,
                init="test",
            )
        # Vector too long
        with self.assertRaises(RuntimeError):
            pygrgl.matmul(
                self.grg,
                np.ones((9, self.grg.num_samples)),
                pygrgl.TraversalDirection.UP,
                init=np.ones(10),
            )

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)


if __name__ == "__main__":
    unittest.main()
