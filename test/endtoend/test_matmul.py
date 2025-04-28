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
        # Create the GRG
        grg_filename = construct_grg("test-200-samples.vcf.gz", "test.matmul.grg")

        grg = pygrgl.load_immutable_grg(grg_filename)

        self.direction_helper(grg, grg.num_samples, pygrgl.TraversalDirection.UP)
        self.direction_helper(grg, grg.num_mutations, pygrgl.TraversalDirection.DOWN)

        if CLEANUP:
            os.remove(grg_filename)

    def test_diploid(self):
        # Test that diploid operations are the same between GRG and numpy
        grg_filename = construct_grg(
            "test-200-samples.vcf.gz", "test.matmul_diploid.grg"
        )

        grg = pygrgl.load_immutable_grg(grg_filename)

        K = 10

        # Test the (KxN)(NxM) multiplication.
        genotype_matrix = grg2X(grg, individual=True)
        rand_matrix = np.random.rand(K, grg.num_individuals)
        np_result = np.matmul(rand_matrix, genotype_matrix)
        grg_result = pygrgl.matmul(
            grg, rand_matrix, pygrgl.TraversalDirection.UP, by_individual=True
        )
        self.assertTrue(np.allclose(grg_result, np_result))

        # Test the (KxM)(MxN) multiplication.
        genotype_matrix = np.transpose(grg2X(grg, individual=True))
        rand_matrix = np.random.rand(K, grg.num_mutations)
        np_result = np.matmul(rand_matrix, genotype_matrix)
        grg_result = pygrgl.matmul(
            grg, rand_matrix, pygrgl.TraversalDirection.DOWN, by_individual=True
        )
        self.assertTrue(np.allclose(grg_result, np_result))

        if CLEANUP:
            os.remove(grg_filename)


if __name__ == "__main__":
    unittest.main()
