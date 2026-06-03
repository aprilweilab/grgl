import os
import shutil
import subprocess
import sys
import tempfile
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
INPUT_DIR = os.path.join(THIS_DIR, "input")
sys.path.insert(0, REPO_ROOT)

import numpy as np
import pygrgl
from pygrgl.clicmd.common import which


class TestPolarize(unittest.TestCase):
    def test_multiallelic_snv_site_from_vcf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_vcf = os.path.join(INPUT_DIR, "multi.vcf")
            ancestral_fasta = os.path.join(tmpdir, "ancestor.fa")
            input_grg = os.path.join(tmpdir, "multi.grg")
            output_grg = os.path.join(tmpdir, "multi.polar.grg")

            seq = ["N"] * 1110700
            seq[14370 - 1] = "G"
            seq[17330 - 1] = "A"
            seq[1110696 - 1] = "T"
            with open(ancestral_fasta, "w") as out:
                out.write(">20\n")
                sequence = "".join(seq)
                for start in range(0, len(sequence), 80):
                    out.write(sequence[start : start + 80] + "\n")

            subprocess.check_call(
                [
                    which("grgl", required=True),
                    input_vcf,
                    "--ignore-missing",
                    "--force",
                    "-o",
                    input_grg,
                ]
            )
            grg_command = ["grg"] if shutil.which("grg") else [sys.executable, "-m", "pygrgl.cli"]
            subprocess.check_call(
                grg_command
                + [
                    "polarize",
                    input_grg,
                    ancestral_fasta,
                    "-o",
                    output_grg,
                    "--map-batch-size",
                    "2",
                ]
            )

            grg = pygrgl.load_mutable_grg(output_grg, load_up_edges=True)
            values = np.eye(grg.num_mutations, dtype=np.int32)
            sample_matrix = pygrgl.matmul(grg, values, pygrgl.TraversalDirection.DOWN)

            observed = {}
            for mut_id in range(grg.num_mutations):
                mutation = grg.get_mutation_by_id(mut_id)
                observed[(int(mutation.position), mutation.ref_allele, mutation.allele)] = tuple(
                    int(sample_id) for sample_id in np.flatnonzero(sample_matrix[mut_id] > 0)
                )

            self.assertEqual(
                observed,
                {
                    (17330, "A", "T"): (0, 1, 2, 4, 5),
                    (1110696, "T", "A"): (3,),
                    (1110696, "T", "G"): (0,),
                },
            )


if __name__ == "__main__":
    unittest.main()
