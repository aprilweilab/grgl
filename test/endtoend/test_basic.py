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

EXPECT_DIR = os.path.join(THIS_DIR, "expect")
INPUT_DIR = os.path.join(THIS_DIR, "input")
SCRIPTS_DIR = os.path.join(THIS_DIR, "..", "..", "scripts")


def get_freqs_unordered(lines: List[str]) -> Dict[Tuple[str, str, str], int]:
    fmap = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        pos, ref, alt, freq, total = line.split("\t")
        if pos == "POSITION":
            continue
        fmap[(float(pos), ref, alt)] = int(freq)
    return fmap


def construct_grg(
    input_file: str, output_file: Optional[str] = None, test_input: bool = True
) -> str:
    cmd = [
        "grg",
        "construct",
        "-p",
        "10",
        "-t",
        "2",
        "-j",
        str(JOBS),
        os.path.join(INPUT_DIR, input_file) if test_input else input_file,
    ]
    if output_file is not None:
        cmd.extend(["-o", output_file])
    else:
        output_file = os.path.basename(input_file) + ".final.grg"
    subprocess.check_call(cmd)
    return output_file


class TestGrgBasic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.grg_filename = construct_grg("test-200-samples.vcf.gz")
        cls.grg = pygrgl.load_immutable_grg(cls.grg_filename)

    def test_construct_allele_freq(self):
        # Use "grg process" to compute allele frequencies.
        af = subprocess.check_output(["grg", "process", "freq", self.grg_filename])
        saw_lines = list(filter(lambda l: l.strip(), af.decode("utf-8").split("\n")))

        # Compare computed allele frequencys against expected.
        with open(os.path.join(EXPECT_DIR, "test-200-samples.freq.txt")) as f:
            expect_lines = list(filter(lambda l: l, [line.strip() for line in f]))
        self.assertEqual(len(saw_lines), len(expect_lines))
        expect_freqs_unordered = get_freqs_unordered(expect_lines)
        # We expect the output to be identical, since the order should be 100% deterministic
        # based on (position, ref, alt) of each mutation.
        for saw, expect in zip(saw_lines, expect_lines):
            self.assertEqual(saw, expect)

        # Now compute the same values using the dot_product Python API, these we don't ensure the
        # order on (though order should match)
        grg = pygrgl.load_immutable_grg(self.grg_filename)
        input_vector = np.ones(grg.num_samples)
        api_result = pygrgl.dot_product(grg, input_vector, pygrgl.TraversalDirection.UP)
        for i in range(len(api_result)):
            m = grg.get_mutation_by_id(i)
            self.assertEqual(
                api_result[i],
                expect_freqs_unordered[(m.position, m.ref_allele, m.allele)],
            )

    def test_split(self):
        """
        Verify that splitting a GRG still captures all of the sample->mutation relationships.
        """
        # Create the GRG
        split_dir = f"{self.grg_filename}.split"
        assert not os.path.exists(split_dir), f"{split_dir} already exists; remove it"

        maxpos = 0
        for mut_id in range(self.grg.num_mutations):
            maxpos = max(maxpos, self.grg.get_mutation_by_id(mut_id).position)
        assert maxpos == 9999126
        assert self.grg.specified_bp_range == (55829, 9999127)
        assert self.grg.specified_bp_range == self.grg.bp_range

        # Split the GRG
        subprocess.check_output(
            ["grg", "split", "-j", str(4), self.grg_filename, "-s", str(100000)]
        )
        saw_freqs = {}
        for grg_part in glob.glob(f"{split_dir}/*.grg"):
            af = subprocess.check_output(["grg", "process", "freq", grg_part])
            saw_freqs.update(get_freqs_unordered(af.decode("utf-8").split("\n")))

        # Compare computed allele frequencys against expected.
        with open(os.path.join(EXPECT_DIR, "test-200-samples.freq.txt")) as f:
            expect_lines = [line for line in f]
        expect_freqs = get_freqs_unordered(expect_lines)
        saw_keys = set(saw_freqs.keys())
        expect_keys = set(expect_freqs.keys())
        self.assertEqual(saw_keys, expect_keys)
        for k in saw_keys:
            self.assertEqual(saw_freqs[k], expect_freqs[k])

    def test_coal_counts(self):
        zyg_info = subprocess.check_output(
            ["grg", "process", "zygosity", self.grg_filename]
        )
        saw_lines = list(
            filter(lambda l: l.strip(), zyg_info.decode("utf-8").split("\n"))
        )

        # Compare computed zygosity information against expected. The expected info was checked against
        # the real values using this pyigd script:
        # https://gist.github.com/dcdehaas/300ca4c97e65362ba6040e3bdcc2d0b9
        with open(os.path.join(EXPECT_DIR, "test-200-samples.zygosity.txt")) as f:
            expect_lines = list(filter(lambda l: l, [line.strip() for line in f]))
        self.assertEqual(len(saw_lines), len(expect_lines))
        for saw, expect in zip(saw_lines, expect_lines):
            self.assertEqual(saw, expect)

    def test_igd_construct(self):
        igd_filename = "test-200-samples.igd"
        cmd = [
            "igdtools",
            os.path.join(INPUT_DIR, "test-200-samples.vcf.gz"),
            "-o",
            igd_filename,
        ]
        subprocess.check_call(cmd)
        grg_filename = construct_grg(igd_filename, test_input=False)
        grg = pygrgl.load_immutable_grg(grg_filename)
        self.assertEqual(grg.num_samples, self.grg.num_samples)

        shape = (1, grg.num_samples)
        from_igd = pygrgl.matmul(grg, np.ones(shape), pygrgl.TraversalDirection.UP)
        from_vcf = pygrgl.matmul(self.grg, np.ones(shape), pygrgl.TraversalDirection.UP)

        self.assertTrue(np.allclose(from_igd, from_vcf))
        if CLEANUP:
            os.remove(grg_filename)

    @classmethod
    def tearDownClass(cls):
        if CLEANUP:
            os.remove(cls.grg_filename)


if __name__ == "__main__":
    unittest.main()
