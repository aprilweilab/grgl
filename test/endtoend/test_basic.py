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


class TestGrgConstruction(unittest.TestCase):
    def test_construct_allele_freq(self):
        # Create the GRG
        grg_filename = construct_grg("test-200-samples.vcf.gz")

        # Use "grg process" to compute allele frequencies.
        af = subprocess.check_output(["grg", "process", "freq", grg_filename])
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
        grg = pygrgl.load_immutable_grg(grg_filename)
        input_vector = np.ones(grg.num_samples)
        api_result = pygrgl.dot_product(grg, input_vector, pygrgl.TraversalDirection.UP)
        for i in range(len(api_result)):
            m = grg.get_mutation_by_id(i)
            self.assertEqual(
                api_result[i],
                expect_freqs_unordered[(m.position, m.ref_allele, m.allele)],
            )

        if CLEANUP:
            os.remove(grg_filename)

    def test_split(self):
        """
        Verify that splitting a GRG still captures all of the sample->mutation relationships.
        """
        # Create the GRG
        grg_filename = construct_grg("test-200-samples.vcf.gz", "test.for_split.grg")
        split_dir = f"{grg_filename}.split"
        assert not os.path.exists(split_dir), f"{split_dir} already exists; remove it"

        check_grg = pygrgl.load_immutable_grg(grg_filename)
        maxpos = 0
        for mut_id in range(check_grg.num_mutations):
            maxpos = max(maxpos, check_grg.get_mutation_by_id(mut_id).position)
        assert maxpos == 9999126
        assert check_grg.specified_bp_range == (55829, 9999127)
        assert check_grg.specified_bp_range == check_grg.bp_range

        # Split the GRG
        subprocess.check_output(
            ["grg", "split", "-j", str(4), grg_filename, "-s", str(100000)]
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
        if CLEANUP:
            os.remove(grg_filename)


class TestCoalescence(unittest.TestCase):
    def test_coal_counts(self):
        grg_filename = construct_grg("test-200-samples.vcf.gz", "test.coal_counts.grg")
        zyg_info = subprocess.check_output(["grg", "process", "zygosity", grg_filename])
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

        if CLEANUP:
            os.remove(grg_filename)


if __name__ == "__main__":
    unittest.main()
