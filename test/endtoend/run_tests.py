import unittest
import subprocess
import os
from typing import List, Dict, Tuple

JOBS = 4
CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

EXPECT_DIR = os.path.join(THIS_DIR, "expect")
INPUT_DIR = os.path.join(THIS_DIR, "input")
SCRIPTS_DIR = os.path.join(THIS_DIR, "..", "..", "scripts")

def get_freqs(lines: List[str]) -> Dict[Tuple[str, str, str], int]:
    fmap = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        pos, ref, alt, freq, total = line.split("\t")
        if pos == "POSITION":
            continue
        fmap[(pos, ref, alt)] = int(freq)
    return fmap


class TestGrgConstruction(unittest.TestCase):
    def construct_grg(self, input_file:str) -> str:
        subprocess.check_call(["grg-construct", "--no-maf-flip", "-p", "10", "-t", "2", "-j", str(JOBS),
                               os.path.join(INPUT_DIR, input_file)])
        return os.path.basename(input_file) + ".final.grg"

    def test_construct_allele_freq(self):
        grg_filename = self.construct_grg("test-200-samples.vcf.gz")
        af = subprocess.check_output(["grgl", "--freq", grg_filename])
        saw_freqs = get_freqs(af.decode("utf-8").split("\n"))
        with open(os.path.join(EXPECT_DIR, "test-200-samples.freq.txt")) as f:
            expect_lines = [line for line in f]
        expect_freqs = get_freqs(expect_lines)

        saw_keys = set(saw_freqs.keys())
        expect_keys = set(expect_freqs.keys())
        self.assertEqual(saw_keys, expect_keys)
        for k in saw_keys:
            self.assertEqual(saw_freqs[k], expect_freqs[k])
        if CLEANUP:
            os.remove(grg_filename)

if __name__ == '__main__':
    unittest.main()
