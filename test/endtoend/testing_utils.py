import subprocess
import os
from typing import Optional

JOBS = 4
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(THIS_DIR, "input")


def construct_grg(
    input_file: str, output_file: Optional[str] = None, test_input: bool = True
) -> str:
    cmd = [
        "grg",
        "construct",
        "-p",
        "10",
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
