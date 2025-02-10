from typing import Optional, List
import subprocess
import time
import sys
import os

THISDIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.join(THISDIR, "..", "..")

def which(exe: str, required=False) -> Optional[str]:
    try:
        result = subprocess.check_output(["which", exe],
                                         stderr=subprocess.STDOUT).decode("utf-8").strip()
    except subprocess.CalledProcessError:
        result = None
    if result is None:
        for p in [*sys.path, REPO_ROOT]:
            p = os.path.join(p, exe)
            if os.path.isfile(p):
                result = p
                break
    if required and result is None:
        raise RuntimeError(f"Could not find executable {exe}")
    return result

def time_call(command: List[str], **kwargs) -> float:
    start_time = time.time()
    subprocess.check_call(command, **kwargs)
    return time.time() - start_time

def is_grg(filename: str) -> bool:
    return filename.endswith(".grg")

def is_trees(filename: str) -> bool:
    return filename.endswith(".trees")

def is_igd(filename: str) -> bool:
    return filename.endswith(".igd")

def is_bgen(filename: str) -> bool:
    return filename.endswith(".bgen")

