# Genotype Representation Graph Library (GRGL)
# Copyright (C) 2024 April Wei
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this program.  If not, see <https://www.gnu.org/licenses/>.
from typing import Optional, List
import subprocess
import time
import sys
import os

THISDIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.join(THISDIR, "..", "..")


def which(exe: str, required=False) -> Optional[str]:
    try:
        result = (
            subprocess.check_output(["which", exe], stderr=subprocess.STDOUT)
            .decode("utf-8")
            .strip()
        )
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
