# Genotype Representation Graph Library (GRGL)
# Copyright (C) 2025 April Wei
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
from .common import which
import subprocess

GRG_MERGE = which("grg-merge", required=True)


def add_options(subparser):
    subparser.add_argument("out_file", help="The output GRG file to create")
    subparser.add_argument("grg_file", nargs="+", help="The input GRG files")
    subparser.add_argument(
        "--use-samples",
        action="store_true",
        help="Use node equivalence via samples beneath each node, instead of just children. "
        "This reuses less hierarchy in the graph, but can produce a smaller graph. "
        "Slower and more RAM intensive.",
    )


def merge_command(args):
    command_args = [GRG_MERGE]
    if args.use_samples:
        command_args.append("--use-samples")
    command_args.append(args.out_file)
    command_args.extend(args.grg_file)
    subprocess.check_call(command_args)
