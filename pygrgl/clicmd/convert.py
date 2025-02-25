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
import subprocess
import sys
from .common import which, is_grg, is_trees, is_igd, is_bgen

GRGL = which("grgl")
GCONVERT = which("gconverter")
GINDEX = which("gindexer")


def add_options(subparser):
    subparser.add_argument(
        "input_file", help="The input file: .vcf, .vcf.gz, .bgen, .igd, or .trees"
    )
    subparser.add_argument("output_file", help="The output file: .igd or .grg")
    subparser.add_argument(
        "--binary-muts",
        "-b",
        action="store_true",
        help="Use binary mutations (don't track specific alternate alleles).",
    )
    subparser.add_argument(
        "--use-node-times",
        "-n",
        action="store_true",
        help="Mark mutations with the time of the node below them, instead of their assigned time.",
    )
    subparser.add_argument(
        "--no-simplify",
        "-l",
        action="store_true",
        help='Do not remove any "information-less" nodes/edges from the graph.',
    )
    subparser.add_argument(
        "--maintain-topo",
        "-t",
        action="store_true",
        help="Maintain all topology below mutations (at the cost of a larger graph).",
    )
    subparser.add_argument(
        "--ts-coals",
        action="store_true",
        help="Compute individual coalescence information at each GRG node (more expensive).",
    )


def convert_command(arguments):
    if is_grg(arguments.output_file):
        if not is_trees(arguments.input_file):
            print(
                f'"grg convert" only supports .trees -> .grg conversion. If you '
                ' want to construct a .grg from something else try "grg construct"'
            )
            exit(2)
        command_args = [GRGL, arguments.input_file, "-o", arguments.output_file]
        if arguments.binary_muts:
            command_args.append("--binary-muts")
        if arguments.use_node_times:
            command_args.append("--ts-node-times")
        if arguments.no_simplify:
            command_args.append("--no-simplify")
        if arguments.maintain_topo:
            command_args.append("--maintain-topo")
        if arguments.ts_coals:
            command_args.append("--ts-coals")
    elif is_igd(arguments.output_file):
        if is_trees(arguments.input_file):
            print("Can only convert .trees files to GRGs", file=sys.stderr)
            exit(2)
        if is_bgen(arguments.input_file):
            subprocess.check_call([GINDEX, arguments.input_file])
        command_args = [GCONVERT, arguments.input_file, arguments.output_file]
    else:
        assert False, "Unexpected output file type (only .grg/.igd supported)"
    subprocess.check_call(command_args)
