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

from pygrgl import load_mutable_grg, save_grg

GRG_MERGE = which("grg-merge", required=True)


def add_options(subparser):
    subparser.add_argument("out_file", help="The output GRG file to create")
    subparser.add_argument("grg_file", nargs="*", help="The input GRG files")
    subparser.add_argument(
        "-C",
        "--no-combine",
        action="store_true",
        help="Don't combine nodes, keep all the nodes/edges from all graphs.",
    )
    subparser.add_argument(
        "-s",
        "--spec-file",
        help="A specification file listing all the GRGs to be merged. Each line of the file contains"
        "a single filename. Each filename can be suffixed with ':<position_offset>` to add the value "
        "<position_offset> to every mutation position in the GRG prior to merging.",
    )
    subparser.add_argument(
        "--use-samples",
        action="store_true",
        help="Use node equivalence via samples beneath each node, instead of just children. "
        "This reuses less hierarchy in the graph, but can produce a smaller graph. "
        "Slower and more RAM intensive.",
    )
    subparser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output.",
    )


def merge_command(args):
    if len(args.grg_file) == 0 and args.spec_file is None:
        raise RuntimeError(
            "You must provide either a list of input GRGs on the command line,"
            " or --spec-file with a file containing a list of GRGs"
        )
    input_files = []
    adjust_positions = []
    if args.spec_file is not None:
        with open(args.spec_file, "r") as f:
            spec = list(filter(lambda s: len(s) > 0, map(str.strip, f)))
        for s in spec:
            s = s.split(":")
            input_files.append(s[0])
            if len(s) > 1:
                adjust_positions.append(int(s[1]))
        if adjust_positions and (len(input_files) != len(adjust_positions)):
            raise RuntimeError(
                "Must provide adjustment positions for all files, if any are provided"
            )
        if adjust_positions and adjust_positions[0] != 0:
            raise RuntimeError(
                "The first to-be-merged graph must have an adjustment position of 0"
            )
    else:
        input_files = args.grg_file
    if len(input_files) <= 1:
        raise RuntimeError("Cannot merge fewer than two graphs together.")
    target = load_mutable_grg(input_files[0])
    target.merge(
        input_files[1:],
        combine_nodes=not args.no_combine,
        use_sample_sets=True if args.use_samples else False,
        verbose=args.verbose,
        position_adjust=adjust_positions[1:],
    )
    save_grg(target, args.out_file)
