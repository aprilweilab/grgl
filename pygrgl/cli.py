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
from .clicmd import construct
from .clicmd import convert
from .clicmd import process
from .clicmd import split
from .clicmd.common import which
import argparse
import subprocess
import sys

CMD_CONVERT = "convert"
CMD_CONSTRUCT = "construct"
CMD_PROCESS = "process"
CMD_SPLIT = "split"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", help="Print software version", action="store_true")
    subparsers = parser.add_subparsers(dest="command")
    convert_parser = subparsers.add_parser(
        CMD_CONVERT,
        help="Convert between tabular formats, or convert ARG (.trees) to GRG.",
    )
    convert.add_options(convert_parser)
    construct_parser = subparsers.add_parser(
        CMD_CONSTRUCT, help="Construct a GRG from a tabular format."
    )
    construct.add_options(construct_parser)
    process_parser = subparsers.add_parser(
        CMD_PROCESS, help="Process a GRG to compute information from it."
    )
    process.add_options(process_parser)
    split_parser = subparsers.add_parser(
        CMD_SPLIT, help="Split a GRG into smaller pieces."
    )
    split.add_options(split_parser)
    args = parser.parse_args()

    if args.version:
        grgl = which("grgl")
        subprocess.check_call([grgl, "--version"])
        exit(0)
    if args.command is None:
        parser.print_help()
        exit(1)
    elif args.command == CMD_CONVERT:
        convert.convert_command(args)
    elif args.command == CMD_CONSTRUCT:
        construct.from_tabular(args)
    elif args.command == CMD_PROCESS:
        process.process_command(args)
    elif args.command == CMD_SPLIT:
        split.do_split(args)
    else:
        print(f"Invalid command {args.command}", file=sys.stderr)
        parser.print_help()
        exit(1)


if __name__ == "__main__":
    main()
