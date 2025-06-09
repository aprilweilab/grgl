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
from .common import which
from enum import Enum
import subprocess
import sys

GRGP = which("grgp", required=True)


class Statistic(Enum):
    GRAPH_STATS = "stats"
    ALLELE_FREQ = "freq"
    GWAS = "gwas"
    ZYGOSITY = "zygosity"

    def __str__(self):
        return self.value


def add_options(subparser):
    subparser.add_argument(
        "operation",
        type=Statistic,
        choices=list(Statistic),
        help="The operation to perform on the GRG file",
    )
    subparser.add_argument("grg_file", help="The input GRG file")
    subparser.add_argument(
        "-p", "--phenotype", help="The phenotype file (for GWAS only)"
    )


def process_command(arguments):
    command_args = [GRGP, arguments.grg_file]
    if arguments.operation == Statistic.GRAPH_STATS:
        command_args.append("-s")
    elif arguments.operation == Statistic.ALLELE_FREQ:
        command_args.append("--freq")
    elif arguments.operation == Statistic.GWAS:
        if arguments.phenotype is None:
            print("You must pass a file via --phenotype", file=sys.stderr)
            exit(1)
        command_args.extend(["--association-study", "--phenotype", arguments.phenotype])
    elif arguments.operation == Statistic.ZYGOSITY:
        command_args.append("--zygosity-info")
    else:
        assert False, "Unexpect operation: " + str(arguments.operation)
    subprocess.check_call(command_args)
