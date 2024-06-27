from enum import Enum
import subprocess
from .common import which

GRGP = which("grgp", required=True)

class Statistic(Enum):
    GRAPH_STATS = "stats"
    ALLELE_FREQ = "freq"
    GWAS = "gwas"

    def __str__(self):
        return self.value

def add_options(subparser):
    subparser.add_argument("operation", type=Statistic, choices=list(Statistic),
        help="The operation to perform on the GRG file")
    subparser.add_argument("grg_file", help="The input GRG file")

def stat_command(arguments):
    command_args = [GRGP, arguments.grg_file]
    if arguments.operation == Statistic.GRAPH_STATS:
        command_args.append("-s")
    elif arguments.operation == Statistic.ALLELE_FREQ:
        command_args.append("--freq")
    elif arguments.operation == Statistic.GWAS:
        assert False, "Not yet implemented"
    else:
        assert False, ("Unexpect operation: " + str(arguments.operation))
    subprocess.check_call(command_args)
