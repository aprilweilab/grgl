import subprocess
import sys
from .common import which, is_grg, is_trees, is_igd, is_bgen

GRGL = which("grgl")
GCONVERT = which("gconverter")
GINDEX = which("gindexer")

def add_options(subparser):
    subparser.add_argument("input_file",
        help="The input file: .vcf, .vcf.gz, .bgen, .igd, or .trees")
    subparser.add_argument("output_file",
        help="The output file: .igd or .grg")
    subparser.add_argument("--binary-muts", "-b", action="store_true",
        help="Use binary mutations (don't track specific alternate alleles).")
    subparser.add_argument("--use-node-times", "-n", action="store_true",
        help="Mark mutations with the time of the node below them, instead of their assigned time.")

def convert_command(arguments):
    if is_grg(arguments.output_file):
        if not is_trees(arguments.input_file):
            print(f"\"grg convert\" only supports .trees -> .grg conversion. If you "
                  " want to construct a .grg from something else try \"grg construct\"")
            exit(2)
        command_args = [GRGL, arguments.input_file,
                        "-o", arguments.output_file]
        if arguments.binary_muts:
            command_args.append("--binary-muts")
        if arguments.use_node_times:
            command_args.append("--ts-node-times")
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