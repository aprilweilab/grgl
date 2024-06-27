from .clicmd import construct
from .clicmd import convert
from .clicmd import process
from .clicmd.common import which
import argparse
import subprocess
import sys

CMD_CONVERT = "convert"
CMD_CONSTRUCT = "construct"
CMD_PROCESS = "process"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", help="Print software version", action="store_true")
    subparsers = parser.add_subparsers(dest="command")
    convert_parser = subparsers.add_parser(CMD_CONVERT,
        help="Convert between tabular formats, or convert ARG (.trees) to GRG.")
    convert.add_options(convert_parser)
    construct_parser = subparsers.add_parser(CMD_CONSTRUCT,
        help="Construct a GRG from a tabular format.")
    construct.add_options(construct_parser)
    process_parser = subparsers.add_parser(CMD_PROCESS,
        help="Process a GRG to compute information from it.")
    process.add_options(process_parser)
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
        process.stat_command(args)
    else:
        print(f"Invalid command {args.command}", file=sys.stderr)
        parser.print_help()
        exit(1)

if __name__ == "__main__":
    main()