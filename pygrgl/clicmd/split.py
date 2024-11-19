from .common import which, time_call

def add_options(subparser):
    subparser.add_argument("input_file", help="The input GRG file")
    subparser.add_argument("size_per_grg", type=float,
        help="The amount of BP or cM per GRG. If --rec-map is specified then this is cM, otherwise BP")
    subparser.add_argument("--jobs", "-j", type=int, default=1,
        help="Number of jobs (threads/cores) to use. Defaults to 1.")
    subparser.add_argument("--rec-map", "-r", default=None,
        help="Use the given HapMap-style recombination map and interpret the size_per_grg as cM.")

grgl_exe = which("grgl")

def do_split(args):
    if grgl_exe is None:
        raise RuntimeError("Could not find 'grgl' executable; please add to your PATH")
    split_arg = f"{args.size_per_grg}"
    if args.rec_map is not None:
        split_arg = f"{args.rec_map}:{split_arg}"
    cmd = [grgl_exe, args.input_file, "--split", split_arg, "-j", str(args.jobs)]
    time_call(cmd)