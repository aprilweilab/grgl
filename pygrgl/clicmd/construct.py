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

# Command-line tool for constructing a GRG. Calls the "grgl" executable under the hood.
import argparse
import os
import sys
import subprocess
from multiprocessing import Pool
from .common import which, time_call


def add_options(subparser):
    subparser.add_argument(
        "input_file", help="The input file: .vcf, .vcf.gz, .igd, or .bgen"
    )
    subparser.add_argument(
        "--range",
        "-r",
        type=str,
        default=None,
        help="Restrict to the given range. Can be absolute (in base-pairs) or relative (0.0 to 1.0).",
    )
    subparser.add_argument(
        "--parts",
        "-p",
        type=int,
        default=8,
        help="The number of parts to split the sequence into; defaults to 8",
    )
    subparser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of jobs (threads/cores) to use. Defaults to 1.",
    )
    subparser.add_argument(
        "--trees",
        "-t",
        type=int,
        default=1,
        help="Number of trees to use during shape construction. Defaults to 1.",
    )
    subparser.add_argument(
        "--binary-muts",
        "-b",
        action="store_true",
        help="Use binary mutations (don't track specific alternate alleles).",
    )
    subparser.add_argument(
        "--no-file-cleanup",
        "-c",
        action="store_true",
        help="Do not cleanup intermediate files (for debugging, e.g.).",
    )
    subparser.add_argument(
        "--maf-flip",
        action="store_true",
        help="Switch the reference allele with the major allele when they differ",
    )
    subparser.add_argument(
        "--shape-lf-filter",
        "-f",
        type=float,
        default=10.0,
        help="During shape construction ignore mutations with counts less than this."
        "If value is <1.0 then it is treated as a frequency. Defaults to 10 (count).",
    )
    subparser.add_argument(
        "--population-ids",
        default=None,
        help='Format: "filename:fieldname". Read population ids from the given '
        "tab-separate file, using the given fieldname.",
    )
    subparser.add_argument(
        "--bs-triplet",
        default=0,
        help="Run the triplet algorithm for this many iterations in BuildShape",
    )
    subparser.add_argument(
        "--out-file",
        "-o",
        default=None,
        help="Specify an output file instead of using the default name.",
    )
    subparser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output, including timing information.",
    )
    subparser.add_argument(
        "--no-merge",
        action="store_true",
        help='Do not merge the resulting GRGs (so if you specified "-p C" there will be C GRGs).',
    )
    subparser.add_argument(
        "--no-indiv-ids",
        action="store_true",
        help="Do not storage individual string identifiers in the GRG.",
    )


grgl_exe = which("grgl")
grg_merge_exe = which("grg-merge")


def log_time(name, time_val, verbose):
    if verbose:
        print(f"{name}={time_val}")


def out_filename_tree(input_file, part, tnum):
    base_name = os.path.basename(input_file)
    return f"{base_name}.part{part}.tree{tnum}.grg"


def out_filename(input_file, part):
    base_name = os.path.basename(input_file)
    return f"{base_name}.part{part}.grg"


def build_shape(range_triple, args, input_file):
    part, lower, upper = range_triple
    assert lower < upper
    span = upper - lower
    pspans = span / args.trees

    for tnum in range(args.trees):
        base = lower + (tnum * pspans)
        command = [grgl_exe, input_file]
        if args.maf_flip:
            command.append("--maf-flip")
        if args.population_ids:
            command.extend(["--population-ids", args.population_ids])
        if args.bs_triplet:
            command.extend(["--bs-triplet", args.bs_triplet])
        if args.no_indiv_ids:
            command.append("--no-indiv-ids")
        command.extend(["--lf-filter", str(args.shape_lf_filter)])
        command.extend(
            [
                "-l",
                "-s",
                "-r",
                f"{base}:{base+pspans}",
                "-o",
                out_filename_tree(input_file, part, tnum),
            ]
        )
        print(command)
        tb_time = time_call(command, stdout=sys.stdout)
        log_time("TREE_BUILD_TIME", tb_time, args.verbose)
    base_name = os.path.basename(input_file)
    shape_filename = f"{base_name}.part{part}.grg"
    command = [
        grg_merge_exe,
        "-l",
        "-s",
        shape_filename,
    ]
    command.extend(
        map(
            lambda tnum: out_filename_tree(input_file, part, tnum), range(0, args.trees)
        )
    )
    print(command)
    tm_time = time_call(command)
    log_time("TREE_MERGE_TIME", tm_time, args.verbose)
    if not args.no_file_cleanup:
        for tnum in range(args.trees):
            os.remove(out_filename_tree(input_file, part, tnum))
    sys.stdout.flush()
    return shape_filename


def build_grg(range_triple, args, input_file):
    shape_grg = build_shape(range_triple, args, input_file)
    part, lower, upper = range_triple
    command = [grgl_exe, shape_grg]
    if args.maf_flip:
        command.append("--maf-flip")
    command.extend(
        [
            "-s",
            "-r",
            f"{lower}:{upper}",
            "-m",
            input_file,
            "-o",
            out_filename(input_file, part),
        ]
    )
    if args.binary_muts:
        command.append("-b")
    print(command)
    map_time = time_call(command, stdout=sys.stdout)
    log_time("MAP_MUTS_TIME", map_time, args.verbose)
    sys.stdout.flush()


def _build_grg(args):
    build_grg(*args)


def from_tabular(args):
    if args.range is not None:
        if ":" not in args.range:
            raise RuntimeError('--range must be specified as "lower:upper"')
        if args.parts != 1:
            print(
                f"WARNING: Cannot specify both --range and --parts. Changing --parts from {args.parts} to 1.",
                file=sys.stderr,
            )
            args.parts = 1

    def verify_file(fn):
        if not os.path.isfile(fn):
            raise RuntimeError(f"File not found: {fn}")

    verify_file(args.input_file)

    if grgl_exe is None:
        raise RuntimeError("Could not find 'grgl' executable; please add to your PATH")
    if args.verbose:
        print(f"Using grgl at: {grgl_exe}")
    if grg_merge_exe is None:
        raise RuntimeError(
            "Could not find 'grg-merge' executable; please add to your PATH"
        )
    if args.verbose:
        print(f"Using grg-merge at: {grg_merge_exe}")

    lz4_exe = which("lz4")
    if args.input_file.endswith(".lz4"):
        assert lz4_exe is not None
        input_file = args.input_file[:-4]
        command = [lz4_exe, "-f", args.input_file, input_file]
        unlz4_elapsed = time_call(command, stdout=sys.stdout)
        print(f"Decompress took {unlz4_elapsed} seconds")
        print("Decompressed size:")
        try:
            subprocess.check_call(["du", "-hs", input_file])
        except subprocess.CalledProcessError as e:
            print(f"Failed to get size of {input_file}!")
            print(e)
    else:
        input_file = args.input_file

    # We used normalized ranges to avoid having to know the sequence length (or the
    # recombination rates, when using non-absolute distances between variants)
    if args.range is None:
        edges = []
        for i in range(args.parts + 1):
            edges.append(i / args.parts)
        ranges = []
        for i in range(1, len(edges)):
            ranges.append((i - 1, edges[i - 1], edges[i]))
    else:
        l = float(args.range.split(":")[0])
        r = float(args.range.split(":")[1])
        ranges = [(0, l, r)]

    base_name = os.path.basename(input_file)

    # Compute the separate GRGs in parallel.
    if len(ranges) == 1:
        build_grg(ranges[0], args, input_file)
    else:
        with Pool(args.jobs) as pool:
            pool.map(_build_grg, [(r, args, input_file) for r in ranges])

    if (not args.no_file_cleanup) and input_file != args.input_file:
        print(f"Removing uncompressed input file {input_file} (done with it)")
        os.remove(input_file)

    if not args.no_merge:
        # Now merge them pairwise.
        print("Merging...")
        if args.out_file is not None:
            final_filename = args.out_file
        else:
            final_filename = f"{base_name}.final.grg"

        command = [
            grg_merge_exe,
            "-s",
            final_filename,
        ]
        command.extend(
            map(lambda part: out_filename(input_file, part), range(0, args.parts))
        )
        print(command)
        final_merge_time = time_call(command)
        log_time("FINAL_MERGE_TIME", final_merge_time, args.verbose)

        if not args.no_file_cleanup:
            for part in range(0, args.parts):
                os.remove(out_filename(input_file, part))


def main():
    parser = argparse.ArgumentParser(description="Construct a GRG from a VCF file.")
    add_options(parser)
    args = parser.parse_args()


if __name__ == "__main__":
    main()
