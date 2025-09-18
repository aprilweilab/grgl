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
import subprocess
import sys
import tqdm
from multiprocessing import Pool
from typing import List, Tuple
from .common import which, time_call, round_up_to


HAP_SEG_LENGTH = 128


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
        default=None,
        help="The number of parts to split the sequence into; defaults to auto-calculate.",
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
        default=None,
        help="Number of trees to use during shape construction. Defaults to auto-calculate.",
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
    subparser.add_argument(
        "--ignore-missing",
        action="store_true",
        help="Do not store missing data in the GRG, pretend missing alleles are the reference allele.",
    )
    # There are 9 compression levels:
    # 1: same as 2
    # 2: single-pass, use FASTER2 to determine # of trees, direct-map muts with count less
    #    than 10.
    # 3: single-pass, use FASTER2 to determine # of trees
    # 4: single-pass, use FASTER1 to determine # of trees
    # 5: single-pass, use OPTIMAL to determine # of trees
    # 6: multi-pass, use OPTIMAL to determine # of trees, direct-map muts with count less
    #    than 10.
    # 7: multi-pass, use OPTIMAL to determine # of trees.
    # 8: same as 7
    # 9: same as 7
    level = subparser.add_mutually_exclusive_group()
    level.add_argument(
        "--level1",
        "-1",
        action="store_true",
        help="The fastest (least compressive) level. The default.",
    )
    level.add_argument("--level2", "-2", action="store_true", help="2nd fastest level.")
    level.add_argument("--level3", "-3", action="store_true", help="3rd fastest level.")
    level.add_argument("--level4", "-4", action="store_true", help="4th fastest level.")
    level.add_argument("--level5", "-5", action="store_true", help="5th fastest level.")
    level.add_argument("--level6", "-6", action="store_true", help="6th fastest level.")
    level.add_argument("--level7", "-7", action="store_true", help="7th fastest level.")
    level.add_argument("--level8", "-8", action="store_true", help="8th fastest level.")
    level.add_argument(
        "--level9",
        "-9",
        action="store_true",
        help="The slowest (most compressive) level.",
    )
    subparser.add_argument(
        "--force",
        action="store_true",
        help="Ignore any warning conditions that cause execution to stop.",
    )


grgl_exe = which("grgl")
grg_merge_exe = which("grg-merge")


# The advantage of a lot of parts is that it can balance the parallelism. The disadvantage,
# for small inputs, is that it can waste a lot of time reading the input file and determining
# that there are no variants in a particular range.
def compute_parts(input_file: str, threads: int):
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    cmd = [grgl_exe, "--count-variants", input_file]
    try:
        output = subprocess.check_output(cmd).decode("utf-8").split("\n")[0]
        num_variants = int(output)
        min_var_per_part = HAP_SEG_LENGTH
        max_parts = num_variants // min_var_per_part
    except (ValueError, subprocess.CalledProcessError):
        print(
            f"Could not count number of variants in {input_file}. Using the default of 100 (use --parts to override).",
            file=sys.stderr,
        )
        # 100 is obviously not optimal for small datasets, but oh well!
        max_parts = 100
    best_parts = max(threads, round_up_to(100, threads))
    return min(max_parts, best_parts)


def log_v(msg, verbose):
    if verbose:
        print(msg, file=sys.stderr)


def log_time(name, time_val, verbose):
    if verbose:
        print(f"{name}={time_val}", file=sys.stderr)


def out_filename(output_file, part):
    base_name = os.path.basename(output_file)
    return f"{base_name}.part{part}.grg"


def get_default_range_suffix(input_file: str) -> str:
    # For IGD or BGEN, split the file by variant counts, because indexing is very easy and
    # the resulting GRG should be better and faster to create.
    suffix = ""
    if input_file.endswith(".igd") or input_file.endswith(".bgen"):
        suffix = "v"
    return suffix


def build_shape(
    range_triple, args, auto_args: List[str], input_file: str, output_file: str
):
    part, lower, upper = range_triple
    assert lower < upper
    suffix = get_default_range_suffix(input_file)

    assert args.trees is not None
    command = [grgl_exe, input_file, "--trees", str(args.trees)] + auto_args
    if args.maf_flip:
        command.append("--maf-flip")
    if args.population_ids:
        command.extend(["--population-ids", args.population_ids])
    if args.no_indiv_ids:
        command.append("--no-indiv-ids")
    if args.verbose:
        command.extend(["--verbose", "-s"])
    if args.ignore_missing:
        command.append("--ignore-missing")
    if args.force:
        command.append("--force")
    shape_filename = out_filename(output_file, part)
    command.extend(
        [
            "-r",
            f"{lower}:{upper}{suffix}",
            "-o",
            shape_filename,
        ]
    )
    log_v(command, args.verbose)
    tb_time = time_call(command, stdout=sys.stdout)
    log_time("TREE_BUILD_TIME", tb_time, args.verbose)
    sys.stdout.flush()
    return shape_filename


def build_grg(
    range_triple: Tuple[int, float, float],
    args,
    auto_args: List[str],
    input_file: str,
    output_file: str,
    do_map_muts: bool,
):
    shape_grg = build_shape(range_triple, args, auto_args, input_file, output_file)
    part, lower, upper = range_triple
    # When we do --fast, no need to call MapMutations
    if do_map_muts:
        command = [grgl_exe, shape_grg]
        suffix = get_default_range_suffix(input_file)
        if args.maf_flip:
            command.append("--maf-flip")
        if args.verbose:
            command.extend(["--verbose", "-s"])
        if args.ignore_missing:
            command.append("--ignore-missing")
        if args.force:
            command.append("--force")
        command.extend(
            [
                "-r",
                f"{lower}:{upper}{suffix}",
                "-m",
                input_file,
                "-o",
                out_filename(output_file, part),
            ]
        )
        if args.binary_muts:
            command.append("-b")
        log_v(command, args.verbose)
        map_time = time_call(command, stdout=sys.stdout)
        log_time("MAP_MUTS_TIME", map_time, args.verbose)
    sys.stdout.flush()


def star_build_grg(args):
    return build_grg(*args)


def from_tabular(args):
    if args.range is not None:
        if ":" not in args.range:
            raise RuntimeError('--range must be specified as "lower:upper"')
        if args.parts is not None:
            raise RuntimeError(f"WARNING: Cannot specify both --range and --parts.")

    def verify_file(fn):
        if not os.path.isfile(fn):
            raise RuntimeError(f"File not found: {fn}")

    verify_file(args.input_file)

    if grgl_exe is None:
        raise RuntimeError("Could not find 'grgl' executable; please add to your PATH")
    log_v(f"Using grgl at: {grgl_exe}", args.verbose)
    if grg_merge_exe is None:
        raise RuntimeError(
            "Could not find 'grg-merge' executable; please add to your PATH"
        )
    log_v(f"Using grg-merge at: {grg_merge_exe}", args.verbose)

    lz4_exe = which("lz4")
    if args.input_file.endswith(".lz4"):
        assert lz4_exe is not None
        input_file = args.input_file[:-4]
        command = [lz4_exe, "-f", args.input_file, input_file]
        unlz4_elapsed = time_call(command, stdout=sys.stdout)
        print(f"Decompress took {unlz4_elapsed} seconds", file=sys.stderr)
        print("Decompressed size:", file=sys.stderr)
        try:
            subprocess.check_call(["du", "-hs", input_file])
        except subprocess.CalledProcessError as e:
            print(f"Failed to get size of {input_file}!", file=sys.stderr)
            print(e, file=sys.stderr)
            exit(2)
    else:
        input_file = args.input_file

    base_name = os.path.basename(input_file)
    if args.out_file is not None:
        final_filename = args.out_file
    else:
        final_filename = f"{base_name}.final.grg"

    if args.parts is None:
        args.parts = compute_parts(input_file, args.jobs)
    print(f"Processing input file in {args.parts} parts.", file=sys.stderr)

    do_map_mutations = False
    # 3: single-pass, use FASTER2 to determine # of trees
    if args.level3:
        auto_args = ["--trees"]
        auto_tree = "faster2"
    # 4: single-pass, use FASTER1 to determine # of trees
    elif args.level4:
        auto_args = []
        auto_tree = "faster1"
    # 5: single-pass, use OPTIMAL to determine # of trees
    elif args.level5:
        auto_args = []
        auto_tree = "optimal"
    # 6: multi-pass, use OPTIMAL to determine # of trees, direct-map muts with count less
    #    than 10.
    elif args.level6:
        auto_args = ["--lf-no-tree", str(10), "--no-tree-map", "--no-simplify"]
        auto_tree = "optimal"
        do_map_mutations = True
    # 7: multi-pass, use OPTIMAL to determine # of trees.
    elif args.level7 or args.level8 or args.level9:
        auto_args = ["--no-tree-map", "--no-simplify"]
        auto_tree = "optimal"
        do_map_mutations = True
    # level1/2 and default are the same
    # 2: single-pass, use FASTER2 to determine # of trees, direct-map muts with count less
    #    than 10.
    else:
        auto_args = ["--lf-no-tree", str(10)]
        auto_tree = "faster2"

    if args.trees is None:
        print(f"Auto-calculating number of trees per part.", file=sys.stderr)
        args.trees = auto_tree
    else:
        print(f"Using {args.trees} trees per part.", file=sys.stderr)

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

    # Compute the separate GRGs in parallel.
    if len(ranges) == 1:
        build_grg(
            ranges[0], args, auto_args, input_file, final_filename, do_map_mutations
        )
    else:
        print("Converting segments of input data to graphs", file=sys.stderr)
        with Pool(args.jobs) as pool:
            list(
                tqdm.tqdm(
                    pool.imap_unordered(
                        star_build_grg,
                        [
                            (
                                r,
                                args,
                                auto_args,
                                input_file,
                                final_filename,
                                do_map_mutations,
                            )
                            for r in ranges
                        ],
                    ),
                    total=len(ranges),
                )
            )

    if (not args.no_file_cleanup) and input_file != args.input_file:
        log_v(
            f"Removing uncompressed input file {input_file} (done with it)",
            args.verbose,
        )
        os.remove(input_file)

    if not args.no_merge:
        # Now merge them pairwise.
        print("Merging...", file=sys.stderr)

        command = [
            grg_merge_exe,
            "-s",
            final_filename,
        ]
        command.extend(
            map(lambda part: out_filename(final_filename, part), range(0, args.parts))
        )
        log_v(command, args.verbose)
        final_merge_time = time_call(command)
        log_time("FINAL_MERGE_TIME", final_merge_time, args.verbose)

        if not args.no_file_cleanup:
            for part in range(0, args.parts):
                os.remove(out_filename(final_filename, part))


def main():
    parser = argparse.ArgumentParser(
        description="Construct a GRG from a VCF/IGD/BGEN file."
    )
    add_options(parser)
    args = parser.parse_args()
    from_tabular(args)


if __name__ == "__main__":
    main()
