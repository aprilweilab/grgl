#!/usr/bin/env python
import argparse
import os
import sys
import time
import subprocess
from multiprocessing import Pool
from typing import Optional, List, Tuple
from dataclasses import dataclass

def which(exe: str) -> Optional[str]:
    try:
        result = subprocess.check_output(["which", exe]).decode("utf-8").strip()
    except subprocess.CalledProcessError:
        result = None
    if result is None:
        for p in sys.path:
            p = os.path.join(p, exe)
            if os.path.isfile(p):
                result = p
                break
    return result

grgl_exe = which("grgl")
if grgl_exe is None:
    raise RuntimeError("Could not find 'grgl' executable; please add to your PATH")
grg_merge_exe = which("grg-merge")
if grg_merge_exe is None:
    raise RuntimeError("Could not find 'grg-merge' executable; please add to your PATH")

@dataclass
class Range:
    start: float
    end: float

@dataclass
class TreeSpec:
    build: Range

def out_filename_tree(input_file: str, tree: TreeSpec):
    base_name = os.path.basename(input_file)
    tree_spec = f"{tree.build.start}_{tree.build.end}"
    return f"{base_name}.tree_{tree_spec}.grg"

def build_shape(input_file: str, tree: TreeSpec) -> str:
    str_range = f"{tree.build.start}:{tree.build.end}"
    out_file = out_filename_tree(input_file, tree)
    command = [grgl_exe, input_file, "-l", "-s", "-r", str_range, "-o", out_file]
    subprocess.check_call(command, stdout=sys.stdout)
    return out_file

def build_shape_wrapper(args):
    return build_shape(*args)

def merge_trees(input_file: str, trees: List[TreeSpec], part: int, cleanup_files: bool) -> Tuple[str, Range]:
    smallest_pos = 2**32
    largest_pos = 0
    base_name = os.path.basename(input_file)
    shape_filename = f"{base_name}.part{part}.grg"
    command = [grg_merge_exe, "-l", "-s", shape_filename, ]
    for tree in trees:
        command.append(out_filename_tree(input_file, tree))
        if tree.build.start < smallest_pos:
            smallest_pos = tree.build.start
        if tree.build.end > largest_pos:
            largest_pos = tree.build.end
    subprocess.check_call(command, stdout=sys.stdout)
    if cleanup_files:
        for tree in trees:
            os.remove(out_filename_tree(input_file, tree))
    return (shape_filename, Range(start=smallest_pos, end=largest_pos))

def map_part(input_file: str, part_file: str, map_range: Range, binary_muts: bool) -> str:
    command = [grgl_exe, part_file, "-s", "-r", f"{map_range.start}:{map_range.end}",
               "-m", input_file, "-o", part_file]
    if binary_muts:
        command.append("-b")
    print(command)
    return subprocess.check_output(command).decode("utf-8")

SEG_SIZE = 500_000 # 500Kbp segment size
TRIALS = reversed([1, 4, 8, 16])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a GRG shape from a VCF/IGD file.")
    parser.add_argument("input_file", help="The input file (.vcf, .vcf.gz, .igd, .bgen)")
    parser.add_argument("start_offset", help="The genome position to start at", type=int)
    parser.add_argument("--jobs", "-j", help="The number of jobs (threads/cores) to use", type=int,
                        default=2)
    args = parser.parse_args()

    def verify_file(fn):
        if not os.path.isfile(fn):
            raise RuntimeError(f"File not found: {fn}")
    verify_file(args.input_file)

    lz4_exe = which("lz4")
    if args.input_file.endswith(".lz4"):
        assert lz4_exe is not None
        start_unlz4 = time.time()
        input_file = args.input_file[:-4]
        command = [lz4_exe, "-f", args.input_file, input_file]
        subprocess.check_call(command, stdout=sys.stdout)
        print(f"Decompress took {time.time() - start_unlz4} seconds")
    else:
        input_file = args.input_file

    fewest_edges = 2**32
    best_tree_count = 0
    for tree_count in TRIALS:
        tree_size = SEG_SIZE / tree_count
        trees = [TreeSpec(build=Range(start=args.start_offset+(c*tree_size), end=args.start_offset+(c*tree_size)+tree_size))
                 for c in range(tree_count)]

        # Build all the trees.
        tree_args = [(input_file, tree) for tree in trees]
        print(tree_args)
        with Pool(args.jobs) as pool:
            pool.map(build_shape_wrapper, tree_args)

        # Merge them as needed.
        print("===> Merging trees")
        part_file, _ = merge_trees(input_file, trees, 0, True)

        print("===> Mapping mutations")
        result_text = map_part(input_file, part_file,
                               Range(start=args.start_offset, end=args.start_offset+SEG_SIZE),
                               False)
        print(result_text)
        for line in result_text.split("\n"):
            if line.startswith("Edges: "):
                edge_ct = int(line[7:])
                if edge_ct < fewest_edges:
                    fewest_edges = edge_ct
                    best_tree_count = tree_count
        
    print()
    print("----------------------------------------------------")
    print(f"Fewest edges: {fewest_edges}")
    print(f"Best tree count: {best_tree_count}")
    print("----------------------------------------------------")