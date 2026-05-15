import os
import sys
from dataclasses import dataclass, field
from typing import List, Tuple

import pygrgl as grgl


MISSING_ALLELE = "."
UNKNOWN_ALLELES = {".", "-", "N"}


@dataclass
class PolarizationStats:
    unpolarized: List[int] = field(default_factory=list)
    total_seen: int = 0
    emitted: int = 0
    already_polarized: int = 0
    swapped: int = 0
    dropped_unknown: int = 0
    inconsistent: int = 0
    after_alignment: int = 0
    no_alignment: int = 0
    indel_flip_skipped: int = 0
    missing_remapped: int = 0


@dataclass(frozen=True)
class Flip:
    mut_id: int
    node_id: int
    position: int
    ref: str
    alt: str
    time: float


def add_options(subparser):
    subparser.add_argument("grg_file", help="Input GRG file to polarize")
    subparser.add_argument("fasta_file", help="FASTA containing the ancestral sequence")
    subparser.add_argument(
        "-o",
        "--output",
        dest="output_file",
        help="Output GRG file (defaults to overwriting the input file)",
    )
    subparser.add_argument(
        "--keep-no-match",
        action="store_true",
        help="Keep mutations with no matching allele in the FASTA (instead of dropping them)",
    )
    subparser.add_argument(
        "--allow-indel-flips",
        action="store_true",
        help="Allow REF/ALT flips for length-changing alleles; default drops ambiguous indel flips",
    )
    subparser.add_argument(
        "--map-batch-size",
        type=int,
        default=4096,
        help="Number of changed mutations to collect before removing and remapping",
    )


def load_numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ValueError("grg polarize requires numpy; install it with `pip install numpy`") from error
    return np


def get_downward_indicators(grg, mutation_ids, use_missingness):
    np = load_numpy()
    values = np.zeros((len(mutation_ids), grg.num_mutations), dtype=np.float64)
    miss = None
    if use_missingness:
        miss = values
    for row, mut_id in enumerate(mutation_ids):
        values[row, mut_id] = 1.0
    if use_missingness:
        values = np.zeros_like(values)
    return grgl.matmul(grg, values, grgl.TraversalDirection.DOWN, miss=miss)


def load_remaps(grg, removals, flips, map_batch_size, stats, force=False):
    if not force and len(removals) + len(flips) < map_batch_size:
        return
    if not removals and not flips:
        return

    np = load_numpy()
    remap_mutations = []
    remap_samples = []

    if flips:
        all_samples = np.arange(grg.num_samples, dtype=np.uint32)
        mutation_ids = [flip.mut_id for flip in flips]
        carrier_matrix = get_downward_indicators(grg, mutation_ids, use_missingness=False)
        missing_matrix = get_downward_indicators(grg, mutation_ids, use_missingness=True)

        for row, flip in enumerate(flips):
            carriers = np.flatnonzero(carrier_matrix[row] > 0).astype(np.uint32)
            missing = np.flatnonzero(missing_matrix[row] > 0).astype(np.uint32)
            unavailable = np.union1d(carriers, missing)
            flipped_carriers = np.setdiff1d(all_samples, unavailable, assume_unique=True)

            removals.append((flip.mut_id, flip.node_id))

            if len(missing) > 0:
                remap_mutations.append(grgl.Mutation(flip.position, MISSING_ALLELE, flip.ref, flip.time))
                remap_samples.append(missing.tolist())
                stats.missing_remapped += 1

            remap_mutations.append(grgl.Mutation(flip.position, flip.ref, flip.alt, flip.time))
            remap_samples.append(flipped_carriers.tolist())

    list(grg.get_node_mutation_miss())
    for mut_id, node_id in removals:
        grg.remove_mutation(mut_id, node_id)
    removals.clear()
    flips.clear()

    if remap_mutations:
        grgl.map_mutations(
            grg,
            remap_mutations,
            remap_samples,
            verbose=False,
            mutation_batch_size=map_batch_size,
        )


def build_mut_lookup(grg):
    mut_lookup = {}
    for mut_id, node_id, missing_node_id in grg.get_mutation_node_miss():
        mut_lookup[int(mut_id)] = (int(node_id), int(missing_node_id))
    return mut_lookup


def polarize_mutations_helper(grg, mut_lookup, batch, stats, allow_indel_flips, map_batch_size):
    results = [False] * len(batch)
    if not batch:
        return results

    removals = []
    flips = []

    for idx, (mut_id, ancestral_allele) in enumerate(batch):
        stats.total_seen += 1

        if mut_id not in mut_lookup:
            continue

        node_id, _missing_node_id = mut_lookup[mut_id]
        if node_id == grgl.INVALID_NODE:
            continue

        mutation = grg.get_mutation_by_id(mut_id)
        position = int(mutation.position)
        ancestral_allele = ancestral_allele.upper()

        if ancestral_allele in UNKNOWN_ALLELES:
            stats.dropped_unknown += 1
            stats.unpolarized.append(position)
            removals.append((mut_id, node_id))
            load_remaps(grg, removals, flips, map_batch_size, stats)
            continue

        ref = mutation.ref_allele
        alt = mutation.allele

        if ancestral_allele == ref.upper():
            stats.already_polarized += 1
            stats.emitted += 1
            results[idx] = True
            continue

        if ancestral_allele == alt.upper():
            if len(ref) != len(alt) and not allow_indel_flips:
                stats.indel_flip_skipped += 1
                stats.unpolarized.append(position)
                removals.append((mut_id, node_id))
                load_remaps(grg, removals, flips, map_batch_size, stats)
                continue

            stats.swapped += 1
            stats.emitted += 1
            results[idx] = True
            flips.append(
                Flip(
                    mut_id=mut_id,
                    node_id=node_id,
                    position=position,
                    ref=ref,
                    alt=alt,
                    time=mutation.time,
                )
            )
            load_remaps(grg, removals, flips, map_batch_size, stats)
            continue

        stats.inconsistent += 1
        stats.unpolarized.append(position)
        removals.append((mut_id, node_id))
        load_remaps(grg, removals, flips, map_batch_size, stats)

    load_remaps(grg, removals, flips, map_batch_size, stats, force=True)
    return results


def polarize_mutations(grg, batch, stats, allow_indel_flips=False, map_batch_size=4096):
    mut_lookup = build_mut_lookup(grg)
    results = polarize_mutations_helper(grg, mut_lookup, batch, stats, allow_indel_flips, map_batch_size)
    grg.sort_mutations()
    return results


def load_fasta(path: str) -> Tuple[object, str]:
    try:
        import pyfaidx
    except ImportError as error:
        raise ValueError("grg polarize requires pyfaidx; install it with `pip install pyfaidx`") from error

    fasta = pyfaidx.Fasta(path, as_raw=False, sequence_always_upper=True)
    contigs = list(fasta.keys())
    if len(contigs) != 1:
        raise ValueError(
            "FASTA must contain exactly one contig for grg polarize. "
            f"Found: {', '.join(contigs)}"
        )
    return fasta, contigs[0]


def fetch_ancestral(fasta, contig, position, length):
    if length <= 0:
        return None
    if position <= 0:
        return None
    if position + length - 1 > len(fasta[contig]):
        return None
    return str(fasta.get_seq(contig, position, position + length - 1)).upper()


def equals_ignore_case(left, right):
    return left.upper() == right.upper()


def polarize_grg_from_fasta(
    grg,
    fasta_file,
    drop_if_no_match=True,
    allow_indel_flips=False,
    map_batch_size=4096,
):
    stats = PolarizationStats()
    fasta, contig = load_fasta(fasta_file)
    mut_lookup = build_mut_lookup(grg)
    batch = []
    total_mutations = grg.num_mutations

    for mut_id in range(total_mutations):
        mutation = grg.get_mutation_by_id(mut_id)
        if mutation.allele == MISSING_ALLELE:
            continue

        ref = mutation.ref_allele
        alt = mutation.allele
        max_allele_len = max(len(ref), len(alt))
        ancestral = fetch_ancestral(fasta, contig, int(mutation.position), max_allele_len)

        if ancestral is None:
            stats.after_alignment += 1
            if drop_if_no_match:
                batch.append((mut_id, "N"))
            continue

        if ancestral[0] in UNKNOWN_ALLELES:
            stats.no_alignment += 1
            if drop_if_no_match:
                batch.append((mut_id, "N"))
            continue

        ref_matches = bool(ref) and equals_ignore_case(ref, ancestral[: len(ref)])
        alt_matches = bool(alt) and equals_ignore_case(alt, ancestral[: len(alt)])

        if ref_matches and alt_matches:
            batch.append((mut_id, ref if len(ref) >= len(alt) else alt))
        elif ref_matches:
            batch.append((mut_id, ref))
        elif alt_matches:
            batch.append((mut_id, alt))
        elif drop_if_no_match:
            batch.append((mut_id, "N"))

        if len(batch) >= map_batch_size:
            polarize_mutations_helper(grg, mut_lookup, batch, stats, allow_indel_flips, map_batch_size)
            batch.clear()

    if batch:
        polarize_mutations_helper(grg, mut_lookup, batch, stats, allow_indel_flips, map_batch_size)

    grg.sort_mutations()
    return stats


def polarize_command(arguments):
    if not os.path.isfile(arguments.grg_file):
        print(f"Input GRG file does not exist: {arguments.grg_file}", file=sys.stderr)
        sys.exit(2)
    if not os.path.isfile(arguments.fasta_file):
        print(f"FASTA file does not exist: {arguments.fasta_file}", file=sys.stderr)
        sys.exit(2)

    output_path = arguments.output_file if arguments.output_file else arguments.grg_file

    grg = grgl.load_mutable_grg(arguments.grg_file, load_up_edges=True)
    if grg is None:
        print(f"Failed to load GRG: {arguments.grg_file}", file=sys.stderr)
        sys.exit(2)

    try:
        stats = polarize_grg_from_fasta(
            grg,
            arguments.fasta_file,
            drop_if_no_match=not arguments.keep_no_match,
            allow_indel_flips=arguments.allow_indel_flips,
            map_batch_size=arguments.map_batch_size,
        )
    except ValueError as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)

    grgl.save_grg(grg, output_path)

    print("Polarization complete")
    print(f"  Total seen:           {stats.total_seen}")
    print(f"  Emitted:              {stats.emitted}")
    print(f"  Already polarized:    {stats.already_polarized}")
    print(f"  Swapped:              {stats.swapped}")
    print(f"  Dropped unknown:      {stats.dropped_unknown}")
    print(f"  Inconsistent:         {stats.inconsistent}")
    print(f"  After alignment end:  {stats.after_alignment}")
    print(f"  Alignment mismatch:   {stats.no_alignment}")
    print(f"  Indel flips skipped:  {stats.indel_flip_skipped}")
    print(f"  Missingness remapped: {stats.missing_remapped}")
