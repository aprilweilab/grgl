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
    missing_node_id: int
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
        help="Number of flipped mutations to process per carrier-read and remap batch",
    )


def load_numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ValueError("grg polarize requires numpy; install it with `pip install numpy`") from error
    return np


def get_downward_indicators(grg, mutation_ids):
    np = load_numpy()
    values = np.zeros((len(mutation_ids), grg.num_mutations), dtype=np.int32)
    for row, mut_id in enumerate(mutation_ids):
        values[row, mut_id] = 1
    return grgl.matmul(grg, values, grgl.TraversalDirection.DOWN)


def get_missingness_indicators(grg, missing_node_ids):
    np = load_numpy()
    if not any(missing_node_id != grgl.INVALID_NODE for missing_node_id in missing_node_ids):
        return np.zeros((len(missing_node_ids), grg.num_samples), dtype=np.int32)

    values = np.zeros((len(missing_node_ids), grg.num_mutations), dtype=np.int32)
    init = np.zeros((len(missing_node_ids), grg.num_nodes), dtype=np.int32)
    for row, missing_node_id in enumerate(missing_node_ids):
        if missing_node_id != grgl.INVALID_NODE:
            init[row, missing_node_id] = 1
    return grgl.matmul(grg, values, grgl.TraversalDirection.DOWN, init=init)


def build_flip_remaps(grg, flips, map_batch_size, stats):
    if not flips:
        return [], []

    np = load_numpy()
    remap_mutations = []
    remap_samples = []

    for start in range(0, len(flips), map_batch_size):
        flip_batch = flips[start : start + map_batch_size]
        all_samples = np.arange(grg.num_samples, dtype=np.uint32)
        mutation_ids = [flip.mut_id for flip in flip_batch]
        missing_node_ids = [flip.missing_node_id for flip in flip_batch]
        carrier_matrix = get_downward_indicators(grg, mutation_ids)
        missing_matrix = get_missingness_indicators(grg, missing_node_ids)

        for row, flip in enumerate(flip_batch):
            carriers = np.flatnonzero(carrier_matrix[row] > 0).astype(np.uint32)
            missing = np.flatnonzero(missing_matrix[row] > 0).astype(np.uint32)
            unavailable = np.union1d(carriers, missing)
            flipped_carriers = np.setdiff1d(all_samples, unavailable, assume_unique=True)

            if len(missing) > 0:
                remap_mutations.append(grgl.Mutation(flip.position, MISSING_ALLELE, flip.ref, flip.time))
                remap_samples.append(missing.tolist())
                stats.missing_remapped += 1

            remap_mutations.append(grgl.Mutation(flip.position, flip.ref, flip.alt, flip.time))
            remap_samples.append(flipped_carriers.tolist())

    return remap_mutations, remap_samples


def apply_remaps(grg, removals, remap_mutations, remap_samples, map_batch_size):
    if not removals and not remap_mutations:
        return

    list(grg.get_node_mutation_miss())
    for mut_id, node_id in removals:
        grg.remove_mutation(mut_id, node_id)

    if remap_mutations:
        grgl.map_mutations(
            grg,
            remap_mutations,
            remap_samples,
            verbose=False,
            mutation_batch_size=map_batch_size,
        )


def build_site_swap_remaps(grg, site_entries, swap_index, ancestral_allele, stats):
    np = load_numpy()
    all_samples = np.arange(grg.num_samples, dtype=np.uint32)
    mutation_ids = [entry[0] for entry in site_entries]
    missing_node_ids = [entry[2] for entry in site_entries]
    carrier_matrix = get_downward_indicators(grg, mutation_ids)
    missing_matrix = get_missingness_indicators(grg, missing_node_ids)

    old_ref = site_entries[0][3].ref_allele
    position = int(site_entries[0][3].position)
    time = site_entries[swap_index][3].time
    missing = np.flatnonzero(missing_matrix[swap_index] > 0).astype(np.uint32)
    unavailable = missing

    remap_mutations = []
    remap_samples = []
    if len(missing) > 0:
        remap_mutations.append(grgl.Mutation(position, MISSING_ALLELE, ancestral_allele, time))
        remap_samples.append(missing.tolist())
        stats.missing_remapped += 1

    for row, (_mut_id, _node_id, _missing_node_id, mutation) in enumerate(site_entries):
        carriers = np.flatnonzero(carrier_matrix[row] > 0).astype(np.uint32)
        unavailable = np.union1d(unavailable, carriers)
        if row != swap_index:
            remap_mutations.append(grgl.Mutation(position, mutation.allele, ancestral_allele, mutation.time))
            remap_samples.append(carriers.tolist())

    old_ref_carriers = np.setdiff1d(all_samples, unavailable, assume_unique=True)
    remap_mutations.append(grgl.Mutation(position, old_ref, ancestral_allele, time))
    remap_samples.append(old_ref_carriers.tolist())
    return remap_mutations, remap_samples


def polarize_site_from_fasta(grg, site_entries, ancestral_allele, stats, drop_if_no_match, map_batch_size):
    if not site_entries:
        return

    position = int(site_entries[0][3].position)
    old_ref = site_entries[0][3].ref_allele
    stats.total_seen += len(site_entries)

    if ancestral_allele is None:
        stats.after_alignment += 1
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals = [(mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries]
            apply_remaps(grg, removals, [], [], map_batch_size)
        return

    ancestral_allele = ancestral_allele.upper()
    if ancestral_allele in UNKNOWN_ALLELES:
        stats.no_alignment += 1
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals = [(mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries]
            apply_remaps(grg, removals, [], [], map_batch_size)
        return

    if ancestral_allele == old_ref.upper():
        stats.already_polarized += len(site_entries)
        stats.emitted += len(site_entries)
        return

    swap_index = -1
    for idx, (_mut_id, _node_id, _missing_node_id, mutation) in enumerate(site_entries):
        if mutation.allele.upper() == ancestral_allele:
            swap_index = idx
            break

    if swap_index < 0:
        stats.inconsistent += len(site_entries)
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals = [(mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries]
            apply_remaps(grg, removals, [], [], map_batch_size)
        return

    stats.swapped += 1
    stats.emitted += len(site_entries)
    removals = [(mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries]
    remap_mutations, remap_samples = build_site_swap_remaps(grg, site_entries, swap_index, ancestral_allele, stats)
    apply_remaps(grg, removals, remap_mutations, remap_samples, map_batch_size)


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

        node_id, missing_node_id = mut_lookup[mut_id]
        if node_id == grgl.INVALID_NODE:
            continue

        mutation = grg.get_mutation_by_id(mut_id)
        position = int(mutation.position)
        ancestral_allele = ancestral_allele.upper()

        if ancestral_allele in UNKNOWN_ALLELES:
            stats.dropped_unknown += 1
            stats.unpolarized.append(position)
            removals.append((mut_id, node_id))
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
                continue

            stats.swapped += 1
            stats.emitted += 1
            results[idx] = True
            flips.append(
                Flip(
                    mut_id=mut_id,
                    node_id=node_id,
                    missing_node_id=missing_node_id,
                    position=position,
                    ref=ref,
                    alt=alt,
                    time=mutation.time,
                )
            )
            continue

        stats.inconsistent += 1
        stats.unpolarized.append(position)
        removals.append((mut_id, node_id))

    remap_mutations, remap_samples = build_flip_remaps(grg, flips, map_batch_size, stats)
    apply_remaps(
        grg,
        removals + [(flip.mut_id, flip.node_id) for flip in flips],
        remap_mutations,
        remap_samples,
        map_batch_size,
    )
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
    if map_batch_size <= 0:
        raise ValueError("--map-batch-size must be greater than zero")

    stats = PolarizationStats()
    fasta, contig = load_fasta(fasta_file)
    mut_lookup = build_mut_lookup(grg)
    total_mutations = grg.num_mutations
    site_entries = []
    site_key = None

    def flush_site():
        if not site_entries:
            return
        position = int(site_entries[0][3].position)
        ancestral = fetch_ancestral(fasta, contig, position, 1)
        polarize_site_from_fasta(grg, site_entries, ancestral, stats, drop_if_no_match, map_batch_size)

    for mut_id in range(total_mutations):
        mutation = grg.get_mutation_by_id(mut_id)
        if mutation.allele == MISSING_ALLELE:
            continue

        if mut_id not in mut_lookup:
            continue

        node_id, missing_node_id = mut_lookup[mut_id]
        if node_id == grgl.INVALID_NODE:
            continue

        key = (int(mutation.position), mutation.ref_allele)
        if site_key is not None and key != site_key:
            flush_site()
            site_entries.clear()
        site_key = key
        site_entries.append((mut_id, node_id, missing_node_id, mutation))

    flush_site()

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
