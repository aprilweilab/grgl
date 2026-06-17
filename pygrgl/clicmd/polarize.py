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
    inconsistent: int = 0
    after_alignment: int = 0
    no_alignment: int = 0
    non_snv_skipped: int = 0
    missing_remapped: int = 0


@dataclass(frozen=True)
class SiteSwap:
    entries: list
    swap_index: int
    ancestral_allele: str


def add_options(subparser):
    subparser.add_argument("grg_file", help="Input GRG file to polarize")
    subparser.add_argument("fasta_file", help="FASTA containing the ancestral sequence")
    subparser.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=True,
        help="Output GRG file",
    )
    subparser.add_argument(
        "--keep-no-match",
        action="store_true",
        help="Keep mutations with no matching allele in the FASTA (instead of dropping them)",
    )
    subparser.add_argument(
        "--map-batch-size",
        type=int,
        default=4096,
        help="Number of flipped mutations to process per carrier-read and remap batch",
    )
    subparser.add_argument(
        "--thread-count",
        type=int,
        default=1,
        help="Number of worker threads to use for candidate processing during remaps",
    )


def load_numpy():
    try:
        import numpy as np
    except ImportError as error:
        raise ValueError("grg polarize requires numpy; install it with `pip install numpy`") from error
    return np


def get_descendant_samples(grg, node_id, np):
    if node_id == grgl.INVALID_NODE:
        return np.array([], dtype=np.uint32)
    descendants = grgl.get_bfs_order(grg, grgl.TraversalDirection.DOWN, [node_id])
    return np.fromiter(
        (node for node in descendants if node < grg.num_samples),
        dtype=np.uint32,
    )


# apply all remaps in a batch (can be from multiple sites)
def apply_remaps(grg, removals, remap_mutations, remap_samples, map_batch_size, thread_count):
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
            thread_count=thread_count,
        )


# compute mutation removals and remaps for one site
def build_site_swap_remaps(grg, site_swaps, map_batch_size, stats):
    if not site_swaps:
        return [], [], []

    np = load_numpy()
    removals = []
    remap_mutations = []
    remap_samples = []

    site_state = []
    flat_entries = []
    for site_index, site_swap in enumerate(site_swaps):
        entries = site_swap.entries
        swap_mutation = entries[site_swap.swap_index][3]
        has_missing = any(entry[2] != grgl.INVALID_NODE for entry in entries)
        site_state.append(
            {
                "ancestral_allele": site_swap.ancestral_allele,
                "old_ref": entries[0][3].ref_allele,
                "position": int(entries[0][3].position),
                "time": swap_mutation.time,
                "unavailable": np.zeros(grg.num_samples, dtype=bool),
                "full_remap": has_missing,
            }
        )
        for row_index, entry in enumerate(entries):
            flat_entries.append((site_index, row_index, site_swap.swap_index, entry))

    for start in range(0, len(flat_entries), map_batch_size):
        entry_batch = flat_entries[start : start + map_batch_size]

        for site_index, row_index, swap_index, entry in entry_batch:
            _mut_id, _node_id, _missing_node_id, mutation = entry
            state = site_state[site_index]
            carriers = get_descendant_samples(grg, _node_id, np)
            state["unavailable"][carriers] = True

            if row_index == swap_index:
                removals.append((_mut_id, _node_id))
                missing = get_descendant_samples(grg, _missing_node_id, np)
                state["unavailable"][missing] = True
                if len(missing) > 0:
                    remap_mutations.append(
                        grgl.Mutation(
                            state["position"],
                            MISSING_ALLELE,
                            state["ancestral_allele"],
                            state["time"],
                        )
                    )
                    remap_samples.append(missing.tolist())
                    stats.missing_remapped += 1
            else:
                updated_mutation = grgl.Mutation(
                    state["position"],
                    mutation.allele,
                    state["ancestral_allele"],
                    mutation.time,
                )
                if state["full_remap"]:
                    removals.append((_mut_id, _node_id))
                    remap_mutations.append(updated_mutation)
                    remap_samples.append(carriers.tolist())
                else:
                    grg.set_mutation_by_id(_mut_id, updated_mutation)

    for state in site_state:
        old_ref_carriers = np.flatnonzero(~state["unavailable"]).astype(np.uint32)
        remap_mutations.append(
            grgl.Mutation(
                state["position"],
                state["old_ref"],
                state["ancestral_allele"],
                state["time"],
            )
        )
        remap_samples.append(old_ref_carriers.tolist())

    return removals, remap_mutations, remap_samples


# determine whether a site is already polarized or needs a swap
def classify_site(
    site_entries,
    ancestral_sequence,
    stats,
    drop_if_no_match,
    removals,
    site_swaps,
):
    if not site_entries:
        return

    position = int(site_entries[0][3].position)
    old_ref = site_entries[0][3].ref_allele
    stats.total_seen += len(site_entries)

    if len(old_ref) != 1 or any(
        len(mutation.allele) != 1 for _mut_id, _node_id, _missing_node_id, mutation in site_entries
    ):
        stats.non_snv_skipped += len(site_entries)
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals.extend((mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries)
        return

    if ancestral_sequence is None:
        stats.after_alignment += 1
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals.extend((mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries)
        return

    ancestral_sequence = ancestral_sequence.upper()
    if ancestral_sequence[0] in UNKNOWN_ALLELES:
        stats.no_alignment += 1
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals.extend((mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries)
        return

    matches = []
    if allele_matches_ancestral(old_ref, ancestral_sequence):
        matches.append((-1, old_ref))

    for idx, (_mut_id, _node_id, _missing_node_id, mutation) in enumerate(site_entries):
        if allele_matches_ancestral(mutation.allele, ancestral_sequence):
            matches.append((idx, mutation.allele))

    if not matches:
        stats.inconsistent += len(site_entries)
        stats.unpolarized.append(position)
        if drop_if_no_match:
            removals.extend((mut_id, node_id) for mut_id, node_id, _missing_node_id, _mutation in site_entries)
        return

    swap_index, ancestral_allele = max(matches, key=lambda item: len(item[1]))

    if swap_index < 0:
        stats.already_polarized += len(site_entries)
        stats.emitted += len(site_entries)
        return

    stats.swapped += 1
    stats.emitted += len(site_entries)
    site_swaps.append(SiteSwap(list(site_entries), swap_index, ancestral_allele.upper()))


def build_mut_lookup(grg):
    mut_lookup = [None] * grg.num_mutations
    for mut_id, node_id, missing_node_id in grg.get_mutation_node_miss():
        mut_lookup[int(mut_id)] = (int(node_id), int(missing_node_id))
    return mut_lookup


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


def fetch_ancestral(sequence, position, length):
    if length <= 0:
        return None
    if position <= 0:
        return None
    if position + length - 1 > len(sequence):
        return None
    return sequence[position - 1 : position - 1 + length]


def equals_ignore_case(left, right):
    return left.upper() == right.upper()


def allele_matches_ancestral(allele, ancestral_sequence):
    return bool(allele) and len(allele) <= len(ancestral_sequence) and equals_ignore_case(
        allele,
        ancestral_sequence[: len(allele)],
    )


def polarize_grg_from_fasta(
    grg,
    fasta_file,
    drop_if_no_match=True,
    map_batch_size=4096,
    thread_count=1,
):
    if map_batch_size <= 0:
        raise ValueError("--map-batch-size must be greater than zero")
    if thread_count <= 0:
        raise ValueError("--thread-count must be greater than zero")

    stats = PolarizationStats()
    fasta, contig = load_fasta(fasta_file)
    ancestral_sequence = str(fasta[contig][:]).upper()
    mut_lookup = build_mut_lookup(grg)
    total_mutations = grg.num_mutations
    site_entries = []
    site_key = None
    pending_removals = []
    pending_site_swaps = []
    pending_swap_entry_count = 0
    pending_map_removals = []
    pending_remap_mutations = []
    pending_remap_samples = []

    def materialize_site_swaps():
        nonlocal pending_swap_entry_count
        if not pending_site_swaps:
            return
        site_removals, remap_mutations, remap_samples = build_site_swap_remaps(
            grg, pending_site_swaps, map_batch_size, stats
        )
        pending_map_removals.extend(site_removals)
        pending_remap_mutations.extend(remap_mutations)
        pending_remap_samples.extend(remap_samples)
        pending_site_swaps.clear()
        pending_swap_entry_count = 0

    def flush_remaps():
        if not pending_map_removals and not pending_remap_mutations:
            return
        apply_remaps(
            grg,
            pending_map_removals,
            pending_remap_mutations,
            pending_remap_samples,
            map_batch_size,
            thread_count,
        )
        pending_map_removals.clear()
        pending_remap_mutations.clear()
        pending_remap_samples.clear()

    def flush_pending(force=False):
        nonlocal pending_swap_entry_count
        if force:
            materialize_site_swaps()
            if pending_remap_mutations:
                pending_map_removals[:0] = pending_removals
                pending_removals.clear()
                flush_remaps()
            elif pending_removals:
                apply_remaps(grg, pending_removals, [], [], map_batch_size, thread_count)
                pending_removals.clear()
            return

        if pending_swap_entry_count >= map_batch_size:
            materialize_site_swaps()
        if len(pending_remap_mutations) >= map_batch_size:
            pending_map_removals[:0] = pending_removals
            pending_removals.clear()
            flush_remaps()
            return
        if pending_removals and not pending_site_swaps and len(pending_removals) >= map_batch_size:
            apply_remaps(grg, pending_removals, [], [], map_batch_size, thread_count)
            pending_removals.clear()

    def flush_site():
        nonlocal pending_swap_entry_count
        if not site_entries:
            return
        previous_swap_count = len(pending_site_swaps)
        position = int(site_entries[0][3].position)
        ancestral = fetch_ancestral(ancestral_sequence, position, 1)
        classify_site(
            site_entries,
            ancestral,
            stats,
            drop_if_no_match,
            pending_removals,
            pending_site_swaps,
        )
        if len(pending_site_swaps) > previous_swap_count:
            pending_swap_entry_count += len(site_entries)
        flush_pending()

    for mut_id in range(total_mutations):
        mutation = grg.get_mutation_by_id(mut_id)
        if mutation.allele == MISSING_ALLELE:
            continue

        lookup_entry = mut_lookup[mut_id]
        if lookup_entry is None:
            continue

        node_id, missing_node_id = lookup_entry
        if node_id == grgl.INVALID_NODE:
            continue

        key = (int(mutation.position), mutation.ref_allele)
        if site_key is not None and key != site_key:
            flush_site()
            site_entries.clear()
        site_key = key
        site_entries.append((mut_id, node_id, missing_node_id, mutation))

    flush_site()
    flush_pending(force=True)

    grg.sort_mutations()
    return stats


def polarize_command(arguments):
    if not os.path.isfile(arguments.grg_file):
        print(f"Input GRG file does not exist: {arguments.grg_file}", file=sys.stderr)
        sys.exit(2)
    if not os.path.isfile(arguments.fasta_file):
        print(f"FASTA file does not exist: {arguments.fasta_file}", file=sys.stderr)
        sys.exit(2)

    grg = grgl.load_mutable_grg(arguments.grg_file, load_up_edges=True)
    if grg is None:
        print(f"Failed to load GRG: {arguments.grg_file}", file=sys.stderr)
        sys.exit(2)

    try:
        stats = polarize_grg_from_fasta(
            grg,
            arguments.fasta_file,
            drop_if_no_match=not arguments.keep_no_match,
            map_batch_size=arguments.map_batch_size,
            thread_count=arguments.thread_count,
        )
    except ValueError as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)

    grgl.save_grg(grg, arguments.output_file)

    print("Polarization complete")
    print(f"  Total seen:           {stats.total_seen}")
    print(f"  Emitted:              {stats.emitted}")
    print(f"  Already polarized:    {stats.already_polarized}")
    print(f"  Swapped:              {stats.swapped}")
    print(f"  Inconsistent:         {stats.inconsistent}")
    print(f"  After alignment end:  {stats.after_alignment}")
    print(f"  Alignment mismatch:   {stats.no_alignment}")
    print(f"  Non-SNVs skipped:     {stats.non_snv_skipped}")
    print(f"  Missingness remapped: {stats.missing_remapped}")
