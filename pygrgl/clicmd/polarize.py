import os
import sys
import pygrgl as grgl


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
        "--one-based-indexing",
        action="store_true",
        help="Treat mutation positions as 1-based when indexing into the FASTA (default is 0-based)",
    )


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

    stats = grgl.polarize_from_fasta(
        grg,
        arguments.fasta_file,
        drop_if_no_match=not arguments.keep_no_match,
        positions_are_one_based=arguments.one_based_indexing,
    )

    grgl.save_grg(grg, output_path)

    print("Polarization complete")
    print(f"  Total seen:        {stats.total_seen}")
    print(f"  Emitted:           {stats.emitted}")
    print(f"  Already polarized: {stats.already_polarized}")
    print(f"  Swapped:           {stats.swapped}")
    print(f"  Dropped unknown:   {stats.dropped_unknown}")
    print(f"  Inconsistent:      {stats.inconsistent}")
