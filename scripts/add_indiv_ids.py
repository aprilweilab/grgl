# Add Individual IDs to an existing GRG that does not have them.
# IDs are provided via a text file, with a single line per ID.
import pygrgl
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: add_indiv_ids.py <GRG filename> <ID text filename>", file=sys.stderr)
        exit(1)

    identifiers = []
    with open(sys.argv[2]) as f:
        for line in f:
            identifiers.append(line.strip())
    
    GRG = sys.argv[1]
    assert GRG.endswith(".grg"), "The GRG file must end with .grg"
    prefix = GRG[:-4]

    grg = pygrgl.load_immutable_grg(GRG)
    assert not grg.has_individual_ids, f"Provided GRG already has individual identifiers"
    assert grg.num_individuals == len(identifiers), f"GRG has {grg.num_individuals} individuals, but the ID text file has {len(identifiers)} lines"

    for ident in identifiers:
        grg.add_individual_id(ident)
    pygrgl.save_grg(grg, f"{prefix}.WITHIDS.grg")
    print(f"Saved result to {prefix}.WITHIDS.grg")
