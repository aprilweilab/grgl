import subprocess
import os
import pygrgl
import numpy
from typing import Optional

JOBS = 4
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(THIS_DIR, "input")


def construct_grg(
    input_file: str,
    output_file: Optional[str] = None,
    test_input: bool = True,
    parts: int = 10,
) -> str:
    cmd = [
        "grg",
        "construct",
        "--force",  # Ignore that we don't have tabix-indexed files
        "-p",
        str(parts),
        "-j",
        str(JOBS),
        os.path.join(INPUT_DIR, input_file) if test_input else input_file,
    ]
    if output_file is not None:
        cmd.extend(["-o", output_file])
    else:
        output_file = os.path.basename(input_file) + ".final.grg"
    subprocess.check_call(cmd)
    return output_file


def allele_frequencies(grg: pygrgl.GRG) -> numpy.typing.NDArray:
    with numpy.errstate(divide="raise"):
        kwargs = {}
        miss = 0
        if grg.has_missing_data:
            miss = numpy.zeros((1, grg.num_mutations), dtype=numpy.int32)
            kwargs["miss"] = miss
        else:
            miss = None
        counts = pygrgl.matmul(
            grg,
            numpy.ones((1, grg.num_samples), dtype=numpy.int32),
            pygrgl.TraversalDirection.UP,
            **kwargs,
        )[0]
        if miss is None:
            miss = 0
        else:
            miss = miss[0]
        denominator = grg.num_samples - miss
        return numpy.divide(
            counts,
            denominator,
            out=numpy.zeros(counts.shape, dtype=numpy.float64),
            where=(denominator != 0),
        )


# It is important that this only uses down edges, so that we can test against
# matmul methods which (should) only require down edges.
def grg2X(grg: pygrgl.GRG, diploid: bool = False):
    samples = grg.num_individuals if diploid else grg.num_samples
    result = numpy.zeros((samples, grg.num_mutations))
    samps_below = [list() for _ in range(grg.num_nodes)]
    for node_id in range(grg.num_nodes):
        sb = []
        if grg.is_sample(node_id):
            sb.append(node_id)
        for child_id in grg.get_down_edges(node_id):
            sb.extend(samps_below[child_id])
        samps_below[node_id] = sb

        muts = grg.get_mutations_for_node(node_id)
        if muts:
            for sample_id in sb:
                indiv = sample_id // grg.ploidy
                for mut_id in muts:
                    if diploid:
                        result[indiv][mut_id] += 1
                    else:
                        result[sample_id][mut_id] = 1
    # Handle missingness afterwards for simplicity (harder to make mistakes this way).
    if grg.has_missing_data:
        freqs = allele_frequencies(grg)
        for mut_id, mut_node, miss_node in grg.get_mutation_node_miss():
            if miss_node != pygrgl.INVALID_NODE:
                for sample_id in samps_below[miss_node]:
                    if diploid:
                        indiv = sample_id // grg.ploidy
                        result[indiv][mut_id] += freqs[mut_id]
                    else:
                        assert (
                            result[sample_id][mut_id] == 0
                        ), f"{result[sample_id][mut_id]}"
                        result[sample_id][mut_id] = freqs[mut_id]
    return result
