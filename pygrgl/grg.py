from _grgl import *
from typing import Dict, Any


def grg_to_cyto_json(grg: GRG, start_from=[], show_mutations=True) -> Dict[str, Any]:
    """
    Create a dictionary that can be used with ipycytoscape to display the graph.
    Right now this only emits down edges (source is older than target).
    """
    def label(node_id: int) -> str:
        label = f"id={node_id}"
        if grg.is_sample(node_id):
            label += f", S"
        mutations = grg.get_mutations_for_node(node_id)
        if show_mutations and mutations:
            label += f", mutations({set(mutations)})"
        return label

    nodeset = start_from if start_from else list(range(0, grg.num_nodes))
    def make_node(n):
        return {"data": {"id": f"n{n}",
                "label": label(n),
                "is_sample": str(grg.is_sample(n))}}

    nodes = {n: make_node(n) for n in nodeset}
    edges = []
    edge_counter = 0
    while nodeset:
        n = nodeset.pop(-1)
        for child_id in grg.get_down_edges(n):
            edges.append({"data": {"id": f"e{edge_counter}", "source": f"n{n}",
                                   "target": f"n{child_id}"}})
            if child_id not in nodes:
                nodes[child_id] = make_node(child_id)
                nodeset.append(child_id)
            edge_counter += 1
    return {"nodes": list(nodes.values()), "edges": edges}
