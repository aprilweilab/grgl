try:
    from ipycytoscape import CytoscapeWidget
except ImportError:
    print("Unable to import ipycytoscape; run 'pip install ipycytoscape' to enable graph display.")
    CytoscapeWidget = None

from .grg import GRG, grg_to_cyto_json

DAG_STYLE = [
    {'selector': 'node',
     'style': {
        'font-family': 'arial',
        'font-size': '10px',
        'label': 'data(label)',
        'background-color': 'blue'}},
    {'selector': 'node[is_sample = "True"]',
     'style': {
        'font-family': 'arial',
        'font-size': '10px',
        'label': 'data(label)',
        'background-color': 'green'}},
    {'selector': 'edge',
     'style': {
        'width': 4,
        'line-color': '#9dbaea',
        'target-arrow-shape': 'triangle',
        'target-arrow-color': '#9dbaea',
        'curve-style': 'bezier'
    }},
    ]

def grg_to_cyto(grg: GRG, start_from=[], show_mutations=True) -> CytoscapeWidget:
    """
    Return a CytoscapeWidget that can be displayed in Jupyter.
    """
    if CytoscapeWidget is None:
        raise RuntimeError("Please install ipycytoscape")
    cyto = CytoscapeWidget()
    # Note: directed=true requires the stylesheet to have arrow styles on edges.
    cyto.graph.add_graph_from_json(
        grg_to_cyto_json(grg, start_from=start_from, show_mutations=show_mutations), directed=True)
    cyto.set_style(DAG_STYLE)
    cyto.set_layout(name="dagre")
    return cyto
