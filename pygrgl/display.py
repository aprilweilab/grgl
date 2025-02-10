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
try:
    from ipycytoscape import CytoscapeWidget
except ImportError:
    print(
        "Unable to import ipycytoscape; run 'pip install ipycytoscape' to enable graph display."
    )
    CytoscapeWidget = None

from .grg import GRG, grg_to_cyto_json

DAG_STYLE = [
    {
        "selector": "node",
        "style": {
            "font-family": "arial",
            "font-size": "10px",
            "label": "data(label)",
            "background-color": "blue",
        },
    },
    {
        "selector": 'node[is_sample = "True"]',
        "style": {
            "font-family": "arial",
            "font-size": "10px",
            "label": "data(label)",
            "background-color": "green",
        },
    },
    {
        "selector": "edge",
        "style": {
            "width": 4,
            "line-color": "#9dbaea",
            "target-arrow-shape": "triangle",
            "target-arrow-color": "#9dbaea",
            "curve-style": "bezier",
        },
    },
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
        grg_to_cyto_json(grg, start_from=start_from, show_mutations=show_mutations),
        directed=True,
    )
    cyto.set_style(DAG_STYLE)
    cyto.set_layout(name="dagre")
    return cyto
