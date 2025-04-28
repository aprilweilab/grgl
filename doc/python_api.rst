.. _python_docs:

Python API
----------

.. automodule:: pygrgl
    :members: grg_to_cyto_json, load_mutable_grg, load_immutable_grg, save_grg, save_subset, grg_from_trees, get_bfs_order, get_dfs_order, get_topo_order, dot_product, matmul, shared_frontier, INVALID_NODE, COAL_COUNT_NOT_SET, GRG_FILE_VERSION

.. automodule:: pygrgl.display
    :members: grg_to_cyto

.. autoclass:: pygrgl.Mutation
    :members: __init__, position, allele, ref_allele, time

.. autoclass:: pygrgl.GRG
    :members: is_sample, num_samples, num_individuals, ploidy, bp_range, specified_bp_range, nodes_are_ordered, mutations_are_ordered, num_nodes, num_edges, num_up_edges, num_down_edges, get_down_edges, get_up_edges, get_sample_nodes, get_root_nodes, get_node_mutation_pairs, get_mutation_node_pairs, get_mutations_for_node, get_mutation_by_id, set_mutation_by_id, node_has_mutations, add_population, get_populations, add_mutation, num_mutations, get_population_id, set_population_id, get_num_individual_coals, set_num_individual_coals

.. autoclass:: pygrgl.MutableGRG
    :members: make_node, connect, disconnect, merge

.. autoclass:: pygrgl.TraversalDirection
    :members: DOWN, UP
