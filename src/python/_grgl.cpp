/* Genotype Representation Graph Library (GRGL)
 * Copyright (C) 2024 April Wei
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "grg_helpers.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/serialize.h"
#include "grgl/ts2grg.h"
#include "grgl/visitor.h"

namespace py = pybind11;

#include <iostream>

py::array_t<double> dotProduct(grgl::GRGPtr& grg, py::array_t<double> input, grgl::TraversalDirection direction) {
    py::buffer_info buffer = input.request();

    const size_t outSize =
        (direction == grgl::TraversalDirection::DIRECTION_DOWN) ? grg->numSamples() : grg->numMutations();
    py::array_t<double> result(outSize);
    py::buffer_info resultBuf = result.request();
    memset(resultBuf.ptr, 0, outSize * sizeof(double));
    grg->dotProduct((const double*)buffer.ptr, (size_t)buffer.size, direction, (double*)resultBuf.ptr, (size_t)outSize);
    return std::move(result);
}

class NodeNumberingIterator : public grgl::GRGVisitor {
public:
    NodeNumberingIterator(grgl::DfsPass pass)
        : m_dfsPass(pass) {}

    bool visit(const grgl::GRGPtr& grg,
               grgl::NodeID nodeId,
               grgl::TraversalDirection direction,
               grgl::DfsPass dfsPass = grgl::DFS_PASS_NONE) override {
        if (dfsPass != m_dfsPass) {
            return true;
        }
        m_nodeIds.push_back(nodeId);
        return true;
    }

    std::vector<grgl::NodeID> m_nodeIds;
    grgl::DfsPass m_dfsPass;
};

std::vector<grgl::NodeID> getBfsOrder(const grgl::GRGPtr& grg,
                                      grgl::TraversalDirection direction,
                                      const grgl::NodeIDList& seedList,
                                      ssize_t maxQueueWidth = -1) {
    NodeNumberingIterator iterator(grgl::DFS_PASS_NONE);
    grg->visitBfs(iterator, direction, seedList, maxQueueWidth);
    return std::move(iterator.m_nodeIds);
}

std::vector<grgl::NodeID> getDfsOrder(const grgl::GRGPtr& grg,
                                      grgl::TraversalDirection direction,
                                      const grgl::NodeIDList& seedList,
                                      bool forwardOnly = false) {
    NodeNumberingIterator iterator(forwardOnly ? grgl::DFS_PASS_THERE : grgl::DFS_PASS_BACK_AGAIN);
    grg->visitDfs(iterator, direction, seedList, forwardOnly);
    return std::move(iterator.m_nodeIds);
}

std::vector<grgl::NodeID>
getTopoOrder(const grgl::GRGPtr& grg, grgl::TraversalDirection direction, const grgl::NodeIDList& seedList) {
    NodeNumberingIterator iterator(grgl::DFS_PASS_NONE);
    grg->visitTopo(iterator, direction, seedList);
    return std::move(iterator.m_nodeIds);
}

size_t hashMutation(const grgl::Mutation* self) { return std::hash<grgl::Mutation>()(*self); }

std::vector<std::pair<grgl::NodeID, grgl::MutationId>> getNodeMutationPairs(const grgl::GRGPtr& grg) {
    std::vector<std::pair<grgl::NodeID, grgl::MutationId>> result;
    for (const auto& nodeAndMutId : grg->getNodeMutationPairs()) {
        result.emplace_back(nodeAndMutId.first, nodeAndMutId.second);
    }
    return std::move(result);
}

PYBIND11_MODULE(_grgl, m) {
    py::class_<grgl::GRG::NodeListIterator>(m, "GRG_NodeListIterator")
        .def(
            "__iter__",
            [](const grgl::GRG::NodeListIterator& nli) { return py::make_iterator(nli.begin(), nli.end()); },
            py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

    py::class_<grgl::NodeData>(m, "NodeData")
        .def_readonly("num_individual_coals", &grgl::NodeData::numIndividualCoals, R"^(
            The number of individuals that coalesce at this node.
            )^")
        .def_readwrite("population_id", &grgl::NodeData::populationId, R"^(
            The population ID (integer) for the population associated with this node.
            )^");

    py::class_<grgl::Mutation>(m, "Mutation")
        .def(py::init<double, std::string, const std::string&, double>(),
             py::arg("position"),
             py::arg("allele"),
             py::arg("ref_allele") = "",
             py::arg("time") = -1.0,
             R"^(
                Construct a new Mutation object, to use as a lookup key or to add to a GRG.
             )^")
        .def_property_readonly("allele", &grgl::Mutation::getAllele, R"^(
            (Read-only) Allele value associated with the Mutation. Can be a single nucleotide
            or a sequence of them.
        )^")
        .def_property_readonly("ref_allele", &grgl::Mutation::getRefAllele, R"^(
            (Read-only) Reference allele at the position that this Mutation occurs. Can be
            empty string if not provided.
        )^")
        .def_property_readonly("position", &grgl::Mutation::getPosition, R"^(
            (Read-only) Position in the genome. Can be absolute or relative (genomic-distance based or
            otherwise normalized).
        )^")
        .def_property("time", &grgl::Mutation::getTime, &grgl::Mutation::setTime, R"^(
            (Read/write) Time value associated with the Mutation, or -1.0 if unused.
        )^")
        .def(pybind11::self == pybind11::self)
        .def(pybind11::self < pybind11::self)
        .def("__hash__", &hashMutation);

    py::class_<grgl::GRG, std::shared_ptr<grgl::GRG>> grgClass(m, "GRG");
    grgClass
        .def("is_sample", &grgl::GRG::isSample, R"^(
                Returns true if the given NodeID is associated with a sample.

                :param node_id: The NodeID to check.
                :type node_id: int
                :return: True iff it is a sample node.
                :rtype: bool
            )^")
        .def_property_readonly("num_samples", &grgl::GRG::numSamples, R"^(
                The number of sample nodes in the GRG.
            )^")
        .def_property_readonly("nodes_are_ordered", &grgl::GRG::nodesAreOrdered, R"^(
                Returns true if the NodeIDs are already in topological order from
                the bottom-up. If this is true, then the first `S` NodeIDs starting
                at 0 will be the sample Nodes, and then the next `N-S` NodeIDs will
                be in order as emitted by a DFS of the GRG starting from the roots
                and emitting the NodeIDs in post-order.

                If this is true, you can often just iterate the NodeIDs from 0...
                num_nodes instead of performing actual graph traversals, depending
                on what you are trying to accomplish.
            )^")
        .def_property_readonly("mutations_are_ordered", &grgl::GRG::mutationsAreOrdered, R"^(
                Returns true if the MutationID order matches the (position, allele)
                sorted order. That is, MutationID of 0 is the lowest position value
                and MutationID of num_mutations-1 is the highest. Ties are broken
                by the lexicographic order of the allele.
            )^")
        .def_property_readonly("num_nodes", &grgl::GRG::numNodes, R"^(
                Get the total number of nodes (including sample and mutation nodes)
                in the GRG.
            )^")
        .def_property_readonly("num_edges", &grgl::GRG::numEdges, R"^(
                Return the total number of down edges in the graph. Down and
                up edges are always symmetric, so the count is the same.
            )^")
        .def("num_up_edges", &grgl::GRG::numUpEdges, R"^(
                Count the number of parents. This can be more efficient than
                getting the list of parents and computing the length.

                :param node_id: The NodeID to get edge count for.
                :type node_id: int
                :return: The number of up edges (parents) for the node..
                :rtype: int
            )^")
        .def("num_down_edges", &grgl::GRG::numDownEdges, R"^(
                Count the number of children. This can be more efficient than
                getting the list of children and computing the length.

                :param node_id: The NodeID to get edge count for.
                :type node_id: int
                :return: The number of down edges (children) for the node..
                :rtype: int
            )^")
        .def("get_down_edges", &grgl::GRG::getDownEdges, R"^(
                Get a list of NodeIDs that are connected to the given NodeID,
                via "down" edges (i.e., children).

                :param node_id: The NodeID to get children for.
                :type node_id: int
                :return: The children of the given node as a list of NodeIDs.
                :rtype: List[int]
            )^")
        .def("get_up_edges", &grgl::GRG::getUpEdges, R"^(
                Get a list of NodeIDs that are connected to the given NodeID,
                via "up" edges (i.e., parents).

                :param node_id: The NodeID to get parents for.
                :type node_id: int
                :return: The parents of the given node as a list of NodeIDs.
                :rtype: List[int]
            )^")
        .def("get_node_data", &grgl::GRG::getNodeData, R"^(
                Get the NodeData object associated with the NodeID.

                :param node_id: The NodeID to get data for.
                :type node_id: int
                :return: The NodeData object.
                :rtype: pygrgl.NodeData
            )^")
        .def("get_sample_nodes", &grgl::GRG::getSampleNodes, R"^(
                Get the NodeIDs for the sample nodes.

                :return: The list of NodeIDs that are sample nodes.
                :rtype: List[int]
            )^")
        .def("get_root_nodes", &grgl::GRG::getRootNodes, R"^(
                Get the NodeIDs for nodes that have no up edges: the roots of the GRG.

                :return: The list of NodeIDs that are root nodes.
                :rtype: List[int]
            )^")
        .def("get_node_mutation_pairs", &getNodeMutationPairs, R"^(
                Get a list of pairs (NodeID, MutationID). Each Mutation typically
                is associated to a single Node, but rarely it can have more than one
                Node, in which case it will show up in more than one pair.
                Results are ordered by NodeID, ascending.

                :return: A list of pairs of NodeID and MutationID.
                :rtype: List[Tuple[int, int]]
            )^")
        .def("get_mutation_node_pairs", &grgl::GRG::getMutationsToNodeOrdered, R"^(
                Get a list of pairs (MutationID, NodeID). Each Mutation typically
                is associated to a single Node, but rarely it can have more than one
                Node, in which case it will show up in more than one pair.
                Results are ordered by MutationID, ascending.

                :return: A list of pairs of MutationID and NodeID.
                :rtype: List[Tuple[int, int]]
            )^")
        .def("get_mutations_for_node", &grgl::GRG::getMutationsForNode, py::arg("node_id"), R"^(
                Get all the (zero or more) Mutations associated with the given NodeID.

                :param node_id: The NodeID to get mutations for.
                :type node_id: int
                :return: A list of MutationIDs.
                :rtype: List[int]
            )^")
        .def("get_mutation_by_id", &grgl::GRG::getMutationById, py::arg("mut_id"), R"^(
                Get the Mutation associated with the given MutationID.

                :param mut_id: The MutationID to get the Mutation for.
                :type mut_id: int
                :return: The mutation.
                :rtype: pygrgl.Mutation
            )^")
        .def("set_mutation_by_id", &grgl::GRG::setMutationById, py::arg("mut_id"), py::arg("mutation"), R"^(
                Set the Mutation associated with the given MutationID.

                :param mut_id: The MutationID to get the Mutation for.
                :type mut_id: int
                :param mutation: The mutation. Users can associate whatever Mutation they want with a particular ID,
                    but usually this is the same as the previous Mutation at this ID, with some non-essential
                    properties changes, like "time".
                :type: pygrgl.Mutation
            )^")

        .def("node_has_mutations", &grgl::GRG::nodeHasMutations, py::arg("node_id"), R"^(
                Return true if there is one or more Mutations associated with the given
                NodeID.

                :param node_id: The NodeID to check for mutations.
                :type node_id: int
                :return: True if the node has at least one mutation.
                :rtype: bool
            )^")
        .def("add_population", &grgl::GRG::addPopulation, py::arg("pop_desc"), R"^(
                Add a new population to the GRG, and return the ID associated with it.

                :param pop_desc: The population description/name.
                :type pop_desc: str
                :return: The PopulationID.
                :rtype: int
            )^")
        .def("get_populations", &grgl::GRG::getPopulations, R"^(
                Get the (possibly empty) list of population descriptions for this GRG.

                :return: The population descriptions.
                :rtype: List[str]
            )^")
        .def("add_mutation", &grgl::GRG::addMutation, py::arg("mutation"), py::arg("node_id"), R"^(
                Add a new Mutation to the GRG, and associate it with the given NodeID.

                :param mutation: The Mutation object.
                :type mutation: pygrgl.Mutation
                :param node_id: The NodeID to attach the Mutation to.
                :type node_id: int
            )^")
        .def_property_readonly("num_mutations", &grgl::GRG::numMutations, R"^(
                Get the total number of mutations in the GRG.

                :return: The mutation count.
                :rtype: int
            )^");
    grgClass.doc() = "A Genotype Representation Graph (GRG) representing a particular dataset. "
                     "This is the immutable portion of the API, so every graph has these operations. "
                     "See MutableGRG for an extension of this that includes the ability to add/remove nodes "
                     "and edges from the graph.";

    py::enum_<grgl::TraversalDirection>(m, "TraversalDirection")
        .value("DOWN", grgl::TraversalDirection::DIRECTION_DOWN, R"^(
            Traverse the graph "down" edges.
        )^")
        .value("UP", grgl::TraversalDirection::DIRECTION_UP, R"^(
            Traverse the graph via "up" edges.
        )^")
        .export_values();

    py::class_<grgl::MutableGRG, std::shared_ptr<grgl::MutableGRG>>(m, "MutableGRG", grgClass)
        .def(py::init<size_t, size_t>(), R"^()^")
        .def("make_node", &grgl::MutableGRG::makeNode, py::arg("count") = 1, R"^(
            Create one or more new nodes in the graph.

            :param count: How many nodes to create (optional, default to 1).
            :type count: int
        )^")
        .def("connect", &grgl::MutableGRG::connect, py::arg("source"), py::arg("target"), R"^(
            Add a down edge from source to target, and an up edge from target to source.

            :param source: The NodeID for the source node (edge starts here).
            :type source: int
            :param target: The NodeID for the target node (edge ends here).
            :type target: int
        )^")
        .def("disconnect", &grgl::MutableGRG::disconnect, py::arg("source"), py::arg("target"), R"^(
            Remove the down edge from source to target, and the up edge from target to source.

            :param source: The NodeID for the source node (edge starts here).
            :type source: int
            :param target: The NodeID for the target node (edge ends here).
            :type target: int
        )^")
        .def("merge", &grgl::MutableGRG::merge, R"^(
            One or more GRG files into the current GRG.
        )^");

    m.def("load_mutable_grg", &grgl::loadMutableGRG, py::arg("filename"), R"^(
        Load a GRG file from disk. Mutable GRGs can have nodes and edges added/removed
        from them.

        :param filename: The file to load.
        :type filename: str
        :return: The GRG.
        :rtype: pygrgl.MutableGRG
    )^");

    m.def("load_immutable_grg",
          &grgl::loadImmutableGRG,
          py::arg("filename"),
          py::arg("load_up_edges") = true,
          py::arg("load_down_edges") = true,
          R"^(
        Load a GRG file from disk. Immutable GRGs are much faster to traverse than mutable
        GRGs and take up less RAM, so this is the preferred method if you are using a GRG
        for calculation or annotation, and not modifying the graph structure itself.

        :param filename: The file to load.
        :type filename: str
        :param load_up_edges: If False, do not load the graph "up" edges (saves RAM).
        :type load_up_edges: bool
        :param load_down_edges: If False, do not load the graph "down" edges (saves RAM).
        :type load_down_edges: bool
        :return: The GRG.
        :rtype: pygrgl.GRG
    )^");

    m.def("save_grg", &grgl::saveGRG, py::arg("grg"), py::arg("filename"), py::arg("allow_simplify") = true, R"^(
        Save the GRG to disk, simplifying it (if possible) in the process.

        :param filename: The file to save to.
        :type filename: str
    )^");

    m.def("grg_from_trees", &grgl::grgFromTrees, py::arg("filename"), py::arg("binary_mutations") = false, R"^(
        Convert a .trees (TSKit tree-sequence) file to a GRG.

        :param filename: The tree-sequence (.trees) file to load.
        :type filename: str
        :param binary_mutations: Set to True to flatten all mutations to be bi-allelic (optional).
        :type binary_mutations: bool
        :return: The GRG.
        :rtype: pygrgl.GRG
    )^");

    m.def("get_bfs_order",
          &getBfsOrder,
          py::arg("grg"),
          py::arg("direction"),
          py::arg("seed_list"),
          py::arg("max_queue_width") = -1,
          R"^(
        Get a list of NodeIDs in breadth-first-search (BFS) order, starting from the given
        seeds and traversing in the provided TraversalDirection (up or down).

        :param grg: The GRG to get nodes for.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param direction: The direction to traverse, up or down.
        :type direction: pygrgl.TraversalDirection
        :param seed_list: The list of NodeIDs that represent the starting place of the traversal. For
            example, if you use pygrgl.GRG.get_sample_nodes() and pygrgl.TraversalDirection.UP then
            the entire graph will be traversed from bottom to top.
        :type seed_list: List[int]
        :param max_queue_width: The maximum width the queue used for bread-first-search. The default
            is -1, which means there is no maximum width. Setting this can help reduce traversal cost
            but will result in an incomplete traversal.
        :type max_queue_width: int
        :return: The ordered list of NodeIDs.
        :rtype: List[int]
    )^");

    m.def("get_dfs_order",
          &getDfsOrder,
          py::arg("grg"),
          py::arg("direction"),
          py::arg("seed_list"),
          py::arg("forward_only") = false,
          R"^(
        Get a list of NodeIDs in depth-first-search (DFS) order, starting from the given
        seeds and traversing in the provided TraversalDirection (up or down).

        :param grg: The GRG to get nodes for.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param direction: The direction to traverse, up or down.
        :type direction: pygrgl.TraversalDirection
        :param seed_list: The list of NodeIDs that represent the starting place of the traversal. For
            example, if you use pygrgl.GRG.get_sample_nodes() and pygrgl.TraversalDirection.UP then
            the entire graph will be traversed from bottom to top.
        :type seed_list: List[int]
        :param forward_only: If True, enumerates nodes in the given direction from seeds and outputs
            them in the order they are first visited. If False, enumerated nodes in the order they are
            visited the _second_ time, i.e. when they are popped off the stack that is used for the
            depth-first-search. This provides a topological order to the nodes when this parameter is
            set to False, and does _not_ when set to True. Default is False.
        :type forward_only: bool
        :return: The ordered list of NodeIDs.
        :rtype: List[int]
    )^");

    m.def("get_topo_order", &getTopoOrder, py::arg("grg"), py::arg("direction"), py::arg("seed_list"), R"^(
        Get a list of NodeIDs in topological order, starting from the given
        seeds and traversing in the provided TraversalDirection (up or down).
        This order is similar to DFS order, but with the seeds being the "endpoints".
        The NodeIDs in this order are guaranteed to have the property: if a NodeID
        X is at position P in the list, then all of the nodes above/below X
        (depending on TraversalDirection) have positions before P.

        :param grg: The GRG to get nodes for.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param direction: The direction to traverse, up or down.
        :type direction: pygrgl.TraversalDirection
        :param seed_list: The list of NodeIDs that represent the starting place of the traversal. For
            example, if you use pygrgl.GRG.get_sample_nodes() and pygrgl.TraversalDirection.UP then
            the entire graph will be traversed from bottom to top.
        :type seed_list: List[int]
        :return: The ordered list of NodeIDs.
        :rtype: List[int]
    )^");

    m.def("dot_product", &dotProduct, py::arg("grg"), py::arg("input"), py::arg("direction"), R"^(
        Compute one of two possible dot products across the entire graph. The input vector :math:`V` can be
        either :math:`1 \times N` (:math:`N` is number of samples) or :math:`1 \times M` (:math:`M` is number of
        mutations). The given direction determines which input vector is expected. Let :math:`X` be the
        :math:`N \times M` genotype matrix.
        For an :math:`1 \times N` input :math:`V`, the product performed is :math:`V \cdot X` which gives a
        :math:`1 \times M` result. I.e., the input vector is a value per sample and the output vector is
        a value per mutation.
        For an :math:`1 \times M` input :math:`V`, the product performed is :math:`V \cdot X^T` which gives a
        :math:`1 \times N` result. I.e., the input vector is a value per mutation and the output vector is
        a value per sample.

        Dot product in the graph works by seeding the input nodes (samples or mutations) with the corresponding
        values from the input vector and then traversing the graph in the relevant direction (up or down). The
        ancestor/descendant values are summed at each node, until the terminal nodes (mutations or samples) are
        reached.

        :param grg: The GRG to perform the computation against.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param input: The numpy array of input values :math:`V`.
        :type seed_list: numpy.array
        :param direction: The direction to traverse, up (input is per sample) or down (input is per mutation).
        :type direction: pygrgl.TraversalDirection
        :return: The numpy array of output values.
        :rtype: numpy.array

    )^");

    m.attr("INVALID_NODE") = INVALID_NODE_ID;
    m.attr("COAL_COUNT_NOT_SET") = grgl::NodeData::COAL_COUNT_NOT_SET;

    m.attr("__version__") = "dev";
}
