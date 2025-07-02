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

#include "common_visitors.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/grg.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/node_data.h"
#include "grgl/serialize.h"
#include "grgl/ts2grg.h"
#include "grgl/version.h"
#include "grgl/visitor.h"

namespace py = pybind11;

#include <iostream>

template <typename T>
inline py::array_t<T> dispatchDotProd(
    grgl::GRGPtr& grg, py::buffer_info& buffer, size_t inCols, size_t outCols, grgl::TraversalDirection direction) {
    py::array_t<T> result(outCols);
    py::buffer_info resultBuf = result.request();
    memset(resultBuf.ptr, 0, outCols * sizeof(T));
    grg->matrixMultiplication<T, T, false>((const T*)buffer.ptr, inCols, 1, direction, (T*)resultBuf.ptr, outCols);
    return std::move(result);
}

py::array dotProduct(grgl::GRGPtr& grg, py::handle input, grgl::TraversalDirection direction) {
    py::array arr = py::array::ensure(input, py::array::c_style | py::array::forcecast);
    py::buffer_info buffer = arr.request();
    if (buffer.ndim != 1) {
        throw grgl::ApiMisuseFailure("dot_product() only support one-dimensional numpy arrays as input.");
    }
    const size_t cols = buffer.shape.at(0);
    if (cols == 0) {
        throw grgl::ApiMisuseFailure("dot_product() requires non-empty array input.");
    }
    const size_t outCols =
        (direction == grgl::TraversalDirection::DIRECTION_DOWN) ? grg->numSamples() : grg->numMutations();
    if (py::isinstance<py::array_t<double>>(input)) {
        return dispatchDotProd<double>(grg, buffer, cols, outCols, direction);
    } else if (py::isinstance<py::array_t<float>>(input)) {
        return dispatchDotProd<float>(grg, buffer, cols, outCols, direction);
    } else if (py::isinstance<py::array_t<std::int64_t>>(input)) {
        return dispatchDotProd<int64_t>(grg, buffer, cols, outCols, direction);
    } else if (py::isinstance<py::array_t<std::int32_t>>(input)) {
        return dispatchDotProd<int32_t>(grg, buffer, cols, outCols, direction);
    }
    std::stringstream ssErr;
    ssErr << "Unsupported numpy dtype: " << arr.dtype();
    throw grgl::ApiMisuseFailure(ssErr.str().c_str());
}

template <typename IOType, typename NodeValueType, bool useBitVector>
inline py::array_t<IOType> dispatchMult(grgl::GRGPtr& grg,
                                        py::buffer_info& buffer,
                                        size_t rows,
                                        size_t inCols,
                                        size_t outCols,
                                        grgl::TraversalDirection direction,
                                        bool emitAllNodes,
                                        bool byIndividual,
                                        py::handle init) {
    const size_t outSize = rows * outCols;
    py::array_t<IOType> result({rows, outCols});
    py::buffer_info resultBuf = result.request();
    memset(resultBuf.ptr, 0, outSize * sizeof(IOType));

    py::buffer_info initBuffer;
    grgl::GRG::NodeInitEnum nodeInitMode = grgl::GRG::NIE_ZERO;
    if (py::isinstance<py::none>(init)) {
        ;
    } else if (py::isinstance<py::str>(init)) {
        if (init.cast<std::string>() == "xtx") {
            nodeInitMode = grgl::GRG::NIE_XTX;
        } else {
            std::stringstream ssErr;
            ssErr << "Unexpected init value: " << init.cast<std::string>();
            throw grgl::ApiMisuseFailure(ssErr.str().c_str());
        }
    } else {
        py::array arr = py::array::ensure(init, py::array::c_style | py::array::forcecast);
        if (py::isinstance<py::array_t<IOType>>(init)) {
            initBuffer = arr.request();
            if (initBuffer.ndim == 1) {
                const size_t initRows = initBuffer.shape.at(0);
                if (initRows != rows) {
                    std::stringstream ssErr;
                    ssErr << "If init has a single dimension, it must match the number of rows in the input matrix";
                    throw grgl::ApiMisuseFailure(ssErr.str().c_str());
                }
                nodeInitMode = grgl::GRG::NIE_VECTOR;
            }
            if (initBuffer.ndim == 2) {
                const size_t initRows = initBuffer.shape.at(0);
                const size_t initCols = initBuffer.shape.at(1);
                if (initRows != rows || initCols != grg->numNodes()) {
                    std::stringstream ssErr;
                    ssErr << "If init is a matrix, it must match the dimensions ROW x NODES";
                    throw grgl::ApiMisuseFailure(ssErr.str().c_str());
                }
                nodeInitMode = grgl::GRG::NIE_MATRIX;
            }
        } else {
            std::stringstream ssErr;
            ssErr << "The init matrix must match the dtype of the input matrix. Got: " << arr.dtype();
            throw grgl::ApiMisuseFailure(ssErr.str().c_str());
        }
    }

    grg->matrixMultiplication<IOType, NodeValueType, useBitVector>((const IOType*)buffer.ptr,
                                                                   inCols,
                                                                   rows,
                                                                   direction,
                                                                   (IOType*)resultBuf.ptr,
                                                                   outSize,
                                                                   emitAllNodes,
                                                                   byIndividual,
                                                                   (const IOType*)initBuffer.ptr,
                                                                   nodeInitMode);
    return std::move(result);
}

py::array matMul(grgl::GRGPtr& grg,
                 py::handle input,
                 grgl::TraversalDirection direction,
                 bool emitAllNodes = false,
                 bool byIndividual = false,
                 py::handle init = nullptr) {
    py::array arr = py::array::ensure(input, py::array::c_style | py::array::forcecast);
    py::buffer_info buffer = arr.request();
    if (buffer.ndim != 2) {
        throw grgl::ApiMisuseFailure("matmul() only support two dimensional numpy arrays as input.");
    }
    const size_t rows = buffer.shape.at(0);
    const size_t cols = buffer.shape.at(1);
    if (rows == 0 || cols == 0) {
        throw grgl::ApiMisuseFailure("matmul() requires non-zero dimensions.");
    }
    const size_t numSamples = byIndividual ? grg->numIndividuals() : grg->numSamples();
    const size_t outCols =
        (emitAllNodes ? grg->numNodes()
                      : ((direction == grgl::TraversalDirection::DIRECTION_DOWN) ? numSamples : grg->numMutations()));

    if (py::isinstance<py::array_t<double>>(input)) {
        return dispatchMult<double, double, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<float>>(input)) {
        return dispatchMult<float, float, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::int64_t>>(input)) {
        return dispatchMult<int64_t, int64_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::int32_t>>(input)) {
        return dispatchMult<int32_t, int32_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::int16_t>>(input)) {
        return dispatchMult<int16_t, int16_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::int8_t>>(input)) {
        return dispatchMult<int8_t, int8_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::uint64_t>>(input)) {
        return dispatchMult<uint64_t, uint64_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::uint32_t>>(input)) {
        return dispatchMult<uint32_t, uint32_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::uint16_t>>(input)) {
        return dispatchMult<uint16_t, uint16_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<std::uint8_t>>(input)) {
        return dispatchMult<uint8_t, uint8_t, false>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    } else if (py::isinstance<py::array_t<bool>>(input)) {
        if (!py::isinstance<py::none>(init)) {
            throw grgl::ApiMisuseFailure("init=... is not supported with boolean dtype");
        }
        return dispatchMult<bool, uint8_t, true>(
            grg, buffer, rows, cols, outCols, direction, emitAllNodes, byIndividual, init);
    }
    std::stringstream ssErr;
    ssErr << "Unsupported numpy dtype: " << arr.dtype();
    throw grgl::ApiMisuseFailure(ssErr.str().c_str());
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

// seedList is assumed to be unique! The frontier will be empty if you include duplicates in the
// seedList.
std::vector<grgl::NodeID>
sharedFrontier(const grgl::GRGPtr& grg, grgl::TraversalDirection direction, const grgl::NodeIDList& seedList) {
    grgl::FrontierVisitor visitor(seedList);
    grg->visitTopo(visitor, direction, seedList);
    return std::move(visitor.m_frontier);
}

PYBIND11_MODULE(_grgl, m) {
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
        .def_property_readonly("num_individuals", &grgl::GRG::numIndividuals, R"^(
                The number of individuals in the GRG. The corresponding samples can be found by
                multiplying the individual index :math:`I` by the ploidy:
                :math:`(ploidy * I) + 0, (ploidy \times I) + 1, ..., (ploidy \times I) + (ploidy - 1)`
            )^")
        .def_property_readonly("ploidy", &grgl::GRG::getPloidy, R"^(
                The ploidy of each individual.
            )^")
        .def_property_readonly("bp_range", &grgl::GRG::getBPRange, R"^(
                The range in base-pair positions that this GRG covers, from its list of mutations.

                A pair (lower, upper) of the range covered by this GRG, where lower is inclusive
                and upper is exclusive.
            )^")
        .def_property_readonly("specified_bp_range", &grgl::GRG::getSpecifiedBPRange, R"^(
                The range in base-pair positions that this GRG covers, as specified during the
                GRG construction. This range may exceed the bp_range if the GRG was constructed
                from a range of the genome that did not immediately start/end with a mutation.

                A pair (lower, upper) of the range covered by this GRG, where lower is inclusive
                and upper is exclusive.
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
        .def("get_node_mutation_pairs", &grgl::GRG::getNodeMutationPairs, R"^(
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
            )^")
        .def("get_population_id", &grgl::GRG::getPopulationId, py::arg("node_id"), R"^(
                Get the population ID for the given node. Can be used to index into the list returned by
                get_populations().

                :param node_id: The node to retrieve.
                :type node_id: int
                :return: The population ID.
                :rtype: int
            )^")
        .def("set_population_id", &grgl::GRG::setPopulationId, py::arg("node_id"), py::arg("population_id"), R"^(
                Set the population ID for the given node.

                :param node_id: The node to access.
                :type node_id: int
                :param population_id: The population ID to associate with the node.
                :type population_id: int
            )^")
        .def("get_num_individual_coals", &grgl::GRG::getNumIndividualCoals, py::arg("node_id"), R"^(
                Get the number of individuals that coalesced at the given node (not below or above).

                :param node_id: The node to retrieve.
                :type node_id: int
                :return: The number of individuals that coalesced, or pygrgl.COAL_COUNT_NOT_SET.
                :rtype: int
            )^")
        .def("set_num_individual_coals",
             &grgl::GRG::setNumIndividualCoals,
             py::arg("node_id"),
             py::arg("num_coals"),
             R"^(
                Set the number of individuals that coalesced at the given node (not below or above).

                :param node_id: The node to access.
                :type node_id: int
                :param num_coals: The number of individuals that coalesced, or pygrgl.COAL_COUNT_NOT_SET.
                :type num_coals: int
            )^")
        .def_property_readonly("has_individual_ids", &grgl::GRG::hasIndividualIds, R"^(
            True if this GRG has string identifiers for each individual. See get_individual_id().
            )^")
        .def("clear_individual_ids", &grgl::GRG::clearIndividualIds, R"^(
            Remove all individual IDs from the current GRG.
            )^")
        .def("add_individual_id", &grgl::GRG::addIndividualId, py::arg("identifier"), R"^(
            Add the next string identifier for an individual in the dataset. If the individual IDs
            are already set, this will throw an exception. This must be called in order, from the
            0th to the (N-1)st individual.

            :param identifier: The string identifier for the next individual.
            :type identifier: str
            )^")
        .def("get_individual_id", &grgl::GRG::getIndividualId, py::arg("individual_index"), R"^(
            Get the string identifiers for each of the N individuals in the dataset, if available.
            These are optional, so the empty list will be returned if the GRG does not have them.
            See has_individual_ids.

            :param individual_index: The individual to retrieve. The individuals are numbered from
                0...(num_individuals-1), and correspond to the sample NodeIDs divided by their ploidy.
            :type individual_index: int
            :return: String identifier for the given individual.
            :rtype: str
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
        .def("make_node", &grgl::MutableGRG::makeNode, py::arg("count") = 1, py::arg("force_ordered") = false, R"^(
            Create one or more new nodes in the graph.

            :param count: How many nodes to create (optional, default to 1).
            :type count: int
            :param force_ordered: Set to True if you are sure adding this node will maintain the topological
                order of the NodeIDs within the graph.
            :type count: bool
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
        .def("merge", &grgl::MutableGRG::merge, py::arg("other_grg_files"), py::arg("combine_nodes"), R"^(
            Merge one or more GRGs into this one. Only succeeds if all GRGs have the same number of

            This assumes that the GRGs were constructed from the same sampleset -- e.g., they
            could be constructed from two subsets of the same sampleset (as long as both were
            constructed with the same sample node numbering) or from a subset of mutations against
            the same sampleset.

            The specified range of the resulting GRG will be (min(range of any input), max(range of any input)),
            even if the provided GRGs do not span a contiguous region. It is up to the caller of this
            method to ensure that either (1) the span is contiguous or (2) they adjust the specified range
            appropriately afterwards.

            :param other_grg_files: A list of filenames for the GRG to merge into the current one.
            :type other_grg_files: List[str]
            :param combine_nodes: True by default. Combine nodes from different GRGs if the node has
                the same samples beneath it.
            :type combine_nodes: bool
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
          R"^(
        Load a GRG file from disk. Immutable GRGs are much faster to traverse than mutable
        GRGs and take up less RAM, so this is the preferred method if you are using a GRG
        for calculation or annotation, and not modifying the graph structure itself.

        :param filename: The file to load.
        :type filename: str
        :param load_up_edges: If False, do not load the graph "up" edges (saves RAM).
        :type load_up_edges: bool
        :return: The GRG.
        :rtype: pygrgl.GRG
    )^");

    m.def("save_grg", &grgl::saveGRG, py::arg("grg"), py::arg("filename"), py::arg("allow_simplify") = true, R"^(
        Save the GRG to disk, simplifying it (if possible) in the process.

        :param grg: The GRG
        :type filename: pygrgl.GRG
        :param filename: The file to save to.
        :type filename: str
        :param allow_simplify: Set to False to disallow removing nodes/edges from the graph that do not
            significantly contribute to the mutation-to-samples mapping.
        :type allow_simplify: bool
    )^");

    m.def("save_subset",
          &grgl::saveGRGSubset,
          py::arg("grg"),
          py::arg("filename"),
          py::arg("direction"),
          py::arg("seed_list"),
          py::arg("bp_range") = std::pair<grgl::BpPosition, grgl::BpPosition>(),
          R"^(
        Save a subset of the GRG to disk, specified by a vector masking either mutation IDs or sample IDs.

        :param grg: The GRG
        :type filename: pygrgl.GRG
        :param filename: The file to save to.
        :type filename: str
        :param direction: Downward means the seeds should be a list of MutationID that should be kept in
            the graph. Upward means the seeds should be a list of sample NodeID that should be kept.
        :type direction: pygrgl.TraversalDirection
        :param seed_list: A list of MutationID or NodeID (see direction parameter).
        :type seed_list: List[int]
        :param bp_range: A pair of integers specifying the base-pair range that this GRG covers. This is just
            meta-data, and does not change the filtering behavior.
        :type bp_range: Tuple[int, int]
    )^");

    m.def("grg_from_trees",
          &grgl::grgFromTrees,
          py::arg("filename"),
          py::arg("binary_mutations") = false,
          py::arg("use_node_times") = false,
          py::arg("maintain_topology") = false,
          py::arg("compute_coals") = false,
          R"^(
        Convert a .trees (TSKit tree-sequence) file to a GRG.

        :param filename: The tree-sequence (.trees) file to load.
        :type filename: str
        :param binary_mutations: Set to True to flatten all mutations to be bi-allelic (optional).
        :type binary_mutations: bool
        :param use_node_times: Mutations will be assigned the time from the node below them, instead of
            the tskit Mutation object.
        :type use_node_times: bool
        :param maintain_topology: Generates a slightly larger GRG, but ensures that we capture all tree
            topology changes induced by recombination, not just the changes that result in a different
            set of samples beneath a node.
        :type maintain_topology: bool
        :param compute_coals: Compute the per-node coalescence counts. I.e., how many individuals coalesced
            exactly at the node (separate children have both haploid copies of the individual). This is an
            expensive computation, slowing down the TS to GRG conversion when there are a lot of samples.
        :type compute_coals: bool
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
        DEPRECATED, use :py:meth:`matmul` instead.

        Special case of :py:meth:`matmul`, where input is a vector (1-dimensional numpy array) instead
        of a matrix. Computes the vector-matrix product between the input vector and the genotype matrix
        :math:`X` (or :math:`X^T`). See :py:meth:`matmul` for more details.

        :param grg: The GRG to perform the computation against.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param input: The numpy 1-dimensional array of input values :math:`V`.
        :type input: numpy.array
        :param direction: The direction to traverse, up (input is per sample) or down (input is per mutation).
        :type direction: pygrgl.TraversalDirection
        :return: The numpy 1-dimensional array of output values.
        :rtype: numpy.array

    )^");

    m.def("matmul",
          &matMul,
          py::arg("grg"),
          py::arg("input"),
          py::arg("direction"),
          py::arg("emit_all_nodes") = false,
          py::arg("by_individual") = false,
          py::arg("init") = nullptr,
          R"^(
        Compute one of two possible matrix multiplications across the entire
        graph. The input matrix :math:`V` can be either :math:`K \times N` (:math:`N`
        is number of samples) or :math:`K \times M` (:math:`M` is number of
        mutations). The given direction determines which input matrix is
        expected. Let :math:`X` be the :math:`N \times M` genotype matrix. For an
        :math:`K \times N` input :math:`V`, the product performed is :math:`V \times X`
        which gives a :math:`K \times M` result. I.e., the input matrix is a column
        per sample and the output matrix is a column per mutation. For an :math:`K
        \times M` input :math:`V`, the product performed is :math:`V \times X^T`
        which gives a :math:`K \times N` result. I.e., the input matrix is a column
        per mutation and the output matrix is a column per sample.

        The simplest case to consider is a vector input (e.g., a :math:`1 \times N`
        matrix). This vector-matrix product in the graph works by seeding the
        input nodes (samples in this example) with the corresponding values from
        the input vector and then traversing the graph in the relevant direction
        (up or down). The ancestor/descendant values are summed at each node,
        until the terminal nodes (mutations in this example) are reached. The
        values at the terminal nodes are then the output vector.
        When a :math:`K`-row matrix is input, instead of a vector, the only
        difference is that each node stores :math:`K` values instead of 1.

        Note: the RAM used will be :math:`O(K * nodes)` where :math:`nodes` is the
        total number of nodes in the graph.

        :param grg: The GRG to perform the computation against.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param input: The numpy 2-dimensional array of input values :math:`V`.
        :type input: numpy.array
        :param direction: The direction to traverse, up (input is per sample) or down (input is per mutation).
        :type direction: pygrgl.TraversalDirection
        :param emit_all_nodes: False by default. Set to True if you want each output row in the matrix to
            have a value for every node, not just every sample/mutation (depending on direction).
        :type emit_all_nodes: bool
        :param by_individual: The dimension that is for samples (either the input or output, depending on the
            direction parameter) uses individuals instead of haploid samples. Instead of outputting vectors
            of :math:`N` (:py:attr:`num_samples`) columns, it is :math:`N / ploidy` (:py:attr:`num_individuals`)
            columns.
        :type by_individual: bool
        :param init: Initialization of the nodes of the graph during matrix multiplication. By default (when this
            is set to None), nodes are initialization to 0. There are three possible types this can take on:
            1. A string "xtx" which means to initialize the nodes with twice their coalesence counts. Using this
            and performing an UP multiplication (with 1s as input) produces the X.T * X product needed for
            GWAS.
            2. A one dimensional numpy array (vector) of length K. The value at position K is assigned to all
            nodes when performing the multiplication for row K from the input matrix.
            3. A two dimensional numpy array (matrix) of size KxT, where T is the total number of nodes in the
            graph (grg.num_nodes). This fully specifies every node value for the entire matrix operation.
        :type init: Union[str, numpy.array]
        :return: The numpy 2-dimensional array of output values.
        :rtype: numpy.array
    )^");

    m.def("shared_frontier", &sharedFrontier, py::arg("grg"), py::arg("direction"), py::arg("seeds"), R"^(
        Get the list of nodes that corresponds to the shared frontier in the graph
        in the given direction. The frontier are the _first_ nodes, along each path,
        that are reached by all seeds nodes in the input.

        Do not pass duplicate Node IDs in "seeds"! There is no checking for duplication,
        but you will get incorrect results if you pass in duplicates.

        :param grg: The GRG to perform the computation against.
        :type grg: pygrgl.GRG or pygrgl.MutableGRG
        :param direction: The direction to traverse, up (terminate at mutations or shared nodes) or
            down (terminate at samples or shared node).
        :type direction: pygrgl.TraversalDirection
        :param seeds: List of node IDs to start the search from.
        :type seeds: List[int]
        :return: A list of node IDs representing the frontier.
        :rtype: List[int]
    )^");

    m.attr("INVALID_NODE") = grgl::INVALID_NODE_ID;
    m.attr("COAL_COUNT_NOT_SET") = grgl::COAL_COUNT_NOT_SET;

    std::stringstream versionString;
    versionString << GRGL_MAJOR_VERSION << "." << GRGL_MINOR_VERSION;
    m.attr("__version__") = versionString.str();
    std::stringstream fileVersionString;
    fileVersionString << grgl::GRG_FILE_MAJOR_VERSION << "." << grgl::GRG_FILE_MINOR_VERSION;
    m.attr("GRG_FILE_VERSION") = fileVersionString.str();
}
