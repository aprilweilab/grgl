#include "grgl/polarization.h"

#include "grgl/map_mutations.h"

#include <vector>

namespace grgl {


// Downward DFS modeled after GRG::visitDfs: stack holds (node, pass) pairs (1=forward, 2=back).
// Graph is acyclic/diamond-free here, so no visited tracking is needed.
static NodeIDList collectSamples(const MutableGRGPtr& grg, NodeID parent, size_t numSamples) {
    NodeIDList samples;
    if (parent == INVALID_NODE_ID) {
        return samples;
    }
    std::vector<std::pair<NodeID, uint8_t>> stack;
    stack.emplace_back(parent, 1);
    while (!stack.empty()) {
        const NodeID nid = stack.back().first;
        const uint8_t pass = stack.back().second;
        if (pass == 1) {
            if (nid < numSamples) {
                samples.push_back(nid);
                stack.pop_back();
                continue;
            }
            stack.back().second = 2;
            for (NodeID child : grg->getDownEdges(nid)) {
                stack.emplace_back(child, 1);
            }
        } else {
            stack.pop_back();
        }
    }
    return samples;
}

std::vector<bool> polarizeMutations(const MutableGRGPtr& grg,
                                    const std::vector<std::pair<MutationId, std::string>>& batch,
                                    PolarizationStats& stats) {
    std::vector<bool> results(batch.size(), false);
    if (batch.empty()) {
        return results;
    }

    const size_t numSamples = grg->numSamples();
    const auto mutNodeList = grg->getMutationsToNodeOrdered<GRG::MutNodeMiss>();

    std::vector<Mutation> remapMutations;
    std::vector<NodeIDList> remapSamples;
    remapMutations.reserve(batch.size());
    remapSamples.reserve(batch.size());

    for (size_t idx = 0; idx < batch.size(); ++idx) {
        stats.totalSeen++;

        const MutationId mutId = batch[idx].first;
        const std::string& ancestralAllele = batch[idx].second;

        NodeID nodeId = INVALID_NODE_ID;
        NodeID missingNodeId = INVALID_NODE_ID;
        for (const auto& entry : mutNodeList) {
            if (std::get<0>(entry) == mutId) {
                nodeId = std::get<1>(entry);
                missingNodeId = std::get<2>(entry);
                break;
            }
        }

        if (nodeId == INVALID_NODE_ID) {
            continue;
        }

        const Mutation& mutation = grg->getMutationById(mutId);
        NodeIDList samples = collectSamples(grg, nodeId, numSamples);
        NodeIDList missingSamples = collectSamples(grg, missingNodeId, numSamples);

        // unknown ancestral allele => drop.
        if (ancestralAllele == "." || ancestralAllele == "-" || ancestralAllele == "N") {
            stats.droppedUnknown++;
            grg->removeMutation(mutId, nodeId);
            continue;
        }

        const std::string ref = mutation.getRefAllele();
        const std::string alt = mutation.getAllele();

        if (ancestralAllele == ref) {
            stats.alreadyPolarized++;
            stats.emitted++;
            results[idx] = true;
            continue;
        }

        if (ancestralAllele == alt) {
            std::vector<bool> sampleSet(numSamples, false);
            for (NodeID sampleId : samples) {
                if (sampleId < numSamples) {
                    sampleSet[sampleId] = true;
                }
            }

            NodeIDList flippedCarriers;
            flippedCarriers.reserve(numSamples - samples.size());
            for (size_t sid = 0; sid < numSamples; sid++) {
                bool isMissing = false;
                for (const NodeID missingId : missingSamples) {
                    if (missingId == sid) {
                        isMissing = true;
                        break;
                    }
                }
                if (!isMissing && !sampleSet[sid]) {
                    flippedCarriers.push_back(static_cast<NodeID>(sid));
                }
            }

            grg->removeMutation(mutId, nodeId);

            stats.swapped++;
            stats.emitted++;
            results[idx] = true;

            remapMutations.emplace_back(mutation.getPosition(), ref, ancestralAllele, mutation.getTime());
            remapSamples.emplace_back(std::move(flippedCarriers));
            continue;
        }

        // inconsistent ancestral allele, drop.
        stats.inconsistent++;
        grg->removeMutation(mutId, nodeId);
    }

    if (!remapMutations.empty()) {
        mapMutations(grg, remapMutations, remapSamples, false, remapMutations.size());
    }

    return results;
}

bool polarizeMutation(const MutableGRGPtr& grg,
                      MutationId mutId,
                      const std::string& ancestralAllele,
                      PolarizationStats& stats) {
    const auto results = polarizeMutations(grg, {{mutId, ancestralAllele}}, stats);
    return !results.empty() && results.front();
}

} // namespace grgl

