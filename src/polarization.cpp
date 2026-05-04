#include "grgl/polarization.h"

#include "grgl/map_mutations.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace grgl {

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

static inline void loadRemaps(const MutableGRGPtr& grg,
                               std::vector<std::pair<MutationId, NodeID>>& removals,
                               std::vector<Mutation>& remapMutations,
                               std::vector<NodeIDList>& remapSamples,
                               size_t flushThreshold,
                               bool force) {
    if (!force && remapMutations.size() < flushThreshold) {
        return;
    }
    grg->getNodesAndMutations();
    for (const auto& removal : removals) {
        grg->removeMutation(removal.first, removal.second);
    }
    removals.clear();
    if (!remapMutations.empty()) {
        mapMutations(grg, remapMutations, remapSamples, false, remapMutations.size());
        remapMutations.clear();
        remapSamples.clear();
    }
}

using MutLookup = std::vector<std::pair<NodeID, NodeID>>;

static MutLookup buildMutLookup(const MutableGRGPtr& grg) {
    const auto mutNodeList = grg->getMutationsToNodeOrdered<GRG::MutNodeMiss>();

    MutLookup mutLookup(grg->numMutations(), {INVALID_NODE_ID, INVALID_NODE_ID});
    for (const auto& entry : mutNodeList) {
        const MutationId mid = std::get<0>(entry);
        if (mid < mutLookup.size()) {
            mutLookup[mid] = {std::get<1>(entry), std::get<2>(entry)};
        }
    }
    return mutLookup;
}

static std::vector<bool> polarizeMutationsHelper(const MutableGRGPtr& grg,
                                                     const MutLookup& mutLookup,
                                                     const std::vector<std::pair<MutationId, std::string>>& batch,
                                                     PolarizationStats& stats) {
    std::vector<bool> results(batch.size(), false);
    if (batch.empty()) {
        return results;
    }

    const size_t numSamples = grg->numSamples();

    std::vector<std::pair<MutationId, NodeID>> removals;
    std::vector<Mutation> remapMutations;
    std::vector<NodeIDList> remapSamples;
    removals.reserve(batch.size());
    remapMutations.reserve(batch.size());
    remapSamples.reserve(batch.size());

    constexpr size_t remapFlushThreshold = 10000;

    for (size_t idx = 0; idx < batch.size(); ++idx) {
        stats.totalSeen++;

        const MutationId mutId = batch[idx].first;
        const std::string& ancestralAllele = batch[idx].second;

        if (mutId >= mutLookup.size()) {
            continue;
        }

        const NodeID nodeId = mutLookup[mutId].first;
        const NodeID missingNodeId = mutLookup[mutId].second;
        if (nodeId == INVALID_NODE_ID) {
            continue;
        }

        const Mutation mutation = grg->getMutationById(mutId);

        // unknown ancestral allele => drop.
        if (ancestralAllele == "." || ancestralAllele == "-" || ancestralAllele == "N") {
            stats.droppedUnknown++;
            removals.emplace_back(mutId, nodeId);
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
            NodeIDList samples = collectSamples(grg, nodeId, numSamples);
            NodeIDList missingSamples = collectSamples(grg, missingNodeId, numSamples);

            std::vector<bool> sampleSet(numSamples, false);
            for (NodeID sampleId : samples) {
                if (sampleId < numSamples) {
                    sampleSet[sampleId] = true;
                }
            }
            std::vector<bool> missingSet(numSamples, false);
            for (NodeID missingId : missingSamples) {
                if (missingId < numSamples) {
                    missingSet[missingId] = true;
                }
            }

            NodeIDList flippedCarriers;
            flippedCarriers.reserve(numSamples - samples.size());
            for (size_t sid = 0; sid < numSamples; sid++) {
                if (!missingSet[sid] && !sampleSet[sid]) {
                    flippedCarriers.push_back(static_cast<NodeID>(sid));
                }
            }

            removals.emplace_back(mutId, nodeId);

            stats.swapped++;
            stats.emitted++;
            results[idx] = true;

            if (!missingSamples.empty()) {
                remapMutations.emplace_back(mutation.getPosition(),
                                            Mutation::ALLELE_MISSING,
                                            mutation.getRefAllele(),
                                            mutation.getTime());
                remapSamples.emplace_back(std::move(missingSamples));
            }

            remapMutations.emplace_back(mutation.getPosition(), ref, ancestralAllele, mutation.getTime());
            remapSamples.emplace_back(std::move(flippedCarriers));
            loadRemaps(grg, removals, remapMutations, remapSamples, remapFlushThreshold, false);
            continue;
        }

        // inconsistent ancestral allele, drop.
        stats.inconsistent++;
        removals.emplace_back(mutId, nodeId);
    }

    loadRemaps(grg, removals, remapMutations, remapSamples, remapFlushThreshold, true);

    return results;
}

std::vector<bool> polarizeMutations(const MutableGRGPtr& grg,
                                    const std::vector<std::pair<MutationId, std::string>>& batch,
                                    PolarizationStats& stats) {
    const MutLookup mutLookup = buildMutLookup(grg);
    auto result = polarizeMutationsHelper(grg, mutLookup, batch, stats);
    grg->sortMutations();
    return result;
}

bool polarizeMutation(const MutableGRGPtr& grg,
                      MutationId mutId,
                      const std::string& ancestralAllele,
                      PolarizationStats& stats) {
    const auto results = polarizeMutations(grg, {{mutId, ancestralAllele}}, stats);
    return !results.empty() && results.front();
}

static std::string loadFasta(const std::string& path) {
    std::ifstream instream(path);
    if (!instream) {
        throw std::runtime_error("Failed to open FASTA: " + path);
    }
    std::string sequence;
    std::string line;
    bool seenHeader = false;
    while (std::getline(instream, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (seenHeader && !sequence.empty()) {
                throw std::runtime_error("Multiple contigs not supported in FASTA: " + path);
            }
            seenHeader = true;
            continue;
        }
        for (char character : line) {
            if (std::isspace(static_cast<unsigned char>(character)) == 0) {
                sequence.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(character))));
            }
        }
    }
    if (sequence.empty()) {
        throw std::runtime_error("No sequence found in FASTA: " + path);
    }
    return sequence;
}

static bool equalsIgnoreCase(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::toupper(static_cast<unsigned char>(a[i])) !=
            std::toupper(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

PolarizationStats polarizeGrgFromFasta(const MutableGRGPtr& grg,
                                       const std::string& fastaPath,
                                       bool dropIfNoMatch,
                                       bool positionsAreOneBased) {
    PolarizationStats stats;
    const std::string fastaSeq = loadFasta(fastaPath);
    const size_t fastaLen = fastaSeq.size();
    const MutLookup mutLookup = buildMutLookup(grg);

    constexpr MutationId chunkSize = 10000;
    std::vector<std::pair<MutationId, std::string>> batch;
    batch.reserve(chunkSize);

    stats.reset();
    const MutationId totalMuts = grg->numMutations();
    for (MutationId chunkStart = 0; chunkStart < totalMuts; chunkStart += chunkSize) {
        batch.clear();
        const MutationId chunkEnd = std::min<MutationId>(chunkStart + chunkSize, totalMuts);

        for (MutationId mutId = chunkStart; mutId < chunkEnd; ++mutId) {
            const Mutation& mutation = grg->getMutationById(mutId);
            const uint64_t pos = mutation.getPosition();
            const std::string& ref = mutation.getRefAllele();
            const std::string& alt = mutation.getAllele();

            const size_t maxAlleleLen = std::max(ref.size(), alt.size());
            if (maxAlleleLen == 0) {
                continue;
            }
            if (positionsAreOneBased && pos == 0) {
                continue;
            }
            const size_t start = positionsAreOneBased ? static_cast<size_t>(pos - 1) : static_cast<size_t>(pos);
            if (start >= fastaLen) {
                continue;
            }
            if (maxAlleleLen > fastaLen - start) {
                continue;
            }

            const std::string ancestral = fastaSeq.substr(start, maxAlleleLen);
            if (ancestral.empty()) {
                continue;
            }

            const char base = ancestral[0];
            if (base == 'N' || base == '.' || base == '-') {
                continue;
            }

            if (!ref.empty() && equalsIgnoreCase(ref, ancestral.substr(0, ref.size()))) {
                batch.emplace_back(mutId, ref);
                continue;
            }
            if (!alt.empty() && equalsIgnoreCase(alt, ancestral.substr(0, alt.size()))) {
                batch.emplace_back(mutId, alt);
                continue;
            }

            if (dropIfNoMatch) {
                batch.emplace_back(mutId, "N");
            }
        }

        if (!batch.empty()) {
            polarizeMutationsHelper(grg, mutLookup, batch, stats);
        }
    }

    grg->sortMutations();

    return stats;
}

} // namespace grgl
