#include "grgl/polarization.h"

#include "grgl/map_mutations.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
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
                               std::vector<Mutation>& remapMutations,
                               std::vector<NodeIDList>& remapSamples,
                               size_t flushThreshold,
                               bool force) {
    if (remapMutations.empty()) {
        return;
    }
    if (force || remapMutations.size() >= flushThreshold) {
        mapMutations(grg, remapMutations, remapSamples, false, remapMutations.size());
        remapMutations.clear();
        remapSamples.clear();
    }
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

    // Build a direct lookup from mutation ID to (nodeId, missingNodeId) to avoid
    // repeatedly scanning the ordered vector (which was O(n^2)).
    std::vector<std::pair<NodeID, NodeID>> mutLookup(grg->numMutations(),
                                                     {INVALID_NODE_ID, INVALID_NODE_ID});
    for (const auto& entry : mutNodeList) {
        const MutationId mid = std::get<0>(entry);
        if (mid < mutLookup.size()) {
            mutLookup[mid] = {std::get<1>(entry), std::get<2>(entry)};
        }
    }

    std::vector<Mutation> remapMutations;
    std::vector<NodeIDList> remapSamples;

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

            if (!missingSamples.empty()) {
                remapMutations.emplace_back(mutation.getPosition(),
                                            Mutation::ALLELE_MISSING,
                                            mutation.getRefAllele(),
                                            mutation.getTime());
                remapSamples.emplace_back(std::move(missingSamples));
            }

            remapMutations.emplace_back(mutation.getPosition(), ref, ancestralAllele, mutation.getTime());
            remapSamples.emplace_back(std::move(flippedCarriers));

            loadRemaps(grg, remapMutations, remapSamples, remapFlushThreshold, false);
            continue;
        }

        // inconsistent ancestral allele, drop.
        stats.inconsistent++;
        grg->removeMutation(mutId, nodeId);
    }

    loadRemaps(grg, remapMutations, remapSamples, remapFlushThreshold, true);

    return results;
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

    std::vector<std::pair<MutationId, std::string>> batch;
    batch.reserve(grg->numMutations());

    for (MutationId mutId = 0; mutId < grg->numMutations(); ++mutId) {
        const Mutation& mutation = grg->getMutationById(mutId);
        const uint64_t pos = mutation.getPosition();
        const std::string& ref = mutation.getRefAllele();
        const std::string& alt = mutation.getAllele();

        const size_t maxAlleleLen = std::max(ref.size(), alt.size());
        size_t start = positionsAreOneBased ? static_cast<size_t>(pos - 1) : static_cast<size_t>(pos);
        if (positionsAreOneBased && pos == 0) {
            continue;
        }
        if (start + maxAlleleLen > fastaLen) {
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
        stats.reset();
        polarizeMutations(grg, batch, stats);
    }

    return stats;
}

} // namespace grgl

