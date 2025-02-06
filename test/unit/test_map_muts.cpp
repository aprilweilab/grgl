#include <gtest/gtest.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <limits>
#include <streambuf>
#include <unistd.h>
#include <vector>

#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/mut_iterator.h"
#include "grgl/map_mutations.h"
#include "common_grgs.h"
#include "grgl/node_data.h"
#include "grgl/serialize.h"
#include "grg_helpers.h"

using namespace grgl;

class TestMutationIterator : public MutationIterator {
public:
    explicit TestMutationIterator(std::vector<Mutation> muts, std::vector<NodeIDList> samples)
            : MutationIterator({0, 1000}, false, false) {
        for (size_t i = 0; i < muts.size(); i++) {
            this->m_alreadyLoaded.push_back({muts.at(i), samples.at(i)});
        }
        std::reverse(this->m_alreadyLoaded.begin(), this->m_alreadyLoaded.end());
    }

    void getMetadata(size_t& ploidy, size_t& numIndividuals, bool& isPhased) override {
        ploidy = 2;
        numIndividuals = 4;
        isPhased = true;
    }

    size_t countMutations() const override {
        return m_alreadyLoaded.size();
    }

    std::vector<std::string> getIndividualIds() override {
        return {"1", "2", "3", "4"};
    }

protected:
    void buffer_next(size_t& totalSamples) override {
    }

    void reset_specific() override {
    }

private:
    size_t m_currentVariant{};
};

// This tests a few different things:
// 1. That MapMutations is adding the mutations once and only once
// 2. That newly added nodes will be re-used
// 3. That coalescences are computed, and that they only represent _newly_ coalesced indivs
TEST(MapMutations, Basic) {
    MutableGRGPtr grg = sample8Grg();
    TestMutationIterator testIt(
        {
            {1, "A", "G", 0.0},
            {2, "G", "C", 0.0},
            {3, "G", "A", 0.0},
            {4, "A", "T", 0.0},
        },
        {
            {0, 1, 3, 7},
            {1, 2},
            {6, 7},
            {4, 5, 6, 7}, // re-uses previously added (6, 7)
        }
    );
    const NodeIDSizeT nodesBefore = grg->numNodes();
    for (NodeID nodeId = 0; nodeId < nodesBefore; nodeId++) {
        if (grg->isSample(nodeId)) {
            ASSERT_EQ(grg->getNumIndividualCoals(nodeId), 0);
        } else {
            ASSERT_EQ(grg->getNumIndividualCoals(nodeId), COAL_COUNT_NOT_SET);
        }
    }

    MutationMappingStats stats = mapMutations(grg, testIt);
    const NodeIDSizeT nodesAfter = grg->numNodes();
    std::cout << stats.totalMutations << "\n";
    std::cout << "Added " << nodesAfter << " nodes\n";
    for (auto mutAndNode : grg->getMutationsToNodeOrdered()) {
        auto mut = grg->getMutationById(mutAndNode.first);
        switch (mut.getPosition()) {
            case 1:
                ASSERT_EQ(grg->getNumIndividualCoals(mutAndNode.second), 0);
                break;
            case 2:
                ASSERT_EQ(grg->getNumIndividualCoals(mutAndNode.second), 0);
                break;
            case 3:
            	// 6, 7 coalesced here
                ASSERT_EQ(grg->getNumIndividualCoals(mutAndNode.second), 1);
                break;
            case 4:
            	// 4, 5 coalesced here (it should reuse the previous (6, 7) node, which
                // will have it's own coal count for 4,5)
                ASSERT_EQ(grg->getNumIndividualCoals(mutAndNode.second), 1);
                break;
            default:
                ASSERT_FALSE(true);
        }
    }
    ASSERT_EQ(grg->numMutations(), 4);

	// We have nodes from the unmapped GRG that we know will coalesce
    ASSERT_EQ(grg->getNumIndividualCoals(8), 1);
    size_t totalCoals = 0;
    for (NodeID nodeId = 0; nodeId < nodesAfter; nodeId++) {
        const NodeIDSizeT coals = grg->getNumIndividualCoals(nodeId);
        if (coals != COAL_COUNT_NOT_SET) {
            totalCoals += grg->getNumIndividualCoals(nodeId);
        }
    }
    ASSERT_EQ(totalCoals, 3);
}