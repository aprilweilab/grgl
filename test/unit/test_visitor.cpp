#include <gtest/gtest.h>

#include "grgl/grg.h"
#include "grgl/visitor.h"

#include "common_grgs.h"

using namespace grgl;

class TopoCheckVisitor : public GRGVisitor {
public:
    TopoCheckVisitor()
        : m_visits(0)
        , m_invalidVisits(0) {
    }

    bool visit(const GRGPtr& grg,
               const NodeID nodeId,
               const TraversalDirection direction,
               const DfsPass dfsPass = DfsPass::DFS_PASS_NONE) override {
        if (m_visited.empty()) {
            m_visited.resize(grg->numNodes());
        }
        if (dfsPass != DfsPass::DFS_PASS_THERE) {
            auto predecessors = (direction == DIRECTION_UP) ? grg->getDownEdges(nodeId)
                                                            : grg->getUpEdges(nodeId);
            for (const auto& predId : predecessors) {
                //ASSERT_LT(pred->getId(), m_visited.size());
                if (!m_visited[predId]) {
                    m_invalidVisits++;
                }
            }
            m_visited[nodeId] = true;
            m_visits++;
        }
        return true;
    }

    std::vector<bool> m_visited;
    size_t m_visits;
    size_t m_invalidVisits;
};

TEST(VisitorsTest, TopologicalVisitor) {
    GRGPtr grg = depth3BinTree();

    TopoCheckVisitor visitorUp;
    grg->visitTopo(visitorUp, TraversalDirection::DIRECTION_UP, grg->getSampleNodes());
    // The number of times the topological order property was broken should be 0.
    ASSERT_EQ(visitorUp.m_invalidVisits, 0);
    // Visit each node once and only once.
    ASSERT_EQ(grg->numNodes(), visitorUp.m_visits);

    TopoCheckVisitor visitorDown;
    grg->visitTopo(visitorDown, TraversalDirection::DIRECTION_DOWN, grg->getRootNodes());
    // The number of times the topological order property was broken should be 0.
    ASSERT_EQ(visitorDown.m_invalidVisits, 0);
    // Visit each node once and only once.
    ASSERT_EQ(grg->numNodes(), visitorDown.m_visits);
}