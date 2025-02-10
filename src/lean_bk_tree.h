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
 * You should have received a copy of the GNU General Public License
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef GRGL_LEAN_BK_TREE_H
#define GRGL_LEAN_BK_TREE_H

#include <cstdlib>
#include <functional>
#include <limits>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

namespace grgl {

template <typename ElementType> class LeanBKTree;

template <typename ElementType> class LeanBKTreeNode {
public:
    explicit LeanBKTreeNode(const ElementType& element)
        : m_elements({element}),
          m_isDeleted(false) {}

    void addElement(ElementType element) {
        if (m_isDeleted) {
            m_elements = {std::move(element)};
            m_isDeleted = false;
        } else {
            m_elements.push_back(std::move(element));
        }
    }

    template <typename Container> void moveElements(Container& result) {
        for (size_t i = 0; i < m_elements.size(); i++) {
            result.push_back(m_elements[i]);
        }
        m_isDeleted = true;
    }

private:
    using ChildPair = std::pair<size_t, std::shared_ptr<LeanBKTreeNode<ElementType>>>;

    // The elements (user data) associated with the above vector.
    std::vector<ElementType> m_elements;
    // The child edges. Likely sparse, depending on how many discrete distances are possible.
    std::vector<ChildPair> m_children;

    bool m_isDeleted;

    friend class LeanBKTree<ElementType>;
};

template <typename ElementType> class LeanBKTree {
public:
    using NodePointer = std::shared_ptr<LeanBKTreeNode<ElementType>>;

    explicit LeanBKTree(std::function<size_t(const ElementType&, const ElementType&)> distFunc)
        : m_distFunc(distFunc) {}

    NodePointer insert(const ElementType& element, size_t& comparisons) {
        if (!m_rootNode) {
            m_rootNode = std::make_shared<LeanBKTreeNode<ElementType>>(element);
            return m_rootNode;
        }
        NodePointer currentNode = m_rootNode;
        while (currentNode) {
            comparisons++;
            const size_t dist = m_distFunc(currentNode->m_elements[0], element);
            if (dist == 0) {
                currentNode->addElement(element);
                return currentNode;
            }
            NodePointer nextNodePtr;
            for (auto& childPair : currentNode->m_children) {
                if (childPair.first == dist) {
                    nextNodePtr = childPair.second;
                }
            }
            if (!nextNodePtr) {
                auto newNode = std::make_shared<LeanBKTreeNode<ElementType>>(element);
                currentNode->m_children.emplace_back(dist, newNode);
                return newNode;
            }
            currentNode = nextNodePtr;
        }
        abort();
    }

    std::vector<NodePointer>
    lookup(const ElementType& queryElement, size_t& nearestDistance, size_t& comparisons, bool collectAll = true) {
        size_t distBest = std::numeric_limits<size_t>::max();
        if (!m_rootNode) {
            nearestDistance = distBest;
            return {};
        }
        std::list<NodePointer> workingList = {m_rootNode};
        std::vector<NodePointer> results;
        while (!workingList.empty()) {
            auto node = workingList.back();
            workingList.pop_back();
            const auto dist = m_distFunc(node->m_elements[0], queryElement);
            comparisons++;
            // We do not return nodes that have no elements.
            bool ignoreThisNode = node->m_isDeleted;
            if (!ignoreThisNode && node->m_elements.size() == 1 && node->m_elements[0] == queryElement) {
                ignoreThisNode = true;
                node->m_isDeleted = true;
            }
            if (!ignoreThisNode) {
                if (dist < distBest) {
                    results = {node};
                    distBest = dist;
                } else if (collectAll && dist == distBest) {
                    results.push_back(node);
                }
            }
            for (auto& pair : node->m_children) {
                const size_t nextDist = pair.first;
                const size_t bound = std::abs((ssize_t)nextDist - (ssize_t)dist);
                if (bound < distBest || (bound == distBest && collectAll)) {
                    workingList.push_back(pair.second);
                }
            }
        }
        nearestDistance = distBest;
        return results;
    }

private:
    NodePointer m_rootNode{};
    std::function<size_t(const ElementType&, const ElementType&)> m_distFunc;
};

} // namespace grgl

#endif /* GRGL_LEAN_BK_TREE_H */
