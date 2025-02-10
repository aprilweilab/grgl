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
#ifndef GRG_POOLED_JOBS_H
#define GRG_POOLED_JOBS_H

#include <iostream>
#include <list>
#include <mutex>
#include <thread>
#include <vector>

namespace grgl {

/**
 * Conceptually similar to the Python multiprocessing.Pool, except using threads instead of processes.
 * All of the work must be added to the queue (via addWork) _prior_ to calling doAllWork. This is essentially
 * just a way to batch a bunch of known work.
 *
 * In order to use this, inherit from it and define the processItem method which do a single item's worth of
 * work. All the synchronization is handled by PooledJobs, _unless_ processing an item requires the use of
 * shared resources that are not encapsulated by the item type T.
 */
template <typename T> class PooledJobs {
public:
    PooledJobs() = default;
    virtual ~PooledJobs() = default;

    PooledJobs(PooledJobs& rhs) = delete;
    PooledJobs(PooledJobs&& rhs) = delete;
    PooledJobs& operator=(PooledJobs& rhs) = delete;
    PooledJobs& operator=(PooledJobs&& rhs) = delete;

    void addWork(T item) {
        std::lock_guard<std::mutex> lock(m_queueMutex);
        m_queue.push_back(std::move(item));
    }

    void doAllWork(size_t jobs) {
        if (jobs == 1) {
            while (!m_queue.empty()) {
                T nextItem = std::move(m_queue.front());
                m_queue.pop_front();
                processItem(std::move(nextItem));
            }
        } else {
            std::vector<std::thread> threadPool;
            for (size_t i = 0; i < jobs; i++) {
                threadPool.emplace_back(&PooledJobs::workerMethod, this);
            }
            for (size_t i = 0; i < jobs; i++) {
                threadPool[i].join();
            }
        }
    }

protected:
    // The worker (a single thread) loops until the queue is empty.
    void workerMethod() {
        while (true) {
            T nextItem; // T needs a default constructor.
            {
                std::lock_guard<std::mutex> lock(m_queueMutex);
                if (m_queue.empty()) {
                    break;
                }
                nextItem = std::move(m_queue.front());
                m_queue.pop_front();
            }
            processItem(std::move(nextItem));
        }
    }

    virtual void processItem(T item) = 0;

private:
    std::mutex m_queueMutex;
    std::list<T> m_queue;
};

} // namespace grgl

#endif /* GRG_POOLED_JOBS_H */
