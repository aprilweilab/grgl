#ifndef GRGL_CUDA_MATMUL_H
#define GRGL_CUDA_MATMUL_H

#ifdef GRGL_CUDA_ENABLED

#include "grgl/grg.h"
#include "grgl/serialize.h"
#include "grg_helpers.h"
#include "grgl/common.h"
#include "grgl/csr_storage.h"
#include "grgl/grgnode.h"
#include "grgl/mutation.h"
#include "grgl/node_data.h"
#include "util.h"
#include <cuda_runtime.h>

namespace grgl {

//typedef uint32_t index_t;
//typedef double data_t;

template<class index_t, class data_t>
struct CudaCSR {
    index_t* row_offsets;    // Device pointer to row offsets (size: num_rows + 1)
    index_t* col_indices;    // Device pointer to column indices (size: nnz)
    index_t* host_height_cutoffs;
    size_t num_rows;         // Number of rows
    size_t num_cols;         // Number of columns  
    size_t nnz;              // Number of non-zero elements
    size_t max_height;
    
    // Constructor
    CudaCSR() : row_offsets(nullptr), col_indices(nullptr), 
        host_height_cutoffs(nullptr),
        num_rows(0), num_cols(0), nnz(0), max_height(0)
    {}
    
    // Destructor
    ~CudaCSR() {
        free();
    }
    
    // Allocate GPU memory
    void allocate(size_t rows, size_t cols, size_t non_zeros, size_t max_h) {
        free(); // Free any existing memory
        num_rows = rows;
        num_cols = cols;
        nnz = non_zeros;
        max_height = max_h;
        
        cudaMalloc(&row_offsets, (num_rows + 1) * sizeof(index_t));
        cudaMalloc(&col_indices, nnz * sizeof(index_t));
        host_height_cutoffs = new index_t[max_height + 1];
    }
    
    // Free GPU memory
    void free() {
        if (row_offsets) {
            cudaFree(row_offsets);
            row_offsets = nullptr;
        }
        if (col_indices) {
            cudaFree(col_indices);
            col_indices = nullptr;
        }
        if (host_height_cutoffs) {
            delete[] host_height_cutoffs;
            host_height_cutoffs = nullptr;
        }
        num_rows = num_cols = nnz = 0;
    }
    
    // Copy data to GPU
    void copyToDevice(const index_t* host_row_offsets, const index_t* host_col_indices) {
        if (row_offsets && col_indices) {
            cudaMemcpy(row_offsets, host_row_offsets, (num_rows + 1) * sizeof(index_t), cudaMemcpyHostToDevice);
            cudaMemcpy(col_indices, host_col_indices, nnz * sizeof(index_t), cudaMemcpyHostToDevice);
        }
    }

    // Debugging: Print CSR structure
    void print() const {
        std::cout << "CudaCSR: " << std::endl;
        std::cout << "  num_rows: " << num_rows << std::endl;
        std::cout << "  num_cols: " << num_cols << std::endl;
        std::cout << "  nnz: " << nnz << std::endl;
        std::cout << "  max_height: " << max_height << std::endl;

        // now copy core data back to cpu and print
        std::vector<index_t> h_row_offsets(num_rows + 1);
        std::vector<index_t> h_col_indices(nnz);
        cudaMemcpy(h_row_offsets.data(), row_offsets, (num_rows + 1) * sizeof(index_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_col_indices.data(), col_indices, nnz * sizeof(index_t), cudaMemcpyDeviceToHost);
        std::cout << "  row_offsets: ";
        for (size_t i = 0; i < num_rows + 1; i++) {
            std::cout << h_row_offsets[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "  col_indices: ";
        for (size_t i = 0; i < nnz; i++) {
            std::cout << h_col_indices[i] << " ";
        }
        std::cout << std::endl;
    }
};

// This visitor is tested ONLY on CSRGRG (immutable, ordered GRG)
// in the first visitor path, we get the height of each node, and thus the new id. new id is of (height, gid)
// we construct both old_id -> new_id map and a new_id -> old_id map
// then, in the epilogue function, we can use the information to have a new CSR
template<class index_t, class data_t>
class GPUCsrVisitor : public GRGVisitor {
public:
    explicit GPUCsrVisitor(size_t num_nodes, size_t num_elements) {
        m_maxHeight = 0;
        m_num_nodes = num_nodes;
        m_num_elements = num_elements;
        m_newId.resize(num_nodes, std::make_pair(0,0));
    }

    bool visit(const grgl::GRGPtr& grg,
               const grgl::NodeID nodeId,
               const grgl::TraversalDirection direction,
               const grgl::DfsPass dfsPass) override {
        
        if (direction != DIRECTION_DOWN) {
            release_assert(false);
            return true;
        }

        if (dfsPass == DfsPass::DFS_PASS_THERE) {
            return true;    // do nothing on the forward pass
        }

        index_t height = 0;

        if (nodeId < grg->numSamples()) {
            height = 0;
        } else {
            auto downEdges = grg->getDownEdges(nodeId);
            for (const auto& child : downEdges) {
                release_assert(child < grg->numNodes());
                auto child_height = m_newId[child].first;
                height = std::max(height, child_height + 1);
            }
        }

        if (height >= m_maxHeight) {
            m_maxHeight = height + 1;
            // m_idCounters.resize(m_maxHeight, 0);
            m_oldId.resize(m_maxHeight);
        }

        // auto id_in_group = m_idCounters[height]++;
        auto id_in_group = m_oldId[height].size();
        m_oldId[height].push_back(nodeId);
        m_newId[nodeId] = std::make_pair(height, id_in_group);
        return true;
    }

    // After visiting all nodes, we can construct the CSR representation
    // This function should be called only once after all visits are done
    void constructCSR(
        GRG* grg,
        CudaCSR<index_t, data_t>& gpu_csr
    ) {
        gpu_csr.allocate(m_num_nodes, m_num_nodes, m_num_elements, m_maxHeight);

        gpu_csr.host_height_cutoffs[0] = 0;
        for (int i = 1; i <= m_maxHeight; i++) {
            gpu_csr.host_height_cutoffs[i] = gpu_csr.host_height_cutoffs[i - 1] + m_oldId[i - 1].size();
        }

        std::vector<index_t> host_row_offsets;
        host_row_offsets.resize(m_num_nodes + 1, 0);
        std::vector<index_t> host_col_indices;

        m_childCounts.resize(m_maxHeight, 0);

        index_t current_col_index = 0;
        for (int h = 0; h < m_maxHeight; h++) {
            for (size_t i = 0; i < m_oldId[h].size(); i++) {
                NodeID old_id = m_oldId[h][i];
                index_t new_id = gpu_csr.host_height_cutoffs[h] + i;
                host_row_offsets[new_id] = current_col_index;
                current_col_index += grg->numDownEdges(old_id);
                
                // for each down edge, we need to find its new id
                auto downEdges = grg->getDownEdges(old_id);
                m_childCounts[h] += downEdges.size();
                for (const auto& child : downEdges) {
                    release_assert(child < grg->numNodes());
                    auto child_new_id = m_newId[child];
                    index_t child_new_flat_id = gpu_csr.host_height_cutoffs[child_new_id.first] + child_new_id.second;
                    host_col_indices.push_back(child_new_flat_id);
                }
            }
        }

        host_row_offsets[m_num_nodes] = current_col_index;
        release_assert(current_col_index == m_num_elements);

        gpu_csr.copyToDevice(host_row_offsets.data(), host_col_indices.data());
    }

    // Flattened mapping from new_id to old_id
    std::vector<index_t> getOldId() const {
        std::vector<index_t> flat_old_id;
        for (const auto& vec : m_oldId) {
            flat_old_id.insert(flat_old_id.end(), vec.begin(), vec.end());
        }
        return flat_old_id;
    }

    // Flattened mapping from old_id to new_id
    std::vector<index_t> getNewId() const {
        std::vector<index_t> flat_new_id;
        flat_new_id.resize(m_num_nodes);
        auto base_id = 0;
        for (size_t height = 0; height < m_oldId.size(); height++) {
            for (size_t i = 0; i < m_oldId[height].size(); i++) {
                NodeID old_id = m_oldId[height][i];
                index_t new_id = base_id + i;
                flat_new_id[old_id] = new_id;
            }
            base_id += m_oldId[height].size();
        }
        return flat_new_id;
    }

    void printChildCounts() const {
        std::cout << "Child counts per height:" << std::endl;
        for (size_t h = 0; h < m_childCounts.size(); h++) {
            std::cout << "Height " << h << ": " << m_childCounts[h] << std::endl;
        }
    }

    void printNodeCounts() const {
        std::cout << "Node counts per height:" << std::endl;
        for (size_t h = 0; h < m_oldId.size(); h++) {
            std::cout << "Height " << h << ": " << m_oldId[h].size() << std::endl;
        }
    }

private:
    index_t m_maxHeight;
    size_t m_num_nodes;
    size_t m_num_elements;

    std::vector<size_t> m_childCounts; // sum of child counts for nodes of the same height
    std::vector<std::vector<index_t>> m_oldId;  // new_id -> old_id mapping
    std::vector<std::pair<index_t, index_t>> m_newId; // old_id -> new_id mapping

};

// CUDA-specific utility functions
bool checkCudaDevice();

// void cudaMatMulInit();
// void cudaMatMulCleanup();

} // namespace grgl

#endif // GRGL_CUDA_ENABLED
#endif // GRGL_CUDA_MATMUL_H