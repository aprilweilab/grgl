#ifndef GRGL_GPU_GRG_H
#define GRGL_GPU_GRG_H

#ifdef GRGL_CUDA_ENABLED

#define CHECK_CUDA_LAST_ERROR() {                                          \
    cudaError_t err = cudaGetLastError();                                 \
    if (err != cudaSuccess) {                                             \
        std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__      \
                  << " code=" << static_cast<int>(err)                    \
                  << " \"" << cudaGetErrorString(err) << "\"" << std::endl; \
        throw std::runtime_error("CUDA kernel launch failed");             \
    }                                                                     \
}

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <chrono>
#include <cuda_runtime.h>
// #include "grgl/cuda_matmul.h"
#include "grgl/common.h"
#include "grgl/grg.h"

namespace grgl {

template<class index_t>
class GPUGRG {
public:

    // these are resident GPU pointers
    index_t* row_offsets;    // Device pointer to row offsets (size: num_rows + 1)
    index_t* old_to_new_mapping;
    index_t* new_to_old_mapping;
    index_t* mutation_and_new_mapping; // Mapping between mutation id and new node id (in pairs)

    // Since col indices may be large, this data may not be resident on GPU
    index_t* col_indices;    // Device pointer to column indices (size: nnz)

    // This is the corresponding host backup data
    std::vector<index_t> host_col_indices;

    std::vector<index_t> host_height_cutoffs;
    std::vector<index_t> host_heavy_cutoffs;
    std::vector<double> host_avg_child_counts;

    // std::vector<index_t> host_old_to_new_mapping;  // old_id -> new_id
    // std::vector<index_t> host_new_to_old_mapping;  // new_id -> old_id

    size_t num_rows;         // Number of rows, i.e. number of nodes
    size_t num_samples;      // Number of samples, i.e. number of leaf nodes
    size_t num_mutations;    // Number of mutations. may NOT correspond to number of non-leaf nodes
    size_t nnz;              // Number of non-zero elements, i.e. number of edges
    size_t max_height;

    cudaStream_t* work_stream_ptr;  // CUDA stream for asynchronous operations

public:
    // Constructor
    GPUGRG() : row_offsets(nullptr), col_indices(nullptr),  old_to_new_mapping(nullptr), new_to_old_mapping(nullptr), mutation_and_new_mapping(nullptr),
        num_rows(0), num_samples(0), nnz(0), max_height(0), num_mutations(0),
        work_stream_ptr(nullptr)
    {}
    
    // Destructor
    ~GPUGRG() {
        free();
    }

    // return the size of GPU memory needed by col indices in bytes
    size_t queryColIndSize() {
        return nnz * sizeof(index_t);
    }

    // This function sets existing GPU memory pointers for col indices
    // will not allocate new memory
    void setColIndBuffer(void* ptr, size_t size) {
        if (size < queryColIndSize()) {
            throw ApiMisuseFailure("Provided buffer size is smaller than required col indices size");
        }
        col_indices = static_cast<index_t*>(ptr);
    }

    // This function copies col indices from Host to GPU
    // Asynchronous copy is supported
    // The caller must ensure the memory space is not used by other operations.
    void copyColIndFromHost() {
        if (col_indices) {
            cudaMemcpyAsync(col_indices, host_col_indices.data(), nnz * sizeof(index_t), cudaMemcpyHostToDevice, *work_stream_ptr);
        }
    }

    // This function will block the calling thread 
    // until all preceding operations related to this GRG is done
    void wait() {
        if (work_stream_ptr) {
            cudaStreamSynchronize(*work_stream_ptr);
        } else {
            throw ApiMisuseFailure("No CUDA stream set for GPUGRG wait operation");
        }
    }

    // This function removes the GPU memory buffer
    // Will not free the memory, just removes the pointers
    void removeColIndBuffer() {
        col_indices = nullptr;
    }

    // This function frees the GPU memory buffer for col indices
    void freeColInd() {
        cudaFree(col_indices);
        col_indices = nullptr;
    }
    
    // Init vars and allocate GPU memory for all structures
    // This function is blocking
    // If no stream_ptr is provided, a new stream will be created
    void init(size_t rows, size_t samples, size_t non_zeros, size_t mutations,size_t max_h, bool allocate_col_indices = true, cudaStream_t* stream_ptr = nullptr) {
        free(); // Free any existing memory
        CHECK_CUDA_LAST_ERROR();
        num_rows = rows;
        num_samples = samples;
        nnz = non_zeros;
        num_mutations = mutations;
        max_height = max_h;
        
        cudaMalloc(&row_offsets, (num_rows + 1) * sizeof(index_t));
        cudaMalloc(&old_to_new_mapping, num_rows * sizeof(index_t));
        cudaMalloc(&new_to_old_mapping, num_rows * sizeof(index_t));
        cudaMalloc(&mutation_and_new_mapping, num_mutations * 2 * sizeof(index_t)); // Allocate for pairs
        CHECK_CUDA_LAST_ERROR();
        if (allocate_col_indices) {
            cudaMalloc(&col_indices, nnz * sizeof(index_t));
        } else {
            col_indices = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();

        host_col_indices.resize(nnz);
        host_height_cutoffs.resize(max_height + 1);
        host_heavy_cutoffs.resize(max_height + 1);
        host_avg_child_counts.resize(max_height + 1);

        if (stream_ptr) {
            work_stream_ptr = stream_ptr;
        } else {
            work_stream_ptr = new cudaStream_t();
            cudaStreamCreate(work_stream_ptr);
        }
        CHECK_CUDA_LAST_ERROR();
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
        if (old_to_new_mapping) {
            cudaFree(old_to_new_mapping);
            old_to_new_mapping = nullptr;
        }
        if (new_to_old_mapping) {
            cudaFree(new_to_old_mapping);
            new_to_old_mapping = nullptr;
        }
        if (mutation_and_new_mapping) {
            cudaFree(mutation_and_new_mapping);
            mutation_and_new_mapping = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();

        host_col_indices.clear();
        host_height_cutoffs.clear();
        host_heavy_cutoffs.clear();
        host_avg_child_counts.clear();

        num_rows = num_samples = nnz = 0;
        num_mutations = 0;
        max_height = 0;

        if (work_stream_ptr) {
            cudaStreamDestroy(*work_stream_ptr);
            delete work_stream_ptr;
            work_stream_ptr = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();
    }
    
    // Copy data to GPU
    // This function is blocking
    void copyToDevice(
            const index_t* row_offsets, 
            const index_t* old_to_new_mapping,
            const index_t* new_to_old_mapping,
            const index_t* mutation_and_new_mapping,
            const index_t* height_cutoffs,
            const index_t* heavy_cutoffs,
            const double* avg_child_counts,
            const index_t* col_indices,
            bool copy_col_indices_to_device = true
    ) {
        
        if (!this->row_offsets || !this->old_to_new_mapping || !this->new_to_old_mapping) {
            throw ApiMisuseFailure("Device memory not allocated yet!");
        }
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->row_offsets, row_offsets, (num_rows + 1) * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->old_to_new_mapping, old_to_new_mapping, num_rows * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->new_to_old_mapping, new_to_old_mapping, num_rows * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->mutation_and_new_mapping, mutation_and_new_mapping, num_mutations * 2 * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();

        // cpu memcpy just use normal memcpys
        memcpy(this->host_height_cutoffs.data(), height_cutoffs, (max_height + 1) * sizeof(index_t));
        memcpy(this->host_heavy_cutoffs.data(), heavy_cutoffs, (max_height + 1) * sizeof(index_t));
        memcpy(this->host_avg_child_counts.data(), avg_child_counts, (max_height + 1) * sizeof(double));
        memcpy(this->host_col_indices.data(), col_indices, nnz * sizeof(index_t));
        CHECK_CUDA_LAST_ERROR();

        if (copy_col_indices_to_device) {
            if (!this->col_indices) {
                throw ApiMisuseFailure("Column indices GPU buffer is not allocated");
            }
            cudaMemcpy(this->col_indices, col_indices, nnz * sizeof(index_t), cudaMemcpyHostToDevice);
        }
        CHECK_CUDA_LAST_ERROR();
    }

    // Debugging: Print CSR structure
    void print() const {
        std::cout << "GPUGRG: " << std::endl;
        std::cout << "  num_rows: " << num_rows << std::endl;
        std::cout << "  nnz: " << nnz << std::endl;
        std::cout << "  max_height: " << max_height << std::endl;
        std::cout << "  num_samples: " << num_samples << std::endl;
        std::cout << "  num_mutations: " << num_mutations << std::endl;

        // now copy core data back to cpu and print
        /*
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
        */
    }
    
    /**
     * Blocking Matrix multiplication on GPU using the GPUGRG structure.
     * 
     * @param inputMatrix The input matrix in row-major order
     * @param numRows Number of rows in the input matrix
     * @param direction Traversal direction (UP or DOWN)
     * @return The resulting matrix as a vector
     */
    template <typename T>
    std::vector<T> matMulBlocking(
        const std::vector<T>& inputMatrix, 
        const size_t numRows, 
        TraversalDirection direction) {
        CHECK_CUDA_LAST_ERROR();
        if (numRows == 0 || (inputMatrix.size() % numRows != 0)) {
            throw ApiMisuseFailure("inputMatrix must be divisible by numRows");
        }
        const size_t numCols = inputMatrix.size() / numRows;
        const size_t outSize =
            (numRows * ((direction == TraversalDirection::DIRECTION_DOWN) ? this->num_samples : this->num_mutations));
        
        T* d_inputMatrix;
        T* d_outputMatrix;
        cudaMalloc(&d_inputMatrix, inputMatrix.size() * sizeof(T));
        cudaMalloc(&d_outputMatrix, outSize * sizeof(T));
        cudaMemcpy(d_inputMatrix, inputMatrix.data(), inputMatrix.size() * sizeof(T), cudaMemcpyHostToDevice);

        CHECK_CUDA_LAST_ERROR();

        T* d_buffer;
        size_t buffer_size_byte = numRows * this->num_rows * sizeof(T);
        cudaMalloc(&d_buffer, buffer_size_byte);

        CHECK_CUDA_LAST_ERROR();

        this->matrixMultiplication<T>(
            d_inputMatrix, 
            numCols, 
            numRows, 
            direction, 
            d_outputMatrix, 
            outSize,
            false,
            d_buffer,
            buffer_size_byte
        );
        cudaDeviceSynchronize();
        CHECK_CUDA_LAST_ERROR();

        std::vector<T> result(outSize);
        cudaMemcpy(result.data(), d_outputMatrix, outSize * sizeof(T), cudaMemcpyDeviceToHost);
        cudaFree(d_inputMatrix);
        cudaFree(d_outputMatrix);
        cudaFree(d_buffer);
        CHECK_CUDA_LAST_ERROR();

        return std::move(result);
    }
    
    /*
     * Raw matrix multiplication implementation.
     * The output matrix is the full matrix.
     * All input and output pointers must be device pointers.
     * The user should guarantee that col_indices pointer is valid, or a copy is in progress with the provided function.
     */
    template <class data_t>
    void matrixMultiplication(const data_t* inputMatrix,
                            size_t inputCols,
                            size_t inputRows,
                            TraversalDirection direction,
                            data_t* outputMatrix,
                            size_t outputSize,
                            bool emitAllNodes,
                            data_t* buffer,
                            size_t buffer_size
                            );
};


namespace {
template<class index_t>
struct nodes {
        index_t id;
        index_t child_count;
    };

template<class index_t>
bool cmp(const nodes<index_t>& a, const nodes<index_t>& b) {
    return a.child_count > b.child_count;
}
}

// This visitor is tested ONLY on CSRGRG (immutable, ordered GRG)
// in the first visitor path, we get the height of each node, and thus the new id. new id is of (height, gid)
// we construct both old_id -> new_id map and a new_id -> old_id map
// then, in the epilogue function, we can use the information to have a new GPUGRG
template<class index_t>
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

    // this function rearranges nodes at each height by their child counts
    void rearrange(GRG* grg) {
        for (int h = 1; h < m_maxHeight; h++) {
            std::cout << "Rearranging height " << h << " with " << m_oldId[h].size() << " nodes." << std::endl;
            std::vector<nodes<index_t>> node_list;
            for (size_t i = 0; i < m_oldId[h].size(); i++) {
                NodeID old_id = m_oldId[h][i];
                index_t child_count = grg->numDownEdges(old_id);
                node_list.push_back({old_id, child_count});
            }
            std::sort(node_list.begin(), node_list.end(), cmp<index_t>);
            for (size_t i = 0; i < node_list.size(); i++) {
                m_oldId[h][i] = node_list[i].id;
                NodeID old_id = node_list[i].id;
                m_newId[old_id] = std::make_pair(h, i);
            }
        }
    }

    // After visiting all nodes, we can construct the CSR representation
    // This function should be called only once after all visits are done
    void constructGPUGRG(
        GRG* grg,
        GPUGRG<index_t>& gpu_grg
    ) {

        std::vector<index_t> host_height_cutoffs;
        host_height_cutoffs.resize(m_maxHeight + 1, 0);
        for (int i = 1; i <= m_maxHeight; i++) {
            host_height_cutoffs[i] = host_height_cutoffs[i - 1] + m_oldId[i - 1].size();
        }

        std::vector<index_t> host_heavy_cutoffs;
        host_heavy_cutoffs.resize(m_maxHeight + 1, 0);

        std::vector<index_t> host_row_offsets;
        host_row_offsets.resize(m_num_nodes + 1, 0);

        std::vector<double> host_avg_child_counts;
        host_avg_child_counts.resize(m_maxHeight, 0.0);

        std::vector<index_t> host_col_indices;

        m_childCounts.resize(m_maxHeight, 0);



        index_t current_col_index = 0;
        for (int h = 0; h < m_maxHeight; h++) {
            for (size_t i = 0; i < m_oldId[h].size(); i++) {
                NodeID old_id = m_oldId[h][i];
                index_t new_id = host_height_cutoffs[h] + i;
                host_row_offsets[new_id] = current_col_index;
                current_col_index += grg->numDownEdges(old_id);
                
                // for each down edge, we need to find its new id
                auto downEdges = grg->getDownEdges(old_id);
                m_childCounts[h] += downEdges.size();
                for (const auto& child : downEdges) {
                    release_assert(child < grg->numNodes());
                    auto child_new_id = m_newId[child];
                    index_t child_new_flat_id = host_height_cutoffs[child_new_id.first] + child_new_id.second;
                    host_col_indices.push_back(child_new_flat_id);
                }
            }
        }

        std::cout << "Host col indices constructed with size: " << host_col_indices.size() << std::endl;

        constexpr double HEAVY_RATIO = 8.0;

        for (int h = 1; h < m_maxHeight; h++) {
            index_t heavy_cutoff_value = HEAVY_RATIO * m_childCounts[h] / m_oldId[h].size();
            host_avg_child_counts[h] = static_cast<double>(m_childCounts[h]) / m_oldId[h].size();
            for (size_t i = 0; i < m_oldId[h].size(); i++) {
                NodeID old_id = m_oldId[h][i];
                index_t child_count = grg->numDownEdges(old_id);
                if (child_count < heavy_cutoff_value) {
                    host_heavy_cutoffs[h] = i;
                    break;
                }
            }
        }

        std::cout << "Host Heavy Cutoffs constructed." << std::endl;

        host_row_offsets[m_num_nodes] = current_col_index;
        release_assert(current_col_index == m_num_elements);

        auto old_to_new_flat = getNewId();

        // calculate the mutation and new id pairs
        release_assert(sizeof(index_t) >= sizeof(MutationId));
        release_assert(sizeof(index_t) >= sizeof(NodeID));

        std::vector<index_t> mutation_new_id_pairs;
        for (const auto& triple : grg->getNodesAndMutations<GRG::MutNodeMiss>()) {
            const NodeID& nodeId = std::get<0>(triple);
            const MutationId& mutId = std::get<1>(triple);

            if (nodeId == INVALID_NODE_ID) {
                mutation_new_id_pairs.push_back((index_t)INVALID_NODE_ID);
                mutation_new_id_pairs.push_back((index_t)INVALID_NODE_ID);
                continue;
            } 
            release_assert(nodeId < m_num_nodes);
            release_assert(mutId < grg->numMutations());
            const index_t& new_id = old_to_new_flat[nodeId];
            mutation_new_id_pairs.push_back((index_t)mutId);
            mutation_new_id_pairs.push_back((index_t)new_id);
        }

        size_t num_valid_mutations = mutation_new_id_pairs.size() / 2;

        std::cout << "Mutation and new id pairs constructed with size: " << mutation_new_id_pairs.size() << std::endl;
        if (num_valid_mutations != grg->numMutations()) {
            std::cout << "Warning: Valid mutations (" << num_valid_mutations << ") is less than total mutations (" << grg->numMutations() << ")." << std::endl;
        }



        std::cout << "Constructing GPUGRG with " << m_num_nodes << " nodes, " << m_num_elements << " elements, max height " << m_maxHeight << std::endl;
        gpu_grg.init(m_num_nodes, grg->numSamples(), m_num_elements, num_valid_mutations, m_maxHeight);

        gpu_grg.copyToDevice(
            host_row_offsets.data(),
            getNewId().data(),
            getOldId().data(),
            mutation_new_id_pairs.data(),
            host_height_cutoffs.data(),
            host_heavy_cutoffs.data(),
            host_avg_child_counts.data(),
            host_col_indices.data()
        );
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


/**
 * Binary file format for GPUGRG storage:
 * 
 * Header (56 bytes):
 * - uint64_t magic_number (8 bytes) - File format identifier
 * - uint64_t num_rows (8 bytes)
 * - uint64_t nnz (8 bytes)
 * - uint64_t max_height (8 bytes)
 * - uint64_t num_samples (8 bytes)
 * - uint64_t num_mutations (8 bytes)
 * - uint64_t index_size (8 bytes) - Size of index_t type in bytes
 * 
 * Data sections:
 * - row_offsets: (num_rows + 1) * sizeof(index_t) bytes
 * - col_indices: nnz * sizeof(index_t) bytes
 * - height_cutoffs: (max_height + 1) * sizeof(index_t) bytes
 * - heavy_cutoffs: (max_height + 1) * sizeof(index_t) bytes
 * - avg_child_counts: (max_height + 1) * sizeof(double) bytes
 * - old_to_new_mapping: num_rows * sizeof(index_t) bytes
 * - new_to_old_mapping: num_rows * sizeof(index_t) bytes
 */

constexpr uint64_t GPUGRG_MAGIC = 0x4752474C47505500ULL; // "GRGLGPU" + null byte

/**
 * Store a GPUGRG structure to disk in binary format.
 * 
 * @param gpu_grg The GPUGRG structure to store (must have valid CPU data)
 * @param filename Path to the output file
 * @throws std::runtime_error if file operations fail
 */
template<class index_t>
void storeGPUGRGToDisk(const GPUGRG<index_t>& gpu_grg, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    
    // Write header
    uint64_t magic = GPUGRG_MAGIC;
    uint64_t num_rows = gpu_grg.num_rows;
    uint64_t nnz = gpu_grg.nnz;
    uint64_t max_height = gpu_grg.max_height;
    uint64_t num_samples = gpu_grg.num_samples;
    uint64_t num_mutations = gpu_grg.num_mutations;
    uint64_t index_size = sizeof(index_t);

    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&num_rows), sizeof(num_rows));
    file.write(reinterpret_cast<const char*>(&nnz), sizeof(nnz));
    file.write(reinterpret_cast<const char*>(&max_height), sizeof(max_height));
    file.write(reinterpret_cast<const char*>(&num_samples), sizeof(num_samples));
    file.write(reinterpret_cast<const char*>(&num_mutations), sizeof(num_mutations));
    file.write(reinterpret_cast<const char*>(&index_size), sizeof(index_size));

    // Copy data from GPU to CPU for storage
    std::vector<index_t> host_row_offsets(num_rows + 1);
    // std::vector<index_t> host_col_indices(nnz);
    std::vector<index_t> host_old_to_new_mapping(num_rows);
    std::vector<index_t> host_new_to_old_mapping(num_rows);
    std::vector<index_t> host_mutation_new_id_pairs(num_mutations * 2);

    cudaMemcpy(host_row_offsets.data(), gpu_grg.row_offsets, 
               (num_rows + 1) * sizeof(index_t), cudaMemcpyDeviceToHost);
    //cudaMemcpy(host_col_indices.data(), gpu_grg.col_indices, 
    //           nnz * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_old_to_new_mapping.data(), gpu_grg.old_to_new_mapping,
               num_rows * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_new_to_old_mapping.data(), gpu_grg.new_to_old_mapping,
               num_rows * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_mutation_new_id_pairs.data(), gpu_grg.mutation_and_new_mapping,
               num_mutations * 2 * sizeof(index_t), cudaMemcpyDeviceToHost);
    
    // Write data sections
    file.write(reinterpret_cast<const char*>(host_row_offsets.data()), 
               (num_rows + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.host_col_indices.data()), 
               nnz * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.host_height_cutoffs.data()), 
               (max_height + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.host_heavy_cutoffs.data()), 
               (max_height + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.host_avg_child_counts.data()), 
               (max_height + 1) * sizeof(double));
    file.write(reinterpret_cast<const char*>(host_old_to_new_mapping.data()), 
               num_rows * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(host_new_to_old_mapping.data()), 
               num_rows * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(host_mutation_new_id_pairs.data()), 
               num_mutations * 2 * sizeof(index_t));
    
    if (!file.good()) {
        throw std::runtime_error("Error writing to file: " + filename);
    }
    
    file.close();
}

/**
 * Load a GPUGRG structure from disk.
 * 
 * @param filename Path to the input file
 * @return GPUGRG structure loaded from disk with GPU memory allocated and data copied
 * @throws std::runtime_error if file operations fail or format is invalid
 */
template<class index_t>
GPUGRG<index_t> loadGPUGRGFromDisk(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }
    
    // Read and validate header
    uint64_t magic, num_rows, nnz, max_height, num_samples, num_mutations, index_size;
    
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    file.read(reinterpret_cast<char*>(&num_rows), sizeof(num_rows));
    file.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));
    file.read(reinterpret_cast<char*>(&max_height), sizeof(max_height));
    file.read(reinterpret_cast<char*>(&num_samples), sizeof(num_samples));
    file.read(reinterpret_cast<char*>(&num_mutations), sizeof(num_mutations));
    file.read(reinterpret_cast<char*>(&index_size), sizeof(index_size));
    
    if (!file.good()) {
        throw std::runtime_error("Error reading file header: " + filename);
    }
    
    if (magic != GPUGRG_MAGIC) {
        throw std::runtime_error("Invalid file format: " + filename);
    }

    if (index_size != sizeof(index_t)) {
        throw std::runtime_error("Incompatible index or data type size in file: " + filename);
    }
    
    // Create and allocate GPUGRG structure
    GPUGRG<index_t> gpu_grg;
    gpu_grg.init(num_rows, num_samples, nnz, num_mutations, max_height);
    
    // Read data sections
    std::vector<index_t> host_row_offsets(num_rows + 1);
    std::vector<index_t> host_col_indices(nnz);
    std::vector<index_t> host_height_cutoffs(max_height + 1);
    std::vector<index_t> host_heavy_cutoffs(max_height + 1);
    std::vector<double> host_avg_child_counts(max_height + 1);
    std::vector<index_t> old_to_new_mapping(num_rows);
    std::vector<index_t> new_to_old_mapping(num_rows);
    std::vector<index_t> host_mutation_new_id_pairs(num_mutations * 2);

    file.read(reinterpret_cast<char*>(host_row_offsets.data()), 
              (num_rows + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(host_col_indices.data()), 
              nnz * sizeof(index_t));
    file.read(reinterpret_cast<char*>(host_height_cutoffs.data()), 
              (max_height + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(host_heavy_cutoffs.data()), 
              (max_height + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(host_avg_child_counts.data()), 
              (max_height + 1) * sizeof(double));
    file.read(reinterpret_cast<char*>(old_to_new_mapping.data()), 
              num_rows * sizeof(index_t));
    file.read(reinterpret_cast<char*>(new_to_old_mapping.data()), 
              num_rows * sizeof(index_t));
    file.read(reinterpret_cast<char*>(host_mutation_new_id_pairs.data()), 
              num_mutations * 2 * sizeof(index_t));
    
    if (!file.good()) {
        throw std::runtime_error("Error reading file data: " + filename);
    }
    
    file.close();
    
    // Copy data to GPU
    gpu_grg.copyToDevice(
        host_row_offsets.data(),
        old_to_new_mapping.data(),
        new_to_old_mapping.data(),
        host_mutation_new_id_pairs.data(),
        host_height_cutoffs.data(),
        host_heavy_cutoffs.data(),
        host_avg_child_counts.data(),
        host_col_indices.data()
    );
    
    return gpu_grg;
}

/**
 * Convert a GRG to GPUGRG format for GPU matrix multiplication.
 * This function performs the same conversion that was previously done inline
 * in matrixMultiplicationGPU, but as a standalone reusable function.
 * 
 * @param grg The GRG to convert (must have ordered nodes)
 * @return GPUGRG structure with GPU memory allocated and data copied
 * @throws std::runtime_error if GRG doesn't meet requirements
 */
template<class index_t>
GPUGRG<index_t> convertGRGToGPUGRG(GRG* grg) {
    // Validate input requirements
    if (!grg->nodesAreOrdered()) {
        throw std::runtime_error("GRG nodes must be topologically ordered for GPU conversion");
    }
    
    if (grg->numNodes() == 0 || grg->numEdges() == 0) {
        throw std::runtime_error("GRG must have nodes and edges for GPU conversion");
    }
    
    // Check for potential overflow with index_t
    if (grg->numNodes() > std::numeric_limits<index_t>::max() || 
        grg->numEdges() > std::numeric_limits<index_t>::max()) {
        throw std::runtime_error("GRG too large for specified index type");
    }
    
    std::cout << "Converting GRG to GPUGRG format..." << std::endl;
    auto time_start = std::chrono::high_resolution_clock::now();
    
    // Create visitor to compute height-based ordering
    GPUCsrVisitor<index_t> visitor(grg->numNodes(), grg->numEdges());
    auto seeds = grg->getRootNodes();
    
    // Perform DFS traversal to compute node heights
    grg->visitDfs(visitor, TraversalDirection::DIRECTION_DOWN, seeds);

    // Rearrange nodes within each height by child counts for better locality
    visitor.rearrange(grg);
    
    // Construct the GPUGRG structure
    GPUGRG<index_t> gpu_grg;
    visitor.constructGPUGRG(grg, gpu_grg);
    
    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    
    std::cout << "Converted GRG to GPUGRG format with " << gpu_grg.num_rows << " rows and " 
              << gpu_grg.nnz << " non-zeros." << std::endl;
    std::cout << "Time taken for conversion: " << duration << " ms" << std::endl;
    
    // Optional: print child counts per height for debugging
    if (std::getenv("GRGL_GPU_DEBUG")) {
        visitor.printChildCounts();
        visitor.printNodeCounts();
    }
    
    return gpu_grg;
}

} // namespace grgl

#endif // GRGL_CUDA_ENABLED
#endif // GRGL_GPU_GRG_H