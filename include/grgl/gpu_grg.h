#ifndef GRGL_GPU_GRG_H
#define GRGL_GPU_GRG_H

#ifndef GRGL_CUDA_ENABLED
#error "Cannot include gpu_grg.h if CUDA is not enabled; rebuild with GRGL_CUDA_ENABLED defined."
#endif

#define CHECK_CUDA_LAST_ERROR()                                                                                        \
    {                                                                                                                  \
        cudaError_t err = cudaGetLastError();                                                                          \
        if (err != cudaSuccess) {                                                                                      \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " code=" << static_cast<int>(err) << " \"" \
                      << cudaGetErrorString(err) << "\"" << std::endl;                                                 \
            throw std::runtime_error("CUDA kernel launch failed");                                                     \
        }                                                                                                              \
    }

#include "grgl/common.h"
#include "grgl/grg.h"
#include <chrono>
#include <cuda_runtime.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace grgl {

template <class index_t> class GPUGRG {
public:
    // these are resident GPU pointers
    index_t* rowOffsets; // Device pointer to row offsets (size: num_rows + 1)
    index_t* oldToNewMapping;
    index_t* newToOldMapping;       // This is not used actually
    index_t* mutationAndNewMapping; // Mapping between mutation id and new node id (in pairs)

    // Since col indices may be large, this data may not be resident on GPU
    index_t* colIndices; // Device pointer to column indices (size: nnz)

    // This is the corresponding host backup data
    std::vector<index_t> hostColIndices;

    std::vector<index_t> hostHeightCutoffs;
    std::vector<index_t> hostHeavyCutoffs;
    std::vector<double> hostAvgChildCounts;

    size_t numRows;      // Number of rows, i.e. number of nodes
    size_t numSamples;   // Number of samples, i.e. number of leaf nodes
    size_t numMutations; // Number of mutations. may NOT correspond to number of non-leaf nodes
    size_t numEdges;           // Number of non-zero elements, i.e. number of edges
    size_t maxHeight;

    cudaStream_t* workStreamPtr; // CUDA stream for asynchronous operations

public:
    // Constructor
    GPUGRG()
        : rowOffsets(nullptr),
          colIndices(nullptr),
          oldToNewMapping(nullptr),
          newToOldMapping(nullptr),
          mutationAndNewMapping(nullptr),
          numRows(0),
          numSamples(0),
          numEdges(0),
          maxHeight(0),
          numMutations(0),
          workStreamPtr(nullptr) {}

    // Destructor
    ~GPUGRG() { free(); }

    // return the size of GPU memory needed by col indices in bytes
    size_t queryColIndSize() { return numEdges * sizeof(index_t); }

    // This function sets existing GPU memory pointers for col indices
    // will not allocate new memory
    void setColIndBuffer(void* ptr, size_t size) {
        if (size < queryColIndSize()) {
            throw ApiMisuseFailure("Provided buffer size is smaller than required col indices size");
        }
        if (colIndices) {
            throw ApiMisuseFailure("Col indices buffer is already set");
        }
        colIndices = static_cast<index_t*>(ptr);
    }

    // This function copies col indices from Host to GPU
    // Asynchronous copy is supported
    // The caller must ensure the memory space is not used by other operations.
    void copyColIndFromHost() {
        if (colIndices) {
            cudaMemcpyAsync(
                colIndices, hostColIndices.data(), numEdges * sizeof(index_t), cudaMemcpyHostToDevice, *workStreamPtr);
        } else {
            throw ApiMisuseFailure("No col indices buffer set for GPUGRG copyColIndFromHost operation");
        }
    }

    // This function will block the calling thread
    // until all preceding operations related to this GRG is done
    void wait() {
        if (workStreamPtr) {
            cudaStreamSynchronize(*workStreamPtr);
        } else {
            throw ApiMisuseFailure("No CUDA stream set for GPUGRG wait operation");
        }
    }

    // This function removes the GPU memory buffer
    // Will not free the memory, just removes the pointers
    void removeColIndBuffer() { colIndices = nullptr; }

    // This function frees the GPU memory buffer for col indices
    void freeColInd() {
        if (colIndices == nullptr) {
            throw ApiMisuseFailure("Col indices buffer is already freed");
        }
        cudaFree(colIndices);
        colIndices = nullptr;
    }

    // Init vars and allocate GPU memory for all structures
    // This function is blocking
    // If no streamPtr is provided, a new stream will be created
    void init(size_t rows,
              size_t samples,
              size_t nonZeros,
              size_t mutations,
              size_t height,
              bool allocateColIndices = true,
              cudaStream_t* streamPtr = nullptr) {
        free(); // Free any existing memory
        CHECK_CUDA_LAST_ERROR();
        numRows = rows;
        numSamples = samples;
        numEdges = nonZeros;
        numMutations = mutations;
        maxHeight = height;

        cudaMalloc(&rowOffsets, (numRows + 1) * sizeof(index_t));
        cudaMalloc(&oldToNewMapping, numRows * sizeof(index_t));
        cudaMalloc(&newToOldMapping, numRows * sizeof(index_t));
        cudaMalloc(&mutationAndNewMapping, numMutations * 2 * sizeof(index_t)); // Allocate for pairs
        CHECK_CUDA_LAST_ERROR();

        if (allocateColIndices) {
            cudaMalloc(&colIndices, numEdges * sizeof(index_t));
        } else {
            colIndices = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();

        hostColIndices.resize(numEdges);
        hostHeightCutoffs.resize(maxHeight + 1);
        hostHeavyCutoffs.resize(maxHeight + 1);
        hostAvgChildCounts.resize(maxHeight + 1);

        if (streamPtr) {
            workStreamPtr = streamPtr;
        } else {
            workStreamPtr = new cudaStream_t();
            cudaStreamCreate(workStreamPtr);
        }
        CHECK_CUDA_LAST_ERROR();
    }

    // Free GPU memory
    void free() {
        if (rowOffsets) {
            cudaFree(rowOffsets);
            rowOffsets = nullptr;
        }
        if (colIndices) {
            cudaFree(colIndices);
            colIndices = nullptr;
        }
        if (oldToNewMapping) {
            cudaFree(oldToNewMapping);
            oldToNewMapping = nullptr;
        }
        if (newToOldMapping) {
            cudaFree(newToOldMapping);
            newToOldMapping = nullptr;
        }
        if (mutationAndNewMapping) {
            cudaFree(mutationAndNewMapping);
            mutationAndNewMapping = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();

        hostColIndices.clear();
        hostHeightCutoffs.clear();
        hostHeavyCutoffs.clear();
        hostAvgChildCounts.clear();

        numRows = numSamples = numEdges = 0;
        numMutations = 0;
        maxHeight = 0;

        if (workStreamPtr) {
            cudaStreamDestroy(*workStreamPtr);
            delete workStreamPtr;
            workStreamPtr = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();
    }

    // Copy data to GPU
    // This function is blocking
    void copyToDevice(const index_t* rowOffsets,
                      const index_t* oldToNewMapping,
                      const index_t* newToOldMapping,
                      const index_t* mutationAndNewMapping,
                      const index_t* heightCutoffs,
                      const index_t* heavyCutoffs,
                      const double* avgChildCounts,
                      const index_t* colIndices,
                      bool copyColIndicesToDevice = true) {

        if (!this->rowOffsets || !this->oldToNewMapping || !this->newToOldMapping) {
            throw ApiMisuseFailure("Device memory not allocated yet!");
        }
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->rowOffsets, rowOffsets, (numRows + 1) * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->oldToNewMapping, oldToNewMapping, numRows * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->newToOldMapping, newToOldMapping, numRows * sizeof(index_t), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->mutationAndNewMapping,
                   mutationAndNewMapping,
                   numMutations * 2 * sizeof(index_t),
                   cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();

        // cpu memcpy just use normal memcpys
        memcpy(this->hostHeightCutoffs.data(), heightCutoffs, (maxHeight + 1) * sizeof(index_t));
        memcpy(this->hostHeavyCutoffs.data(), heavyCutoffs, (maxHeight + 1) * sizeof(index_t));
        memcpy(this->hostAvgChildCounts.data(), avgChildCounts, (maxHeight + 1) * sizeof(double));
        memcpy(this->hostColIndices.data(), colIndices, numEdges * sizeof(index_t));
        CHECK_CUDA_LAST_ERROR();

        if (copyColIndicesToDevice) {
            if (!this->colIndices) {
                throw ApiMisuseFailure("Column indices GPU buffer is not allocated");
            }
            cudaMemcpy(this->colIndices, colIndices, numEdges * sizeof(index_t), cudaMemcpyHostToDevice);
        }
        CHECK_CUDA_LAST_ERROR();
    }

    // Debugging: Print CSR structure
    void print() const {
        std::cout << "GPUGRG: " << std::endl;
        std::cout << "  num_rows: " << numRows << std::endl;
        std::cout << "  nnz: " << numEdges << std::endl;
        std::cout << "  max_height: " << maxHeight << std::endl;
        std::cout << "  num_samples: " << numSamples << std::endl;
        std::cout << "  num_mutations: " << numMutations << std::endl;

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
    std::vector<T>
    matMulBlocking(const std::vector<T>& inputMatrix, const size_t numRows, TraversalDirection direction) {
        CHECK_CUDA_LAST_ERROR();
        if (numRows == 0 || (inputMatrix.size() % numRows != 0)) {
            throw ApiMisuseFailure("inputMatrix must be divisible by numRows");
        }
        const size_t numCols = inputMatrix.size() / numRows;
        const size_t outSize =
            (numRows * ((direction == TraversalDirection::DIRECTION_DOWN) ? this->numSamples : this->numMutations));

        T* d_inputMatrix;
        T* d_outputMatrix;
        cudaMalloc(&d_inputMatrix, inputMatrix.size() * sizeof(T));
        cudaMalloc(&d_outputMatrix, outSize * sizeof(T));
        cudaMemcpy(d_inputMatrix, inputMatrix.data(), inputMatrix.size() * sizeof(T), cudaMemcpyHostToDevice);

        CHECK_CUDA_LAST_ERROR();

        T* d_buffer;
        size_t bufferSizeByte = numRows * this->numRows * sizeof(T);
        cudaMalloc(&d_buffer, bufferSizeByte);

        CHECK_CUDA_LAST_ERROR();

        this->matrixMultiplication<T>(
            d_inputMatrix, numCols, numRows, direction, d_outputMatrix, outSize, false, d_buffer, bufferSizeByte);
        this->wait();
        CHECK_CUDA_LAST_ERROR();

        std::vector<T> result(outSize);
        cudaMemcpy(result.data(), d_outputMatrix, outSize * sizeof(T), cudaMemcpyDeviceToHost);
        cudaFree(d_inputMatrix);
        cudaFree(d_outputMatrix);
        cudaFree(d_buffer);
        CHECK_CUDA_LAST_ERROR();

        return std::move(result);
    }

    template <typename T>
    std::vector<T> matMulPerf(const std::vector<T>& inputMatrix,
                              const size_t numRows,
                              TraversalDirection direction,
                              int iterations = 5) {
        CHECK_CUDA_LAST_ERROR();
        if (numRows == 0 || (inputMatrix.size() % numRows != 0)) {
            throw ApiMisuseFailure("inputMatrix must be divisible by numRows");
        }
        const size_t numCols = inputMatrix.size() / numRows;
        const size_t outSize =
            (numRows * ((direction == TraversalDirection::DIRECTION_DOWN) ? this->numSamples : this->numMutations));

        T* d_inputMatrix;
        T* d_outputMatrix;
        cudaMalloc(&d_inputMatrix, inputMatrix.size() * sizeof(T));
        cudaMalloc(&d_outputMatrix, outSize * sizeof(T));
        cudaMemcpy(d_inputMatrix, inputMatrix.data(), inputMatrix.size() * sizeof(T), cudaMemcpyHostToDevice);

        CHECK_CUDA_LAST_ERROR();

        T* d_buffer;
        size_t bufferSizeByte = numRows * this->numRows * sizeof(T);
        cudaMalloc(&d_buffer, bufferSizeByte);

        CHECK_CUDA_LAST_ERROR();

        auto start = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; it++) {
            this->matrixMultiplication<T>(
                d_inputMatrix, numCols, numRows, direction, d_outputMatrix, outSize, false, d_buffer, bufferSizeByte);
            this->wait();
            CHECK_CUDA_LAST_ERROR();
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Average end-to-end runtime for matrixMultiplication function (GPU) over " << iterations
                  << " iterations: " << (elapsed.count() / iterations) << " ms" << std::endl;

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
                              size_t buffer_size);
};

namespace {
template <class index_t> struct nodes {
    index_t id;
    index_t childCounts;
};

template <class index_t> bool cmp(const nodes<index_t>& a, const nodes<index_t>& b) {
    return a.childCounts > b.childCounts;
}
} // namespace

// This visitor is tested ONLY on CSRGRG (immutable, ordered GRG)
// in the first visitor path, we get the height of each node, and thus the new id. new id is of (height, gid)
// we construct both old_id -> new_id map and a new_id -> old_id map
// then, in the epilogue function, we can use the information to have a new GPUGRG
template <class index_t> class GPUCsrVisitor : public GRGVisitor {
public:
    explicit GPUCsrVisitor(size_t num_nodes, size_t num_elements) {
        h_maxHeight = 0;
        h_numNodes = num_nodes;
        h_numEdges = num_elements;
        h_newID.resize(num_nodes, std::make_pair(0, 0));
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
            return true; // do nothing on the forward pass
        }

        index_t height = 0;

        if (nodeId < grg->numSamples()) {
            height = 0;
        } else {
            auto downEdges = grg->getDownEdges(nodeId);
            for (const auto& child : downEdges) {
                release_assert(child < grg->numNodes());
                auto childHeight = h_newID[child].first;
                height = std::max(height, childHeight + 1);
            }
        }

        if (height >= h_maxHeight) {
            h_maxHeight = height + 1;
            // m_idCounters.resize(m_maxHeight, 0);
            h_oldId.resize(h_maxHeight);
        }

        // auto id_in_group = m_idCounters[height]++;
        auto idInGroup = h_oldId[height].size();
        h_oldId[height].push_back(nodeId);
        h_newID[nodeId] = std::make_pair(height, idInGroup);
        return true;
    }

    // Flattened mapping from new_id to old_id
    std::vector<index_t> getOldId() const {
        std::vector<index_t> flatOldID;
        for (const auto& vec : h_oldId) {
            flatOldID.insert(flatOldID.end(), vec.begin(), vec.end());
        }
        return flatOldID;
    }

    // Flattened mapping from old_id to new_id
    std::vector<index_t> getNewId() const {
        std::vector<index_t> flatNewID;
        flatNewID.resize(h_numNodes);
        auto baseID = 0;
        for (size_t height = 0; height < h_oldId.size(); height++) {
            for (size_t i = 0; i < h_oldId[height].size(); i++) {
                NodeID old_id = h_oldId[height][i];
                index_t new_id = baseID + i;
                flatNewID[old_id] = new_id;
            }
            baseID += h_oldId[height].size();
        }
        return flatNewID;
    }

    // this function rearranges nodes at each height by their child counts
    void rearrange(GRG* grg) {
        for (int h = 1; h < h_maxHeight; h++) {
            std::cout << "Rearranging height " << h << " with " << h_oldId[h].size() << " nodes." << std::endl;
            std::vector<nodes<index_t>> nodeList;
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                NodeID oldID = h_oldId[h][i];
                index_t childCount = grg->numDownEdges(oldID);
                nodeList.push_back({oldID, childCount});
            }
            std::sort(nodeList.begin(), nodeList.end(), cmp<index_t>);
            for (size_t i = 0; i < nodeList.size(); i++) {
                h_oldId[h][i] = nodeList[i].id;
                NodeID old_id = nodeList[i].id;
                h_newID[old_id] = std::make_pair(h, i);
            }
        }
    }

    // After visiting all nodes, we can construct the CSR representation
    // This function should be called only once after all visits are done
    void constructGPUGRG(GRG* grg, GPUGRG<index_t>& gpu_grg) {

        std::vector<index_t> h_heightCutoffs;
        h_heightCutoffs.resize(h_maxHeight + 1, 0);
        for (int i = 1; i <= h_maxHeight; i++) {
            h_heightCutoffs[i] = h_heightCutoffs[i - 1] + h_oldId[i - 1].size();
        }

        std::vector<index_t> h_heavyCutoffs;
        h_heavyCutoffs.resize(h_maxHeight + 1, 0);

        std::vector<index_t> h_rowOffsets;
        h_rowOffsets.resize(h_numNodes + 1, 0);

        std::vector<double> h_avgChildCounts;
        h_avgChildCounts.resize(h_maxHeight, 0.0);

        std::vector<index_t> h_colIndices;

        h_childCounts.resize(h_maxHeight, 0);

        index_t currrentColIndex = 0;
        for (int h = 0; h < h_maxHeight; h++) {
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                NodeID oldID = h_oldId[h][i];
                index_t newID = h_heightCutoffs[h] + i;
                h_rowOffsets[newID] = currrentColIndex;
                currrentColIndex += grg->numDownEdges(oldID);

                // for each down edge, we need to find its new id
                auto downEdges = grg->getDownEdges(oldID);
                h_childCounts[h] += downEdges.size();
                for (const auto& child : downEdges) {
                    release_assert(child < grg->numNodes());
                    auto childNewID = h_newID[child];
                    index_t childNewFlatID = h_heightCutoffs[childNewID.first] + childNewID.second;
                    h_colIndices.push_back(childNewFlatID);
                }
            }
        }

        std::cout << "Host col indices constructed with size: " << h_colIndices.size() << std::endl;

        constexpr double HEAVY_RATIO = 8.0; // TODO: define this in another place

        for (int h = 1; h < h_maxHeight; h++) {
            index_t heavyCutoffValue = HEAVY_RATIO * h_childCounts[h] / h_oldId[h].size();
            h_avgChildCounts[h] = static_cast<double>(h_childCounts[h]) / h_oldId[h].size();
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                NodeID oldID = h_oldId[h][i];
                index_t childCount = grg->numDownEdges(oldID);
                if (childCount < heavyCutoffValue) {
                    h_heavyCutoffs[h] = i;
                    break;
                }
            }
        }

        std::cout << "Host Heavy Cutoffs constructed." << std::endl;

        h_rowOffsets[h_numNodes] = currrentColIndex;
        release_assert(currrentColIndex == h_numEdges);

        auto oldToNewFlat = getNewId();

        // calculate the mutation and new id pairs
        release_assert(sizeof(index_t) >= sizeof(MutationId));
        release_assert(sizeof(index_t) >= sizeof(NodeID));

        std::vector<index_t> mutationNewPairs;
        for (const auto& triple : grg->getNodesAndMutations<GRG::MutNodeMiss>()) {
            const NodeID& nodeId = std::get<0>(triple);
            const MutationId& mutId = std::get<1>(triple);

            if (nodeId == INVALID_NODE_ID) {
                mutationNewPairs.push_back((index_t)INVALID_NODE_ID);
                mutationNewPairs.push_back((index_t)INVALID_NODE_ID);
                continue;
            }
            release_assert(nodeId < h_numNodes);
            release_assert(mutId < grg->numMutations());
            const index_t& newID = oldToNewFlat[nodeId];
            mutationNewPairs.push_back((index_t)mutId);
            mutationNewPairs.push_back((index_t)newID);
        }

        size_t numValidMutations = mutationNewPairs.size() / 2;

        std::cout << "Mutation and new id pairs constructed with size: " << mutationNewPairs.size() << std::endl;
        if (numValidMutations != grg->numMutations()) {
            std::cout << "Warning: Valid mutations (" << numValidMutations << ") is less than total mutations ("
                      << grg->numMutations() << ")." << std::endl;
        }

        std::cout << "Constructing GPUGRG with " << h_numNodes << " nodes, " << h_numEdges
                  << " elements, max height " << h_maxHeight << std::endl;
        gpu_grg.init(h_numNodes, grg->numSamples(), h_numEdges, numValidMutations, h_maxHeight);

        gpu_grg.copyToDevice(h_rowOffsets.data(),
                             getNewId().data(),
                             getOldId().data(),
                             mutationNewPairs.data(),
                             h_heightCutoffs.data(),
                             h_heavyCutoffs.data(),
                             h_avgChildCounts.data(),
                             h_colIndices.data());
    }

    void printChildCounts() const {
        std::cout << "Child counts per height:" << std::endl;
        for (size_t h = 0; h < h_childCounts.size(); h++) {
            std::cout << "Height " << h << ": " << h_childCounts[h] << std::endl;
        }
    }

    void printNodeCounts() const {
        std::cout << "Node counts per height:" << std::endl;
        for (size_t h = 0; h < h_oldId.size(); h++) {
            std::cout << "Height " << h << ": " << h_oldId[h].size() << std::endl;
        }
    }

private:
    index_t h_maxHeight;
    size_t h_numNodes;
    size_t h_numEdges;

    std::vector<size_t> h_childCounts;                // sum of child counts for nodes of the same height
    std::vector<std::vector<index_t>> h_oldId;        // new_id -> old_id mapping
    std::vector<std::pair<index_t, index_t>> h_newID; // old_id -> new_id mapping
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

constexpr uint64_t GPUGRG_MAGIC = 0x4752474C47505500ULL; // "GRGLGPU" + version byte 00

/**
 * Store a GPUGRG structure to disk in binary format.
 *
 * @param gpu_grg The GPUGRG structure to store (must have valid CPU data)
 * @param filename Path to the output file
 * @throws std::runtime_error if file operations fail
 */
template <class index_t> void storeGPUGRGToDisk(const GPUGRG<index_t>& gpu_grg, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }

    // Write header
    uint64_t magic = GPUGRG_MAGIC;
    uint64_t numRows = gpu_grg.numRows;
    uint64_t numEdges = gpu_grg.numEdges;
    uint64_t maxHeight = gpu_grg.maxHeight;
    uint64_t numSamples = gpu_grg.numSamples;
    uint64_t numMutations = gpu_grg.numMutations;
    uint64_t indexSize = sizeof(index_t);

    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&numRows), sizeof(numRows));
    file.write(reinterpret_cast<const char*>(&numEdges), sizeof(numEdges));
    file.write(reinterpret_cast<const char*>(&maxHeight), sizeof(maxHeight));
    file.write(reinterpret_cast<const char*>(&numSamples), sizeof(numSamples));
    file.write(reinterpret_cast<const char*>(&numMutations), sizeof(numMutations));
    file.write(reinterpret_cast<const char*>(&indexSize), sizeof(indexSize));

    // Copy data from GPU to CPU for storage
    std::vector<index_t> h_RowOffsets(numRows + 1);
    // std::vector<index_t> host_col_indices(nnz);
    std::vector<index_t> h_oldToNewMapping(numRows);
    std::vector<index_t> h_newToOldMapping(numRows);
    std::vector<index_t> h_mutationNewPairs(numMutations * 2);

    cudaMemcpy(h_RowOffsets.data(), gpu_grg.rowOffsets, (numRows + 1) * sizeof(index_t), cudaMemcpyDeviceToHost);
    // cudaMemcpy(host_col_indices.data(), gpu_grg.col_indices,
    //            nnz * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(
        h_oldToNewMapping.data(), gpu_grg.oldToNewMapping, numRows * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(
        h_newToOldMapping.data(), gpu_grg.newToOldMapping, numRows * sizeof(index_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mutationNewPairs.data(),
               gpu_grg.mutationAndNewMapping,
               numMutations * 2 * sizeof(index_t),
               cudaMemcpyDeviceToHost);

    // Write data sections
    file.write(reinterpret_cast<const char*>(h_RowOffsets.data()), (numRows + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostColIndices.data()), numEdges * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostHeightCutoffs.data()), (maxHeight + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostHeavyCutoffs.data()), (maxHeight + 1) * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostAvgChildCounts.data()), (maxHeight + 1) * sizeof(double));
    file.write(reinterpret_cast<const char*>(h_oldToNewMapping.data()), numRows * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(h_newToOldMapping.data()), numRows * sizeof(index_t));
    file.write(reinterpret_cast<const char*>(h_mutationNewPairs.data()), numMutations * 2 * sizeof(index_t));

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
template <class index_t> GPUGRG<index_t> loadGPUGRGFromDisk(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading: " + filename);
    }

    // Read and validate header
    uint64_t magic, numRows, numEdges, maxHeight, numSamples, numMutations, indexSize;

    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    file.read(reinterpret_cast<char*>(&numRows), sizeof(numRows));
    file.read(reinterpret_cast<char*>(&numEdges), sizeof(numEdges));
    file.read(reinterpret_cast<char*>(&maxHeight), sizeof(maxHeight));
    file.read(reinterpret_cast<char*>(&numSamples), sizeof(numSamples));
    file.read(reinterpret_cast<char*>(&numMutations), sizeof(numMutations));
    file.read(reinterpret_cast<char*>(&indexSize), sizeof(indexSize));

    if (!file.good()) {
        throw std::runtime_error("Error reading file header: " + filename);
    }

    if (magic != GPUGRG_MAGIC) {
        throw std::runtime_error("Invalid file format: " + filename);
    }

    if (indexSize != sizeof(index_t)) {
        throw std::runtime_error("Incompatible index or data type size in file: " + filename);
    }

    // Create and allocate GPUGRG structure
    GPUGRG<index_t> gpuGRG;
    gpuGRG.init(numRows, numSamples, numEdges, numMutations, maxHeight);

    // Read data sections
    std::vector<index_t> h_rowOffsets(numRows + 1);
    std::vector<index_t> h_colIndices(numEdges);
    std::vector<index_t> h_heightCutoffs(maxHeight + 1);
    std::vector<index_t> h_heavyCutoffs(maxHeight + 1);
    std::vector<double> h_avgChildCounts(maxHeight + 1);
    std::vector<index_t> h_oldToNewMapping(numRows);
    std::vector<index_t> h_newToOldMapping(numRows);
    std::vector<index_t> h_mutationNewPairs(numMutations * 2);

    file.read(reinterpret_cast<char*>(h_rowOffsets.data()), (numRows + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_colIndices.data()), numEdges * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_heightCutoffs.data()), (maxHeight + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_heavyCutoffs.data()), (maxHeight + 1) * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_avgChildCounts.data()), (maxHeight + 1) * sizeof(double));
    file.read(reinterpret_cast<char*>(h_oldToNewMapping.data()), numRows * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_newToOldMapping.data()), numRows * sizeof(index_t));
    file.read(reinterpret_cast<char*>(h_mutationNewPairs.data()), numMutations * 2 * sizeof(index_t));

    if (!file.good()) {
        throw std::runtime_error("Error reading file data: " + filename);
    }

    file.close();

    // Copy data to GPU
    gpuGRG.copyToDevice(h_rowOffsets.data(),
                         h_oldToNewMapping.data(),
                         h_newToOldMapping.data(),
                         h_mutationNewPairs.data(),
                         h_heightCutoffs.data(),
                         h_heavyCutoffs.data(),
                         h_avgChildCounts.data(),
                         h_colIndices.data());

    CHECK_CUDA_LAST_ERROR();

    return gpuGRG;
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
template <class index_t> GPUGRG<index_t> convertGRGToGPUGRG(GRG* grg) {
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

    std::cout << "Converted GRG to GPUGRG format with " << gpu_grg.numRows << " rows and " << gpu_grg.numEdges
              << " non-zeros." << std::endl;
    std::cout << "Time taken for conversion: " << duration << " ms" << std::endl;

    // Optional: print child counts per height for debugging
    if (std::getenv("GRGL_GPU_DEBUG")) {
        visitor.printChildCounts();
        visitor.printNodeCounts();
    }

    return gpu_grg;
}

bool hasCudaSupport();

} // namespace grgl

#endif // GRGL_GPU_GRG_H