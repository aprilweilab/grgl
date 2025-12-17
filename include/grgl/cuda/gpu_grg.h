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
#include "grgl/grgnode.h"
#include <chrono>
#include <cuda_runtime.h>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace grgl {

template <typename T> class CudaBuffer {
public:
    CudaBuffer()
        : m_buffer(nullptr) {}

    CudaBuffer(size_t dim1) {
        // Check for potential overflow before multiplication
        if (dim1 > std::numeric_limits<size_t>::max() / sizeof(T)) {
            throw std::runtime_error("CudaBuffer allocation size would overflow");
        }
        const size_t allocSize = dim1 * sizeof(T);
        cudaMalloc(&m_buffer, allocSize);
        CHECK_CUDA_LAST_ERROR();
    }

    CudaBuffer(size_t dim1, size_t dim2) {
        // Check for potential overflow before multiplication
        if (dim1 == 0 || dim2 > std::numeric_limits<size_t>::max() / sizeof(T) / dim1) {
            throw std::runtime_error("CudaBuffer allocation size would overflow");
        }
        const size_t allocSize = dim1 * dim2 * sizeof(T);
        cudaMalloc(&m_buffer, allocSize);
        CHECK_CUDA_LAST_ERROR();
    }

    ~CudaBuffer() {
        if (m_buffer != nullptr) {
            cudaFree(m_buffer);
            m_buffer = nullptr;
        }
    }

    T* get() { return m_buffer; }

    // Disable all copying of this type.
    CudaBuffer(const CudaBuffer& rhs) = delete;
    CudaBuffer& operator=(const CudaBuffer& rhs) = delete;
    CudaBuffer& operator=(CudaBuffer&& rhs) noexcept = delete;
    CudaBuffer(CudaBuffer&& rhs) noexcept = delete;

private:
    T* m_buffer;
};

typedef std::shared_ptr<CudaBuffer<NodeIDSizeT>> CudaNodeIDBufferPtr;

constexpr uint64_t GPUGRG_MAGIC = 0x4752474C47505501ULL; // "GRGLGPU" + version byte 01

class GPUGRG {
public:
    // friend GPUGRG loadGPUGRGFromDisk(const std::string& filename);
    // friend void storeGPUGRGToDisk(const GPUGRG& gpu_grg, const std::string& filename);
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

    // Init vars and allocate GPU memory for all structures
    // This function is blocking
    // If no streamPtr is provided, a new stream will be created
    void init(size_t rows,
              size_t samples,
              size_t nonZeros,
              size_t mutations,
              size_t height,
              bool allocateColIndices = true,
              std::shared_ptr<cudaStream_t> streamPtr = nullptr) {
        free(); // Free any existing memory
        numRows = rows;
        numSamples = samples;
        numEdges = nonZeros;
        numMutations = mutations;
        maxHeight = height;

        rowOffsets = std::make_shared<CudaBuffer<NodeIDSizeT>>(numRows + 1);
        oldToNewMapping = std::make_shared<CudaBuffer<NodeIDSizeT>>(numRows);
        newToOldMapping = std::make_shared<CudaBuffer<NodeIDSizeT>>(numRows);
        mutationAndNewMapping = std::make_shared<CudaBuffer<NodeIDSizeT>>(numMutations, 2);

        if (allocateColIndices) {
            colIndices = std::make_shared<CudaBuffer<NodeIDSizeT>>(numEdges);
        } else {
            colIndices = nullptr;
        }

        hostColIndices.resize(numEdges);
        hostHeightCutoffs.resize(maxHeight + 1);
        hostHeavyCutoffs.resize(maxHeight + 1);
        hostAvgChildCounts.resize(maxHeight + 1);

        if (streamPtr) {
            workStreamPtr = streamPtr;
        } else {
            workStreamPtr = std::make_shared<cudaStream_t>();
            cudaStreamCreate(workStreamPtr.get());
        }
        CHECK_CUDA_LAST_ERROR();
    }

    // Copy data to GPU
    // This function is blocking
    void copyToDevice(const NodeIDSizeT* rowOffsets,
                      const NodeIDSizeT* oldToNewMapping,
                      const NodeIDSizeT* newToOldMapping,
                      const NodeIDSizeT* mutationAndNewMapping,
                      const NodeIDSizeT* heightCutoffs,
                      const NodeIDSizeT* heavyCutoffs,
                      const double* avgChildCounts,
                      const NodeIDSizeT* colIndices,
                      bool copyColIndicesToDevice = true) {

        if (!this->rowOffsets || !this->oldToNewMapping || !this->newToOldMapping) {
            throw ApiMisuseFailure("Device memory not allocated yet!");
        }
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->getRowOffsets(), rowOffsets, (numRows + 1) * sizeof(NodeIDSizeT), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->getOldToNewMapping(), oldToNewMapping, numRows * sizeof(NodeIDSizeT), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->getNewToOldMapping(), newToOldMapping, numRows * sizeof(NodeIDSizeT), cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();
        cudaMemcpy(this->getMutationAndNewMapping(),
                   mutationAndNewMapping,
                   numMutations * 2 * sizeof(NodeIDSizeT),
                   cudaMemcpyHostToDevice);
        CHECK_CUDA_LAST_ERROR();

        // cpu memcpy just use normal memcpys
        release_assert(this->hostHeightCutoffs.size() == maxHeight + 1);
        memcpy(this->hostHeightCutoffs.data(), heightCutoffs, (maxHeight + 1) * sizeof(NodeIDSizeT));
        release_assert(this->hostHeavyCutoffs.size() == maxHeight + 1);
        memcpy(this->hostHeavyCutoffs.data(), heavyCutoffs, (maxHeight + 1) * sizeof(NodeIDSizeT));
        release_assert(this->hostAvgChildCounts.size() == maxHeight + 1);
        memcpy(this->hostAvgChildCounts.data(), avgChildCounts, (maxHeight + 1) * sizeof(double));
        release_assert(this->hostColIndices.size() == numEdges);
        memcpy(this->hostColIndices.data(), colIndices, numEdges * sizeof(NodeIDSizeT));

        if (copyColIndicesToDevice) {
            if (!this->colIndices) {
                throw ApiMisuseFailure("Column indices GPU buffer is not allocated");
            }
            cudaMemcpy(this->getColIndices(), colIndices, numEdges * sizeof(NodeIDSizeT), cudaMemcpyHostToDevice);
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

        CudaBuffer<T> d_inputMatrix(inputMatrix.size());
        CudaBuffer<T> d_outputMatrix(outSize);
        cudaMemcpy(d_inputMatrix.get(), inputMatrix.data(), inputMatrix.size() * sizeof(T), cudaMemcpyHostToDevice);

        CHECK_CUDA_LAST_ERROR();

        CudaBuffer<T> d_buffer(numRows, this->numRows);
        const size_t bufferSizeByte = numRows * this->numRows * sizeof(T);

        CHECK_CUDA_LAST_ERROR();

        this->matrixMultiplication<T>(d_inputMatrix.get(),
                                      numCols,
                                      numRows,
                                      direction,
                                      d_outputMatrix.get(),
                                      outSize,
                                      false,
                                      d_buffer.get(),
                                      bufferSizeByte);
        this->wait();
        CHECK_CUDA_LAST_ERROR();

        std::vector<T> result(outSize);
        cudaMemcpy(result.data(), d_outputMatrix.get(), outSize * sizeof(T), cudaMemcpyDeviceToHost);
        CHECK_CUDA_LAST_ERROR();

        return std::move(result);
    }

    // This function is for performance testing.
    // Will run the same operation multiple times and report average time.
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

        CudaBuffer<T> d_inputMatrix(inputMatrix.size());
        CudaBuffer<T> d_outputMatrix(outSize);
        cudaMemcpy(d_inputMatrix.get(), inputMatrix.data(), inputMatrix.size() * sizeof(T), cudaMemcpyHostToDevice);

        CHECK_CUDA_LAST_ERROR();

        CudaBuffer<T> d_buffer(numRows, this->numRows);
        const size_t bufferSizeByte = numRows * this->numRows * sizeof(T);

        CHECK_CUDA_LAST_ERROR();

        auto start = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; it++) {
            this->matrixMultiplication<T>(d_inputMatrix.get(),
                                          numCols,
                                          numRows,
                                          direction,
                                          d_outputMatrix.get(),
                                          outSize,
                                          false,
                                          d_buffer.get(),
                                          bufferSizeByte);
            this->wait();
            CHECK_CUDA_LAST_ERROR();
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;

        std::cout << "Average end-to-end runtime for matrixMultiplication function (GPU) over " << iterations
                  << " iterations: " << (elapsed.count() / iterations) << " ms" << std::endl;

        std::vector<T> result(outSize);
        cudaMemcpy(result.data(), d_outputMatrix.get(), outSize * sizeof(T), cudaMemcpyDeviceToHost);
        CHECK_CUDA_LAST_ERROR();

        return std::move(result);
    }

    NodeIDSizeT* getRowOffsets() {
        if (!rowOffsets) {
            throw ApiMisuseFailure("Row offsets GPU buffer is not allocated");
        }
        return rowOffsets->get();
    }

    NodeIDSizeT* getOldToNewMapping() {
        if (!oldToNewMapping) {
            throw ApiMisuseFailure("Old to new mapping GPU buffer is not allocated");
        }
        return oldToNewMapping->get();
    }

    NodeIDSizeT* getNewToOldMapping() {
        if (!newToOldMapping) {
            throw ApiMisuseFailure("New to old mapping GPU buffer is not allocated");
        }
        return newToOldMapping->get();
    }

    NodeIDSizeT* getMutationAndNewMapping() {
        if (!mutationAndNewMapping) {
            throw ApiMisuseFailure("Mutation and new mapping GPU buffer is not allocated");
        }
        return mutationAndNewMapping->get();
    }

    NodeIDSizeT* getColIndices() {
        if (!colIndices) {
            throw ApiMisuseFailure("Col indices GPU buffer is not allocated");
        }
        return colIndices->get();
    }

    // return the size of GPU memory needed by col indices in bytes
    size_t queryColIndSize() { return numEdges * sizeof(NodeIDSizeT); }

    // This function sets existing GPU memory pointers for col indices
    // will not allocate new memory
    // the size should be in **BYTES**
    void setColIndBuffer(CudaNodeIDBufferPtr buffer, size_t size) {
        if (size < queryColIndSize()) {
            throw ApiMisuseFailure("Provided buffer size is smaller than required col indices size");
        }
        if (colIndices) {
            throw ApiMisuseFailure("Col indices buffer is already set");
        }
        colIndices = buffer;
    }

    // This function copies col indices from Host to GPU
    // Asynchronous copy is supported
    // The caller must ensure the memory space is not used by other operations.
    void copyColIndFromHost() {
        if (colIndices) {
            cudaMemcpyAsync(getColIndices(),
                            hostColIndices.data(),
                            numEdges * sizeof(NodeIDSizeT),
                            cudaMemcpyHostToDevice,
                            *workStreamPtr);
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

    // Free GPU memory
    void free() {
        rowOffsets = nullptr;
        colIndices = nullptr;
        oldToNewMapping = nullptr;
        newToOldMapping = nullptr;
        mutationAndNewMapping = nullptr;
        CHECK_CUDA_LAST_ERROR();

        hostColIndices.clear();
        hostHeightCutoffs.clear();
        hostHeavyCutoffs.clear();
        hostAvgChildCounts.clear();

        numRows = numSamples = numEdges = 0;
        numMutations = 0;
        maxHeight = 0;

        if (workStreamPtr) {
            // This may not work if we are multithreading
            if (workStreamPtr.use_count() == 1) {
                cudaStreamDestroy(*workStreamPtr);
            }
            workStreamPtr = nullptr;
        }
        CHECK_CUDA_LAST_ERROR();
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

    // This is the corresponding host backup data
    std::vector<NodeIDSizeT> hostColIndices;

    std::vector<NodeIDSizeT> hostHeightCutoffs;
    std::vector<NodeIDSizeT> hostHeavyCutoffs;
    std::vector<double> hostAvgChildCounts;

    size_t numRows;      // Number of rows, i.e. number of nodes
    size_t numSamples;   // Number of samples, i.e. number of leaf nodes
    size_t numMutations; // Number of mutations. may NOT correspond to number of non-leaf nodes
    size_t numEdges;     // Number of non-zero elements, i.e. number of edges
    size_t maxHeight;

    std::shared_ptr<cudaStream_t> workStreamPtr; // CUDA stream for asynchronous operations
private:
    // these are resident GPU pointers
    CudaNodeIDBufferPtr rowOffsets; // Device pointer to row offsets (size: num_rows + 1)
    CudaNodeIDBufferPtr oldToNewMapping;
    CudaNodeIDBufferPtr newToOldMapping;       // This is not used actually
    CudaNodeIDBufferPtr mutationAndNewMapping; // Mapping between mutation id and new node id (in pairs)

    // Since col indices may be large, this data may not be resident on GPU
    CudaNodeIDBufferPtr colIndices; // Device pointer to column indices (size: nnz)
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
 * - uint64_t index_size (8 bytes) - Size of NodeIDSizeT type in bytes
 *
 * Data sections:
 * - row_offsets: (num_rows + 1) * sizeof(NodeIDSizeT) bytes
 * - col_indices: nnz * sizeof(NodeIDSizeT) bytes
 * - height_cutoffs: (max_height + 1) * sizeof(NodeIDSizeT) bytes
 * - heavy_cutoffs: (max_height + 1) * sizeof(NodeIDSizeT) bytes
 * - avg_child_counts: (max_height + 1) * sizeof(double) bytes
 * - old_to_new_mapping: num_rows * sizeof(NodeIDSizeT) bytes
 * - new_to_old_mapping: num_rows * sizeof(NodeIDSizeT) bytes
 */

/**
 * Store a GPUGRG structure to disk in binary format.
 *
 * @param gpu_grg The GPUGRG structure to store (must have valid CPU data)
 * @param filename Path to the output file
 * @throws std::runtime_error if file operations fail
 */
void storeGPUGRGToDisk(GPUGRG& gpu_grg, const std::string& filename);

/**
 * Load a GPUGRG structure from disk.
 *
 * @param filename Path to the input file
 * @return GPUGRG structure loaded from disk with GPU memory allocated and data copied
 * @throws std::runtime_error if file operations fail or format is invalid
 */
GPUGRG loadGPUGRGFromDisk(const std::string& filename);

namespace {
struct nodes {
    NodeIDSizeT id;
    NodeIDSizeT childCounts;
};

bool cmp(const nodes& a, const nodes& b) { return a.childCounts > b.childCounts; }
} // namespace

// This visitor is tested ONLY on CSRGRG (immutable, ordered GRG)
// in the first visitor path, we get the height of each node, and thus the new id. new id is of (height, gid)
// we construct both old_id -> new_id map and a new_id -> old_id map
// then, in the epilogue function, we can use the information to have a new GPUGRG
class GPUCsrVisitor : public GRGVisitor {
public:
    explicit GPUCsrVisitor(const size_t num_nodes, const size_t num_elements) {
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

        NodeIDSizeT height = 0;

        if (nodeId < grg->numSamples()) {
            height = 0;
        } else {
            const auto downEdges = grg->getDownEdges(nodeId);
            for (const auto& child : downEdges) {
                release_assert(child < grg->numNodes());
                const auto childHeight = h_newID[child].first;
                height = std::max(height, childHeight + 1);
            }
        }

        if (height >= h_maxHeight) {
            h_maxHeight = height + 1;
            // m_idCounters.resize(m_maxHeight, 0);
            h_oldId.resize(h_maxHeight);
        }

        // auto id_in_group = m_idCounters[height]++;
        const auto idInGroup = h_oldId[height].size();
        h_oldId[height].push_back(nodeId);
        h_newID[nodeId] = std::make_pair(height, idInGroup);
        return true;
    }

    // Flattened mapping from new_id to old_id
    std::vector<NodeIDSizeT> getOldId() const {
        std::vector<NodeIDSizeT> flatOldID;
        for (const auto& vec : h_oldId) {
            flatOldID.insert(flatOldID.end(), vec.begin(), vec.end());
        }
        return std::move(flatOldID);
    }

    // Flattened mapping from old_id to new_id
    std::vector<NodeIDSizeT> getNewId() const {
        std::vector<NodeIDSizeT> flatNewID;
        flatNewID.resize(h_numNodes);
        auto baseID = 0;
        for (size_t height = 0; height < h_oldId.size(); height++) {
            for (size_t i = 0; i < h_oldId[height].size(); i++) {
                NodeID old_id = h_oldId[height][i];
                NodeIDSizeT new_id = baseID + i;
                flatNewID[old_id] = new_id;
            }
            baseID += h_oldId[height].size();
        }
        return std::move(flatNewID);
    }

    // this function rearranges nodes at each height by their child counts
    void rearrange(GRG* grg) {
        for (int h = 1; h < h_maxHeight; h++) {
#ifdef GRGL_GPU_DEBUG
            std::cout << "Rearranging height " << h << " with " << h_oldId[h].size() << " nodes." << std::endl;
#endif
            std::vector<nodes> nodeList;
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                const NodeID oldID = h_oldId[h][i];
                const NodeIDSizeT childCount = grg->numDownEdges(oldID);
                nodeList.push_back({oldID, childCount});
            }
            std::sort(nodeList.begin(), nodeList.end(), cmp);
            for (size_t i = 0; i < nodeList.size(); i++) {
                h_oldId[h][i] = nodeList[i].id;
                const NodeID old_id = nodeList[i].id;
                h_newID[old_id] = std::make_pair(h, i);
            }
        }
    }

    // After visiting all nodes, we can construct the CSR representation
    // This function should be called only once after all visits are done
    void constructGPUGRG(GRG* grg, GPUGRG& gpu_grg) {

        std::vector<NodeIDSizeT> h_heightCutoffs;
        h_heightCutoffs.resize(h_maxHeight + 1, 0);
        for (int i = 1; i <= h_maxHeight; i++) {
            h_heightCutoffs.at(i) = h_heightCutoffs.at(i - 1) + h_oldId.at(i - 1).size();
        }

        std::vector<NodeIDSizeT> h_heavyCutoffs;
        h_heavyCutoffs.resize(h_maxHeight + 1, 0);

        std::vector<NodeIDSizeT> h_rowOffsets;
        h_rowOffsets.resize(h_numNodes + 1, 0);

        std::vector<double> h_avgChildCounts;
        h_avgChildCounts.resize(h_maxHeight, 0.0);

        std::vector<NodeIDSizeT> h_colIndices;

        h_childCounts.resize(h_maxHeight, 0);

        NodeIDSizeT currrentColIndex = 0;
        for (int h = 0; h < h_maxHeight; h++) {
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                const NodeID oldID = h_oldId[h][i];
                const NodeIDSizeT newID = h_heightCutoffs[h] + i;
                h_rowOffsets[newID] = currrentColIndex;
                currrentColIndex += grg->numDownEdges(oldID);

                // for each down edge, we need to find its new id
                const auto downEdges = grg->getDownEdges(oldID);
                h_childCounts[h] += downEdges.size();
                for (const auto& child : downEdges) {
                    release_assert(child < grg->numNodes());
                    const auto childNewID = h_newID[child];
                    const NodeIDSizeT childNewFlatID = h_heightCutoffs[childNewID.first] + childNewID.second;
                    h_colIndices.push_back(childNewFlatID);
                }
            }
        }

#ifdef GRGL_GPU_DEBUG
        std::cout << "Host col indices constructed with size: " << h_colIndices.size() << std::endl;
#endif

        constexpr double HEAVY_RATIO = 8.0; // TODO: define this in another place

        for (int h = 1; h < h_maxHeight; h++) {
            const NodeIDSizeT heavyCutoffValue = HEAVY_RATIO * h_childCounts.at(h) / h_oldId.at(h).size();
            h_avgChildCounts.at(h) = static_cast<double>(h_childCounts.at(h)) / h_oldId.at(h).size();
            for (size_t i = 0; i < h_oldId[h].size(); i++) {
                const NodeID oldID = h_oldId[h][i];
                const NodeIDSizeT childCount = grg->numDownEdges(oldID);
                if (childCount < heavyCutoffValue) {
                    h_heavyCutoffs.at(h) = i;
                    break;
                }
            }
        }

#ifdef GRGL_GPU_DEBUG
        std::cout << "Host Heavy Cutoffs constructed." << std::endl;
#endif

        h_rowOffsets[h_numNodes] = currrentColIndex;
        release_assert(currrentColIndex == h_numEdges);

        const auto oldToNewFlat = getNewId();

        // calculate the mutation and new id pairs
        release_assert(sizeof(NodeIDSizeT) >= sizeof(MutationId));
        release_assert(sizeof(NodeIDSizeT) >= sizeof(NodeID));

        std::vector<NodeIDSizeT> mutationNewPairs;
        for (const auto& triple : grg->getNodesAndMutations<GRG::MutNodeMiss>()) {
            const NodeID& nodeId = std::get<0>(triple);
            const MutationId& mutId = std::get<1>(triple);

            if (nodeId == INVALID_NODE_ID) {
                mutationNewPairs.push_back((NodeIDSizeT)INVALID_NODE_ID);
                mutationNewPairs.push_back((NodeIDSizeT)INVALID_NODE_ID);
                continue;
            }
            release_assert(nodeId < h_numNodes);
            release_assert(mutId < grg->numMutations());
            const NodeIDSizeT& newID = oldToNewFlat[nodeId];
            mutationNewPairs.push_back((NodeIDSizeT)mutId);
            mutationNewPairs.push_back((NodeIDSizeT)newID);
        }

        const size_t numValidMutations = mutationNewPairs.size() / 2;

#ifdef GRGL_GPU_DEBUG
        std::cout << "Mutation and new id pairs constructed with size: " << mutationNewPairs.size() << std::endl;
#endif
        if (numValidMutations != grg->numMutations()) {
            std::cout << "Warning: Valid mutations (" << numValidMutations << ") is less than total mutations ("
                      << grg->numMutations() << ")." << std::endl;
        }
#ifdef GRGL_GPU_DEBUG
        std::cout << "Constructing GPUGRG with " << h_numNodes << " nodes, " << h_numEdges << " elements, max height "
                  << h_maxHeight << std::endl;
#endif
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
    NodeIDSizeT h_maxHeight;
    size_t h_numNodes;
    size_t h_numEdges;

    std::vector<size_t> h_childCounts;                        // sum of child counts for nodes of the same height
    std::vector<std::vector<NodeIDSizeT>> h_oldId;            // new_id -> old_id mapping
    std::vector<std::pair<NodeIDSizeT, NodeIDSizeT>> h_newID; // old_id -> new_id mapping
};

/**
 * Convert a GRG to GPUGRG format for GPU matrix multiplication.
 *
 * @param grg The GRG to convert (must have ordered nodes)
 * @return GPUGRG structure with GPU memory allocated and data copied
 * @throws std::runtime_error if GRG doesn't meet requirements
 */
GPUGRG convertGRGToGPUGRG(GRGPtr grg);

bool hasCudaSupport();

} // namespace grgl

#endif // GRGL_GPU_GRG_H
