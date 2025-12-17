// #define GRGL_CUDA_ENABLED
#ifndef GRGL_CUDA_MATMUL_CU
#define GRGL_CUDA_MATMUL_CU

#include "grgl/common.h"
#ifndef GRGL_CUDA_ENABLED
#error "Cannot build cuda_matmul.cu if CUDA is not enabled; rebuild with GRGL_CUDA_ENABLED defined."
#endif

#include <chrono>
#include <cstdint>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

#include "grgl/cuda/cuda_kernels.cuh"
#include "grgl/cuda/gpu_grg.h"
#include "grgl/grgnode.h"

namespace grgl {

#define MIN_BLOCK_SIZE           64  // set to 64 for better occupancy
#define CONCURRENT_BLOCKS_PER_SM 32  // 2048 / MIN_BLOCK_SIZE
#define NUM_SMS                  108 // for A100 GPU
#define MIN_LUNCHED_BLOCKS       (CONCURRENT_BLOCKS_PER_SM * NUM_SMS)

template <typename data_t>
void cudaTraverseUpSingleElementLauncher(NodeIDSizeT* rowOffsets,
                                         NodeIDSizeT* colIndices,
                                         data_t* values,
                                         size_t rowStart,
                                         size_t rowEnd,
                                         size_t heavyCutoff,
                                         double avgEdgePerNode,
                                         cudaStream_t stream) {

    if (heavyCutoff != 0) {
        const int blockSize = MIN_BLOCK_SIZE;
        const int rowsPerBlock = 2;
        const int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
        cudaTraversalUpSingleElementKernel<data_t, 32><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowStart + heavyCutoff, rowsPerBlock);
        rowStart += heavyCutoff;
    }

    const int blockSize = MIN_BLOCK_SIZE;
    const int numRows = rowEnd - rowStart;
    const int rowsPerBlockMaxOccupancy = (numRows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rowsPerBlock = 32;
    if (rowsPerBlockMaxOccupancy < rowsPerBlock) {
        rowsPerBlock = rowsPerBlockMaxOccupancy;
    }
    const int numBlocks = (numRows + rowsPerBlock - 1) / rowsPerBlock;

    if (avgEdgePerNode <= 8) {
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 16) {
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 32) {
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    }
}

template <typename data_t>
void cudaTraverseDownSingleElementLauncher(NodeIDSizeT* rowOffsets,
                                           NodeIDSizeT* colIndices,
                                           data_t* values,
                                           size_t rowStart,
                                           size_t rowEnd,
                                           size_t heavyCutoff,
                                           double avgEdgePerNode,
                                           cudaStream_t stream) {

    if (heavyCutoff != 0) {
        const int blockSize = MIN_BLOCK_SIZE;
        const int rowsPerBlock = 1;
        const int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
        cudaTraversalDownSingleElementKernel<data_t, 32><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowStart + heavyCutoff, rowsPerBlock);
        rowStart += heavyCutoff;
    }

    const int blockSize = MIN_BLOCK_SIZE;
    const int numRows = rowEnd - rowStart;
    const int rowsPerBlockMaxOccupancy = (numRows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rowsPerBlock = 16;
    if (rowsPerBlockMaxOccupancy < rowsPerBlock) {
        rowsPerBlock = rowsPerBlockMaxOccupancy;
    }

    const int numBlocks = (numRows + rowsPerBlock - 1) / rowsPerBlock;

    if (avgEdgePerNode <= 8) {
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 16) {
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 32) {
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW>
            <<<numBlocks, blockSize, 0, stream>>>(rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    }
}

template <typename data_t>
void cudaTraverseUP(GPUGRG& gpuGRG, data_t* d_inoutValues, cudaStream_t& stream, NodeIDSizeT unit = 1) {

#ifdef GRGL_GPU_DEBUG
    std::cout << "Starting GPU traversal..." << std::endl;
#endif

    // launch the kernels level by level
    for (size_t level = 1; level < gpuGRG.maxHeight(); level++) {
        const size_t rowStart = gpuGRG.getHostHeightCutoff(level);
        const size_t rowEnd = gpuGRG.getHostHeightCutoff(level + 1);
        const NodeIDSizeT heavyCutoff = gpuGRG.getHostHeavyCutoff(level);

        if (unit == 1) {
            cudaTraverseUpSingleElementLauncher<data_t>(gpuGRG.getRowOffsets(),
                                                        gpuGRG.getColIndices(),
                                                        d_inoutValues,
                                                        rowStart,
                                                        rowEnd,
                                                        heavyCutoff,
                                                        gpuGRG.getHostAvgChildCount(level),
                                                        stream);

        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;
            if (heavyCutoff != 0) {
                const int blockSize = HEAVY_THREADS_PER_UNIT * unit;
                const int rowsPerBlock = 1;
                const int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

                cudaTraversalUpMultiElementKernel<data_t>
                    <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * HEAVY_THREADS_PER_UNIT, stream>>>(
                        gpuGRG.getRowOffsets(),
                        gpuGRG.getColIndices(),
                        d_inoutValues,
                        rowStart,
                        rowStart + heavyCutoff,
                        rowsPerBlock,
                        unit);
            }

            const int blockSize = unit * THREADS_PER_UNIT;
            const int rowsPerBlock = 1;
            const int numBlocks = (rowEnd - rowStart - heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

            cudaTraversalUpMultiElementKernel<data_t>
                <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * THREADS_PER_UNIT, stream>>>(
                    gpuGRG.getRowOffsets(),
                    gpuGRG.getColIndices(),
                    d_inoutValues,
                    rowStart + heavyCutoff,
                    rowEnd,
                    rowsPerBlock,
                    unit);
        }
    }
}

template <typename data_t>
void cudaTraverseDOWN(GPUGRG& gpuGRG, data_t* d_inoutValues, cudaStream_t stream, NodeIDSizeT unit = 1) {

#ifdef GRGL_GPU_DEBUG
    std::cout << "Starting GPU traversal..." << std::endl;
#endif

    // launch the kernels level by level
    for (int level = gpuGRG.maxHeight() - 1; level > 0; level--) {
        const size_t rowStart = gpuGRG.getHostHeightCutoff(level);
        const size_t rowEnd = gpuGRG.getHostHeightCutoff(level + 1);
        const double avgEdgePerNode = gpuGRG.getHostAvgChildCount(level);
        const size_t heavyCutoff = gpuGRG.getHostHeavyCutoff(level);

        if (unit == 1) {
            cudaTraverseDownSingleElementLauncher<data_t>(gpuGRG.getRowOffsets(),
                                                          gpuGRG.getColIndices(),
                                                          d_inoutValues,
                                                          rowStart,
                                                          rowEnd,
                                                          heavyCutoff,
                                                          avgEdgePerNode,
                                                          stream);
        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;

            if (heavyCutoff != 0) {
                const int blockSize = HEAVY_THREADS_PER_UNIT * unit;
                const int rowsPerBlock = 1;
                const int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
                cudaTraversalDownMultiElementKernel<data_t>
                    <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * HEAVY_THREADS_PER_UNIT, stream>>>(
                        gpuGRG.getRowOffsets(),
                        gpuGRG.getColIndices(),
                        d_inoutValues,
                        rowStart,
                        rowStart + heavyCutoff,
                        rowsPerBlock,
                        unit);
            }

            const int blockSize = unit * THREADS_PER_UNIT;
            const int rowsPerBlock = 1;
            const int numBlocks = (rowEnd - rowStart - heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

            cudaTraversalDownMultiElementKernel<data_t>
                <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * THREADS_PER_UNIT, stream>>>(
                    gpuGRG.getRowOffsets(),
                    gpuGRG.getColIndices(),
                    d_inoutValues,
                    rowStart + heavyCutoff,
                    rowEnd,
                    rowsPerBlock,
                    unit);
        }
    }
}

template <typename data_t>
void GPUGRG::matrixMultiplication(const data_t* inputMatrix,
                                  size_t inputCols,
                                  size_t inputRows,
                                  TraversalDirection direction,
                                  data_t* outputMatrix,
                                  size_t outputSize,
                                  bool emitAllNodes,
                                  bool byIndividual,
                                  data_t* buffer,
                                  size_t buffer_size) {

    // make sure IOType and NodeValueType are the same type
    // static_assert(std::is_same<IOType, NodeValueType>::value, "IOType and NodeValueType must be the same type for GPU
    // matmul."); make data_t an alias for IOType/NodeValueType
    api_exc_check(!byIndividual, "by_individual is not yet supported: TODO");

    release_assert(inputCols > 0);
    release_assert(inputRows > 0);
    const size_t outputCols = outputSize / inputRows;
    validateMatMulInputs<GPUGRG>(
        this, inputCols, inputRows, direction, outputSize, emitAllNodes, byIndividual, outputCols);

    release_assert(buffer_size >= this->numNodes() * inputRows *
                                      sizeof(data_t)); // buffer should be large enough to hold all node values
    cudaMemsetAsync(buffer, 0, this->numNodes() * inputRows * sizeof(data_t), *(this->workStreamPtr));

    const auto numFeatures = inputRows;

    if (direction == DIRECTION_DOWN) {

        // Downward, we are calculating "how do the mutations impact the samples?"

        const int nodePerBlock = 128 / numFeatures;
        const int blockSize = nodePerBlock * numFeatures;
        int numBlocks = (this->numMutations() + nodePerBlock - 1) / nodePerBlock;

        cudaReorderPermutationPairKernel<data_t, false, true, true>
            <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(buffer,
                                                                  inputMatrix,
                                                                  this->getMutationAndNewMapping(),
                                                                  0,
                                                                  this->numMutations(),
                                                                  numFeatures,
                                                                  this->numMutations(),
                                                                  nodePerBlock);
        CHECK_CUDA_LAST_ERROR();

        cudaTraverseDOWN<data_t>(*this, buffer, *(this->workStreamPtr), numFeatures);
        CHECK_CUDA_LAST_ERROR();

        if (emitAllNodes) {
            numBlocks = (this->numNodes() + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                      reinterpret_cast<data_t*>(buffer),
                                                                      this->getOldToNewMapping(),
                                                                      0,
                                                                      this->numNodes(),
                                                                      numFeatures,
                                                                      this->numNodes(),
                                                                      nodePerBlock);
        } else {
            cudaMemsetAsync(outputMatrix, 0, this->numSamples() * numFeatures * sizeof(data_t), *(this->workStreamPtr));
            numBlocks = (this->numSamples() + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                      reinterpret_cast<data_t*>(buffer),
                                                                      this->getOldToNewMapping(),
                                                                      0,
                                                                      this->numSamples(),
                                                                      numFeatures,
                                                                      this->numSamples(),
                                                                      nodePerBlock);
        }
        CHECK_CUDA_LAST_ERROR();

    } else {
        // Upward, we are calculating "how do the samples impact the mutations?"
        const int nodePerBlock = 128 / numFeatures;
        const int blockSize = nodePerBlock * numFeatures;
        int numBlocks = (this->numSamples() + nodePerBlock - 1) / nodePerBlock;

        cudaReorderMapKernel<data_t, false, true>
            <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(reinterpret_cast<data_t*>(buffer),
                                                                  inputMatrix,
                                                                  this->getOldToNewMapping(),
                                                                  0,
                                                                  this->numSamples(),
                                                                  numFeatures,
                                                                  this->numSamples(),
                                                                  nodePerBlock);

        CHECK_CUDA_LAST_ERROR();

        cudaTraverseUP<data_t>(*this, buffer, *(this->workStreamPtr), numFeatures);
        CHECK_CUDA_LAST_ERROR();

        // Since cudaTraverseUP is currently blocking, we do not need to do stream synchronization here.

        if (emitAllNodes) {
            numBlocks = (this->numNodes() + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                      reinterpret_cast<data_t*>(buffer),
                                                                      this->getOldToNewMapping(),
                                                                      0,
                                                                      this->numNodes(),
                                                                      numFeatures,
                                                                      this->numNodes(),
                                                                      nodePerBlock);
        } else {
            cudaMemsetAsync(
                outputMatrix, 0, this->numMutations() * numFeatures * sizeof(data_t), *(this->workStreamPtr));
            numBlocks = (this->numMutations() + nodePerBlock - 1) / nodePerBlock;
            cudaReorderPermutationPairKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                      reinterpret_cast<data_t*>(buffer),
                                                                      this->getMutationAndNewMapping(),
                                                                      0,
                                                                      this->numMutations(),
                                                                      numFeatures,
                                                                      this->numMutations(),
                                                                      nodePerBlock);
        }
        CHECK_CUDA_LAST_ERROR();
    }
    CHECK_CUDA_LAST_ERROR();
}

#define MATMUL_TEMPLATE(T)                                                                                             \
    template void GPUGRG::matrixMultiplication<T>(                                                                     \
        const T*, size_t, size_t, TraversalDirection, T*, size_t, bool, bool, T*, size_t)

MATMUL_TEMPLATE(float);
MATMUL_TEMPLATE(double);
MATMUL_TEMPLATE(uint32_t);
MATMUL_TEMPLATE(int32_t);

// Current CUDA implementation does not support these types for multiplication, because there is no
// atomic addition operation for them.
// MATMUL_TEMPLATE(uint64_t);
// MATMUL_TEMPLATE(int64_t);
// MATMUL_TEMPLATE(uint8_t);
// MATMUL_TEMPLATE(uint16_t);
// MATMUL_TEMPLATE(int8_t);
// MATMUL_TEMPLATE(int16_t);

bool hasCudaSupport() {
    int deviceCount;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    return (error == cudaSuccess && deviceCount > 0);
}

} // namespace grgl

#endif // GRGL_CUDA_MATMUL_CU
