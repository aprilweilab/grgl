// #define GRGL_CUDA_ENABLED
#ifndef GRGL_CUDA_MATMUL_CU
#define GRGL_CUDA_MATMUL_CU

#ifndef GRGL_CUDA_ENABLED
#error "Cannot build cuda_matmul.cu if CUDA is not enabled; rebuild with GRGL_CUDA_ENABLED defined."
#endif

#include <chrono>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
// #include "grgl/grg.h"
// #include "grgl/cuda_matmul.h"
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
        // std::cout << "Warning: heavyCutoff is not zero in cudaTraverseUpSingleElementLauncher. This may lead to
        // performance issues." << std::endl;
        int blockSize = MIN_BLOCK_SIZE;
        int rowsPerBlock = 2;
        int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
        cudaTraversalUpSingleElementKernel<data_t, 32><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowStart + heavyCutoff, rowsPerBlock);
        rowStart += heavyCutoff;
    }

    int blockSize = MIN_BLOCK_SIZE;
    int numRows = rowEnd - rowStart;
    int rowsPerBlockMaxOccupancy = (numRows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rowsPerBlock = 32;
    if (rowsPerBlockMaxOccupancy < rowsPerBlock) {
        rowsPerBlock = rowsPerBlockMaxOccupancy;
    }
    int numBlocks = (numRows + rowsPerBlock - 1) / rowsPerBlock;

    /*
    std::cout << "Launching cudaTraverseUpSingleElementLauncher with blockSize " << blockSize
              << ", rowsPerBlock " << rowsPerBlock
              << ", numBlocks " << numBlocks
              << ", avgEdgePerNode " << avgEdgePerNode
              << std::endl;
    */

    if (avgEdgePerNode <= 8) {
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 16) {
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 32) {
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalUpSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
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
        // std::cout << "Warning: heavyCutoff is not zero in cudaTraverseDownSingleElementLauncher. This may lead to
        // performance issues." << std::endl;
        int blockSize = MIN_BLOCK_SIZE;
        int rowsPerBlock = 1;
        int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
        cudaTraversalDownSingleElementKernel<data_t, 32><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowStart + heavyCutoff, rowsPerBlock);
        rowStart += heavyCutoff;
    }

    int blockSize = MIN_BLOCK_SIZE;
    int numRows = rowEnd - rowStart;
    int rowsPerBlockMaxOccupancy = (numRows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rowsPerBlock = 16;
    if (rowsPerBlockMaxOccupancy < rowsPerBlock) {
        rowsPerBlock = rowsPerBlockMaxOccupancy;
    }

    int numBlocks = (numRows + rowsPerBlock - 1) / rowsPerBlock;

    // std::cout << "Launching cudaTraverseDownSingleElementLauncher with " << numBlocks << " blocks, block size " <<
    // blockSize << ", rows per block " << rowsPerBlock << ", avgEdgePerNode " << avgEdgePerNode << std::endl;

    if (avgEdgePerNode <= 8) {
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 16) {
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else if (avgEdgePerNode <= 32) {
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalDownSingleElementKernel<data_t, THREADS_PER_ROW><<<numBlocks, blockSize, 0, stream>>>(
            rowOffsets, colIndices, values, rowStart, rowEnd, rowsPerBlock);
    }
}

template <typename data_t>
void cudaTraverseUP(GPUGRG& gpuGRG, data_t* d_inoutValues, cudaStream_t& stream, NodeIDSizeT unit = 1) {
    std::cout << "Starting GPU traversal..." << std::endl;

    // launch the kernels level by level
    for (size_t level = 1; level < gpuGRG.maxHeight; level++) {
        size_t rowStart = gpuGRG.hostHeightCutoffs[level];
        size_t rowEnd = gpuGRG.hostHeightCutoffs[level + 1];
        // size_t numRows = rowEnd - rowStart;
        NodeIDSizeT heavyCutoff = gpuGRG.hostHeavyCutoffs[level];

        if (unit == 1) {
            cudaTraverseUpSingleElementLauncher<data_t>(gpuGRG.rowOffsets,
                                                                 gpuGRG.colIndices,
                                                                 d_inoutValues,
                                                                 rowStart,
                                                                 rowEnd,
                                                                 heavyCutoff,
                                                                 gpuGRG.hostAvgChildCounts[level],
                                                                 stream);

        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;
            if (heavyCutoff != 0) {
                int blockSize = HEAVY_THREADS_PER_UNIT * unit;
                int rowsPerBlock = 1;
                int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

                cudaTraversalUpMultiElementKernel<data_t>
                    <<<numBlocks,
                       blockSize,
                       COL_BUFFER_ITER * sizeof(NodeIDSizeT) * HEAVY_THREADS_PER_UNIT,
                       stream>>>(gpuGRG.rowOffsets,
                                 gpuGRG.colIndices,
                                 d_inoutValues,
                                 rowStart,
                                 rowStart + heavyCutoff,
                                 rowsPerBlock,
                                 unit);
            }

            int blockSize = unit * THREADS_PER_UNIT;
            int rowsPerBlock = 1;
            int numBlocks = (rowEnd - rowStart - heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

            cudaTraversalUpMultiElementKernel<data_t>
                <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * THREADS_PER_UNIT, stream>>>(
                    gpuGRG.rowOffsets,
                    gpuGRG.colIndices,
                    d_inoutValues,
                    rowStart + heavyCutoff,
                    rowEnd,
                    rowsPerBlock,
                    unit);
        }
        // std::cout << "Launched kernel for level " << level << " with rowStart " << rowStart << " and rowEnd " <<
        // rowEnd << std::endl;
    }
}

template <typename data_t>
void cudaTraverseDOWN(GPUGRG& gpuGRG, data_t* d_inoutValues, cudaStream_t stream, NodeIDSizeT unit = 1) {
    std::cout << "Starting GPU traversal..." << std::endl;

    // launch the kernels level by level
    for (int level = gpuGRG.maxHeight - 1; level > 0; level--) {
        size_t rowStart = gpuGRG.hostHeightCutoffs[level];
        size_t rowEnd = gpuGRG.hostHeightCutoffs[level + 1];
        double avgEdgePerNode = gpuGRG.hostAvgChildCounts[level];
        // size_t numRows = rowEnd - rowStart;
        size_t heavyCutoff = gpuGRG.hostHeavyCutoffs[level];

        if (unit == 1) {
            cudaTraverseDownSingleElementLauncher<data_t>(gpuGRG.rowOffsets,
                                                                   gpuGRG.colIndices,
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
                int blockSize = HEAVY_THREADS_PER_UNIT * unit;
                int rowsPerBlock = 1;
                int numBlocks = (heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;
                cudaTraversalDownMultiElementKernel<data_t>
                    <<<numBlocks,
                       blockSize,
                       COL_BUFFER_ITER * sizeof(NodeIDSizeT) * HEAVY_THREADS_PER_UNIT,
                       stream>>>(gpuGRG.rowOffsets,
                                 gpuGRG.colIndices,
                                 d_inoutValues,
                                 rowStart,
                                 rowStart + heavyCutoff,
                                 rowsPerBlock,
                                 unit);
            }

            int blockSize = unit * THREADS_PER_UNIT;
            int rowsPerBlock = 1;
            int numBlocks = (rowEnd - rowStart - heavyCutoff + rowsPerBlock - 1) / rowsPerBlock;

            cudaTraversalDownMultiElementKernel<data_t>
                <<<numBlocks, blockSize, COL_BUFFER_ITER * sizeof(NodeIDSizeT) * THREADS_PER_UNIT, stream>>>(
                    gpuGRG.rowOffsets,
                    gpuGRG.colIndices,
                    d_inoutValues,
                    rowStart + heavyCutoff,
                    rowEnd,
                    rowsPerBlock,
                    unit);
        }
        // std::cout << "Launched kernel for level " << level << " with rowStart " << rowStart << " and rowEnd " <<
        // rowEnd << std::endl;
    }

    std::cout << "Finished GPU traversal DOWN." << std::endl;
}

template <typename data_t>
void GPUGRG::matrixMultiplication(const data_t* inputMatrix,
                                           size_t inputCols,
                                           size_t inputRows,
                                           TraversalDirection direction,
                                           data_t* outputMatrix,
                                           size_t outputSize,
                                           bool emitAllNodes,
                                           data_t* buffer,
                                           size_t buffer_size) {

    // make sure IOType and NodeValueType are the same type
    // static_assert(std::is_same<IOType, NodeValueType>::value, "IOType and NodeValueType must be the same type for GPU
    // matmul."); make data_t an alias for IOType/NodeValueType

    // using data_t = T;
    // using NodeIDSizeT = uint32_t;

    release_assert(inputCols > 0);
    release_assert(inputRows > 0);
    release_assert(buffer_size >= this->numRows * inputRows *
                                      sizeof(data_t)); // buffer should be large enough to hold all node values
    cudaMemsetAsync(buffer, 0, this->numRows * inputRows * sizeof(data_t), *(this->workStreamPtr));
    // release_assert(useBitVector == false); // Bitvector not supported on GPU yet
    // cols are num of samples or mutations
    // rows are num of different input vectors per sample/mutation
    const size_t outputCols = outputSize / inputRows;
    release_assert(outputSize % inputRows == 0);
    release_assert(outputCols > 0);
    // other validations are skipped

    // create cuda streams
    /*
    const int num_streams = 4;
    std::vector<cudaStream_t> streams(num_streams);
    for (int i = 0; i < num_streams; i++) {
        cudaStreamCreate(&streams[i]);
    }
    */
    auto numFeatures = inputRows;

    if (direction == DIRECTION_DOWN) {

        // Downward, we are calculating "how do the mutations impact the samples?"

        int nodePerBlock = 128 / numFeatures;
        int blockSize = nodePerBlock * numFeatures;
        int numBlocks = (this->numMutations + nodePerBlock - 1) / nodePerBlock;

        cudaReorderPermutationPairKernel<data_t, false, true, true>
            <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(buffer,
                                                                             inputMatrix,
                                                                             this->mutationAndNewMapping,
                                                                             0,
                                                                             this->numMutations,
                                                                             numFeatures,
                                                                             this->numMutations,
                                                                             nodePerBlock);
        CHECK_CUDA_LAST_ERROR();

        cudaTraverseDOWN<data_t>(*this, buffer, *(this->workStreamPtr), numFeatures);
        CHECK_CUDA_LAST_ERROR();

        if (emitAllNodes) {
            numBlocks = (this->numRows + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                                 reinterpret_cast<data_t*>(buffer),
                                                                                 this->oldToNewMapping,
                                                                                 0,
                                                                                 this->numRows,
                                                                                 numFeatures,
                                                                                 this->numRows,
                                                                                 nodePerBlock);
        } else {
            cudaMemsetAsync(
                outputMatrix, 0, this->numSamples * numFeatures * sizeof(data_t), *(this->workStreamPtr));
            numBlocks = (this->numSamples + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                                 reinterpret_cast<data_t*>(buffer),
                                                                                 this->oldToNewMapping,
                                                                                 0,
                                                                                 this->numSamples,
                                                                                 numFeatures,
                                                                                 this->numSamples,
                                                                                 nodePerBlock);
        }
        CHECK_CUDA_LAST_ERROR();

    } else {
        // Upward, we are calculating "how do the samples impact the mutations?"
        int nodePerBlock = 128 / numFeatures;
        int blockSize = nodePerBlock * numFeatures;
        int numBlocks = (this->numSamples + nodePerBlock - 1) / nodePerBlock;

        cudaReorderMapKernel<data_t, false, true>
            <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(reinterpret_cast<data_t*>(buffer),
                                                                             inputMatrix,
                                                                             this->oldToNewMapping,
                                                                             0,
                                                                             this->numSamples,
                                                                             numFeatures,
                                                                             this->numSamples,
                                                                             nodePerBlock);

        CHECK_CUDA_LAST_ERROR();

        /*
        cudaEvent_t event;
        cudaEventCreate(&event);
        cudaEventRecord(event, *(this->work_stream_ptr));
        // techinically we need to destroy the events
        // however since it may cause sync problems they're not destroyed for now

        std::vector<cudaStream_t> streams;
        for (size_t i = 0; i < 2; i++) {
            cudaStream_t stream;
            cudaStreamCreate(&stream);
            streams.push_back(stream);
            cudaStreamWaitEvent(stream, event, 0);

        }

        CHECK_CUDA_LAST_ERROR();
        */

        cudaTraverseUP<data_t>(*this, buffer, *(this->workStreamPtr), numFeatures);
        CHECK_CUDA_LAST_ERROR();

        // Since cudaTraverseUP is currently blocking, we do not need to do stream synchronization here.
        /*
        std::vector<cudaEvent_t> end_events;
        for (size_t i = 0; i < streams.size(); i++) {
            cudaEvent_t end_event;
            cudaEventCreate(&end_event);
            cudaEventRecord(end_event, streams[i]);
            end_events.push_back(end_event);
            cudaStreamWaitEvent(*(this->work_stream_ptr), end_event, 0);
        }
        */

        if (emitAllNodes) {
            numBlocks = (this->numRows + nodePerBlock - 1) / nodePerBlock;
            cudaReorderMapKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                                 reinterpret_cast<data_t*>(buffer),
                                                                                 this->oldToNewMapping,
                                                                                 0,
                                                                                 this->numRows,
                                                                                 numFeatures,
                                                                                 this->numRows,
                                                                                 nodePerBlock);
        } else {
            cudaMemsetAsync(
                outputMatrix, 0, this->numMutations * numFeatures * sizeof(data_t), *(this->workStreamPtr));
            numBlocks = (this->numMutations + nodePerBlock - 1) / nodePerBlock;
            cudaReorderPermutationPairKernel<data_t, true, false>
                <<<numBlocks, blockSize, 0, *(this->workStreamPtr)>>>(outputMatrix,
                                                                                 reinterpret_cast<data_t*>(buffer),
                                                                                 this->mutationAndNewMapping,
                                                                                 0,
                                                                                 this->numMutations,
                                                                                 numFeatures,
                                                                                 this->numMutations,
                                                                                 nodePerBlock);
        }
        CHECK_CUDA_LAST_ERROR();
    }
    CHECK_CUDA_LAST_ERROR();
}

/*
template void GRG::matrixMultiplicationGPU<float, float, false>(
    const float*, size_t, size_t, TraversalDirection, float*, size_t,
    bool, bool, const float*, NodeInitEnum, float*);
*/
/*
template void GRG::matrixMultiplicationGPU<double, double, false>(
    const double*, size_t, size_t, TraversalDirection, double*, size_t,
    bool, bool, const double*, NodeInitEnum, double*);

template void GRG::matrixMultiplicationGPU<int, int, false>(
    const int*, size_t, size_t, TraversalDirection, int*, size_t,
    bool, bool, const int*, NodeInitEnum, int*);

template void GRG::matrixMultiplicationGPU<uint32_t, uint32_t, false>(
    const uint32_t*, size_t, size_t, TraversalDirection, uint32_t*, size_t,
    bool, bool, const uint32_t*, NodeInitEnum, uint32_t*);
    */
template void GPUGRG::matrixMultiplication<float>(
    const float*, size_t, size_t, TraversalDirection, float*, size_t, bool, float*, size_t);

template void GPUGRG::matrixMultiplication<double>(
    const double*, size_t, size_t, TraversalDirection, double*, size_t, bool, double*, size_t);

template void GPUGRG::matrixMultiplication<int>(
    const int*, size_t, size_t, TraversalDirection, int*, size_t, bool, int*, size_t);

/*
template void GPUGRG<uint32_t>::matrixMultiplication<long long>(
    const long long*, size_t, size_t, TraversalDirection,
    long long*, size_t, bool, long long*, size_t
);
*/

bool hasCudaSupport() {
    int deviceCount;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    return (error == cudaSuccess && deviceCount > 0);
}

/*
template void GPUGRG<uint64_t>::matrixMultiplication<long long>(
    const long long*, size_t, size_t, TraversalDirection,
    long long*, size_t, bool, long long*, size_t
);
*/

} // namespace grgl

#endif // GRGL_CUDA_MATMUL_CU
