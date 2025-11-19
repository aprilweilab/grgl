#ifndef GRGL_CUDA_KERNELS_CUH
#define GRGL_CUDA_KERNELS_CUH

#ifndef GRGL_CUDA_ENABLED
#error "Cannot include cuda_kernels.cuh if CUDA is not enabled; rebuild with GRGL_CUDA_ENABLED defined."
#endif

#include <chrono>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include "grgnode.h"

namespace grgl {

#define COL_BUFFER_ITER 16

template <class index_t, class data_t, bool INPUT_NODE_MAJOR, bool PERMUTATION_IS_SRC_TO_DST, bool ATOMIC_ADD = false>
__global__ void cudaReorderMapKernel(data_t* dst,
                                     const data_t* src,
                                     const index_t* permutation,
                                     size_t st,
                                     size_t ed,
                                     size_t numUnits,
                                     size_t numNodes,
                                     index_t nodePerBlock) {
    static_assert(INPUT_NODE_MAJOR == !PERMUTATION_IS_SRC_TO_DST,
                  "For multi-element reorder, INPUT_NODE_MAJOR must be opposite to PERMUTATION_IS_SRC_TO_DST");
    index_t myNode = st + blockIdx.x * nodePerBlock + threadIdx.x / numUnits;
    if (myNode >= ed) {
        return;
    }
    index_t myElementPos = threadIdx.x % numUnits;

    index_t srcID;
    index_t dstID;

    if (INPUT_NODE_MAJOR) {
        srcID = permutation[myNode] * numUnits + myElementPos;
        dstID = myNode + myElementPos * numNodes;
    } else {
        srcID = myElementPos * numNodes + myNode;
        dstID = permutation[myNode] * numUnits + myElementPos;
    }

    if (!ATOMIC_ADD)
        dst[dstID] = src[srcID];
    else
        atomicAdd(&dst[dstID], src[srcID]);
}

template <class index_t, class data_t, bool INPUT_NODE_MAJOR, bool PERMUTATION_FWD_ORDER, bool ATOMIC_ADD = false>
__global__ void cudaReorderPermutationPairKernel(data_t* dst,
                                                 const data_t* src,
                                                 const index_t* permutationPairs,
                                                 size_t st,
                                                 size_t ed,
                                                 size_t numUnits,
                                                 size_t numNodes,
                                                 index_t nodePerBlock) {
    static_assert(INPUT_NODE_MAJOR == !PERMUTATION_FWD_ORDER,
                  "For reorder, INPUT_NODE_MAJOR must be opposite to PERMUTATION_FWD_ORDER");
    index_t myNode = st + blockIdx.x * nodePerBlock + threadIdx.x / numUnits;
    if (myNode >= ed) {
        return;
    }
    index_t myElementPos = threadIdx.x % numUnits;

    index_t srcID;
    index_t dstID;

    index_t ele1 = permutationPairs[myNode * 2];
    index_t ele2 = permutationPairs[myNode * 2 + 1];

    if (ele1 == INVALID_NODE_ID || ele2 == INVALID_NODE_ID) {
        // invalid mapping, skip
        return;
    }

    if (INPUT_NODE_MAJOR) {
        srcID = ele2 * numUnits + myElementPos;
        dstID = ele1 + myElementPos * numNodes;
    } else {
        srcID = ele1 + myElementPos * numNodes;
        dstID = ele2 * numUnits + myElementPos;
    }

    if (!ATOMIC_ADD) {
        dst[dstID] = src[srcID];
    } else {
        atomicAdd(&dst[dstID], src[srcID]);
    }
}

template <class index_t, class data_t, int THREADS_PER_ROW>
__global__ void cudaTraversalUpSingleElementKernel(index_t* rowOffsets,
                                                   index_t* colIndices,
                                                   data_t* values,
                                                   size_t rowStart,
                                                   size_t rowEnd,
                                                   size_t rowsPerBlock) {
    index_t myStartRow = rowStart + (blockIdx.x * rowsPerBlock);
    index_t myEndRow = myStartRow + rowsPerBlock;
    if (myEndRow > rowEnd) {
        myEndRow = rowEnd;
    }
    int myGroupID = threadIdx.x / THREADS_PER_ROW;
    int myThreadID = threadIdx.x % THREADS_PER_ROW;
    int stride = blockDim.x / THREADS_PER_ROW;
    for (index_t row = myStartRow + myGroupID; row < myEndRow; row += stride) {
        index_t curStartRow = rowOffsets[row];
        index_t curEndRow = rowOffsets[row + 1];
        data_t sum = 0;
        for (index_t idx = curStartRow + myThreadID; idx < curEndRow; idx += THREADS_PER_ROW) {
            index_t col = colIndices[idx];
            sum += values[col];
        }
        // Reduce within the group. Use Warp shuffle.
        for (int offset = THREADS_PER_ROW / 2; offset > 0; offset /= 2) {
            sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
        }
        if (myThreadID == 0) {
            values[row] = sum;
        }
    }
}

template <class index_t, class data_t, int THREADS_PER_ROW>
__global__ void cudaTraversalDownSingleElementKernel(index_t* rowOffsets,
                                                     index_t* colIndices,
                                                     data_t* values,
                                                     size_t rowStart,
                                                     size_t rowEnd,
                                                     size_t rowsPerBlock) {
    index_t myStartRow = rowStart + (blockIdx.x * rowsPerBlock);
    index_t myEndRow = myStartRow + rowsPerBlock;
    if (myEndRow > rowEnd) {
        myEndRow = rowEnd;
    }
    int myGroupID = threadIdx.x / THREADS_PER_ROW;
    int myThreadID = threadIdx.x % THREADS_PER_ROW;
    int stride = blockDim.x / THREADS_PER_ROW;
    for (index_t row = myStartRow + myGroupID; row < myEndRow; row += stride) {
        index_t curStartRow = rowOffsets[row];
        index_t curEndRow = rowOffsets[row + 1];
        data_t myVal = values[row];
        for (index_t idx = curStartRow + myThreadID; idx < curEndRow; idx += THREADS_PER_ROW) {
            index_t col = colIndices[idx];
            atomicAdd(&values[col], myVal);
            // printf("Thread %d adding %f to row %d. New value is %f. My id is %d\n", threadIdx.x, myVal, col,
            // values[col], row);
        }
    }
}

template <class index_t, class data_t>
__global__ void cudaTraversalUpMultiElementKernel(index_t* rowOffsets,
                                                  index_t* colIndices,
                                                  data_t* values,
                                                  size_t rowStart,
                                                  size_t rowEnd,
                                                  size_t rowsPerBlock,
                                                  index_t unit = 1) {
    index_t myStartRow = rowStart + (blockIdx.x * rowsPerBlock);
    index_t myEndRow = myStartRow + rowsPerBlock;
    if (myEndRow > rowEnd) {
        myEndRow = rowEnd;
    }
    int myElementID = threadIdx.x % unit;
    int myOffset = threadIdx.x / unit;
    int stride = blockDim.x / unit;

    extern __shared__ char buf[];
    index_t* colBuffer = (index_t*)(buf);

    for (index_t row = myStartRow; row < myEndRow; row++) {
        index_t curStartRow = rowOffsets[row];
        index_t curEndRow = rowOffsets[row + 1];
        index_t rowNNZ = curEndRow - curStartRow;
        index_t totalIters = (rowNNZ + stride - 1) / stride;
        data_t sum = 0;

        for (index_t i = 0; i < totalIters; i++) {
            // preload colBuffer if necessary
            // need to avoid deadlock due to sync threads
            index_t curStep = i % COL_BUFFER_ITER;

            if (curStep == 0) {
                __syncthreads();
                for (int j = threadIdx.x; j < stride * COL_BUFFER_ITER; j += blockDim.x) {
                    index_t idx = curStartRow + i * stride + j;
                    if (idx >= curEndRow) {
                        break;
                    }
                    colBuffer[j] = colIndices[idx];
                }
                __syncthreads();
            }

            index_t idx = curStartRow + i * stride + myOffset;
            if (idx >= curEndRow) {
                continue;
            }

            // index_t col = colIndices[idx];
            index_t col =
                colBuffer[myOffset + curStep * stride]; // colIndices[curStep * stride + myOffset];
            sum += values[col * unit + myElementID];
            // __syncthreads();
        }
        // Use atomicAdd to update the value
        atomicAdd(&values[row * unit + myElementID], sum);
    }
}

template <class index_t, class data_t>
__global__ void cudaTraversalDownMultiElementKernel(index_t* rowOffsets,
                                                    index_t* colIndices,
                                                    data_t* values,
                                                    size_t rowStart,
                                                    size_t rowEnd,
                                                    size_t rowsPerBlock,
                                                    index_t unit = 1) {
    index_t myStartRow = rowStart + (blockIdx.x * rowsPerBlock);
    index_t myEndRow = myStartRow + rowsPerBlock;
    if (myEndRow > rowEnd) {
        myEndRow = rowEnd;
    }
    int myElementID = threadIdx.x % unit;
    int myOffset = threadIdx.x / unit;
    int stride = blockDim.x / unit;

    extern __shared__ char buf[];
    index_t* colBuffer = (index_t*)(buf);

    for (index_t row = myStartRow; row < myEndRow; row++) {
        index_t curStartRow = rowOffsets[row];
        index_t curEndRow = rowOffsets[row + 1];
        index_t rowNNZ = curEndRow - curStartRow;
        index_t totalIters = (rowNNZ + stride - 1) / stride;
        data_t sum = 0;

        auto myVal = values[row * unit + myElementID];

        for (index_t i = 0; i < totalIters; i++) {
            // preload colBuffer if necessary
            // need to avoid deadlock due to sync threads
            index_t curStep = i % COL_BUFFER_ITER;

            if (curStep == 0) {
                __syncthreads();
                for (int j = threadIdx.x; j < stride * COL_BUFFER_ITER; j += blockDim.x) {
                    index_t idx = curStartRow + i * stride + j;
                    if (idx >= curEndRow) {
                        break;
                    }
                    colBuffer[j] = colIndices[idx];
                }
                __syncthreads();
            }

            index_t idx = curStartRow + i * stride + myOffset;
            if (idx >= curEndRow) {
                continue;
            }

            // index_t col = colIndices[idx];
            index_t col =
                colBuffer[myOffset + curStep * stride]; // colIndices[curStep * stride + myOffset];
            atomicAdd(&values[col * unit + myElementID], myVal);
            // __syncthreads();
        }
    }
}

} // namespace grgl

#endif // GRGL_CUDA_KERNELS_CUH