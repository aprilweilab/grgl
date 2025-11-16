#ifndef GRGL_CUDA_KERNELS_CUH
#define GRGL_CUDA_KERNELS_CUH

#ifdef GRGL_CUDA_ENABLED
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <chrono>

namespace grgl {

#define COL_BUFFER_ITER 16

template <class index_t, class data_t, bool INPUT_NODE_MAJOR, bool PERMUTATION_IS_SRC_TO_DST, bool ATOMIC_ADD=false>
__global__ void cudaReorderMapKernel(
    data_t* dst,
    const data_t* src,
    const index_t* permutation,
    size_t st,
    size_t ed,
    size_t num_unit,
    size_t num_node,
    index_t node_per_block
) {
    static_assert(INPUT_NODE_MAJOR == !PERMUTATION_IS_SRC_TO_DST, "For multi-element reorder, INPUT_NODE_MAJOR must be opposite to PERMUTATION_IS_SRC_TO_DST");
    index_t my_node = st + blockIdx.x * node_per_block + threadIdx.x / num_unit;
    if (my_node >= ed) {
        return;
    }
    index_t my_element_pos = threadIdx.x % num_unit;

    index_t src_id;
    index_t dst_id;

    if (INPUT_NODE_MAJOR) {
        src_id = permutation[my_node] * num_unit + my_element_pos;
        dst_id = my_node + my_element_pos * num_node;
    } else {
        src_id = my_element_pos * num_node + my_node;
        dst_id = permutation[my_node] * num_unit + my_element_pos;
    }

    if (!ATOMIC_ADD)
        dst[dst_id] = src[src_id];
    else
        atomicAdd(&dst[dst_id], src[src_id]);
}

template <class index_t, class data_t, bool INPUT_NODE_MAJOR, bool PERMUTATION_FWD_ORDER, bool ATOMIC_ADD=false>
__global__ void cudaReorderPermutationPairKernel(
    data_t* dst,
    const data_t* src,
    const index_t* permutation_pair,
    size_t st,
    size_t ed,
    size_t num_unit,
    size_t num_node,
    index_t node_per_block
) {
    static_assert(INPUT_NODE_MAJOR == !PERMUTATION_FWD_ORDER, "For reorder, INPUT_NODE_MAJOR must be opposite to PERMUTATION_FWD_ORDER");
    index_t my_node = st + blockIdx.x * node_per_block + threadIdx.x / num_unit;
    if (my_node >= ed) {
        return;
    }
    index_t my_element_pos = threadIdx.x % num_unit;

    index_t src_id;
    index_t dst_id;

    index_t ele1 = permutation_pair[my_node * 2];
    index_t ele2 = permutation_pair[my_node * 2 + 1];

    if (ele1 == INVALID_NODE_ID || ele2 == INVALID_NODE_ID) {
        // invalid mapping, skip
        return;
    }

    if (INPUT_NODE_MAJOR) {
        src_id = ele2 * num_unit + my_element_pos;
        dst_id = ele1 + my_element_pos * num_node;
    } else {
        src_id = ele1 + my_element_pos * num_node;
        dst_id = ele2 * num_unit + my_element_pos;
    }

    if (!ATOMIC_ADD) {
        dst[dst_id] = src[src_id];
    } else {
        atomicAdd(&dst[dst_id], src[src_id]);
    }
}


template<class index_t, class data_t, int THREADS_PER_ROW>
__global__ void cudaTraversalUpSingleElementKernel(
        index_t* row_offsets,
        index_t* col_indices,
        data_t* values,
        size_t start_row,
        size_t end_row,
        size_t rows_per_block) {
    index_t my_start_row = start_row + (blockIdx.x * rows_per_block);
    index_t my_end_row = my_start_row + rows_per_block;
    if (my_end_row > end_row) {
        my_end_row = end_row;
    }
    int my_group_id = threadIdx.x / THREADS_PER_ROW;
    int my_thread_id = threadIdx.x % THREADS_PER_ROW;
    int stride = blockDim.x / THREADS_PER_ROW;
    for (index_t row = my_start_row + my_group_id; row < my_end_row; row += stride ) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        data_t sum = 0;
        for (index_t idx = row_start + my_thread_id; idx < row_end; idx += THREADS_PER_ROW) {
            index_t col = col_indices[idx];
            sum += values[col];
        }
        // Reduce within the group. Use Warp shuffle.
        for (int offset = THREADS_PER_ROW / 2; offset > 0; offset /= 2) {
            sum += __shfl_down_sync(0xFFFFFFFF, sum, offset);
        }
        if (my_thread_id == 0) {
            values[row] = sum;
        }
    }
}

template<class index_t, class data_t, int THREADS_PER_ROW>
__global__ void cudaTraversalDownSingleElementKernel(
        index_t* row_offsets,
        index_t* col_indices,
        data_t* values,
        size_t start_row,
        size_t end_row,
        size_t rows_per_block) {
    index_t my_start_row = start_row + (blockIdx.x * rows_per_block);
    index_t my_end_row = my_start_row + rows_per_block;
    if (my_end_row > end_row) {
        my_end_row = end_row;
    }
    int my_group_id = threadIdx.x / THREADS_PER_ROW;
    int my_thread_id = threadIdx.x % THREADS_PER_ROW;
    int stride = blockDim.x / THREADS_PER_ROW;
    for (index_t row = my_start_row + my_group_id; row < my_end_row; row += stride ) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        data_t my_val = values[row];
        for (index_t idx = row_start + my_thread_id; idx < row_end; idx += THREADS_PER_ROW) {
            index_t col = col_indices[idx];
            atomicAdd(&values[col], my_val);
            // printf("Thread %d adding %f to row %d. New value is %f. My id is %d\n", threadIdx.x, my_val, col, values[col], row);
        }
    }
}

template<class index_t, class data_t>
__global__ void cudaTraversalUpMultiElementKernel(
        index_t* row_offsets,
        index_t* col_indices,
        data_t* values,
        size_t start_row,
        size_t end_row,
        size_t rows_per_block,
        index_t unit=1
    ) {
    index_t my_start_row = start_row + (blockIdx.x * rows_per_block);
    index_t my_end_row = my_start_row + rows_per_block;
    if (my_end_row > end_row) {
        my_end_row = end_row;
    }
    int my_element_id = threadIdx.x % unit;
    int my_offset = threadIdx.x / unit;
    int stride = blockDim.x / unit;

    extern __shared__ char buf[];
    index_t* col_buffer = (index_t*)(buf);

    for (index_t row = my_start_row; row < my_end_row; row++) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        index_t row_nnz = row_end - row_start;
        index_t total_iters = (row_nnz + stride - 1) / stride;
        data_t sum = 0;

        for (index_t i = 0; i < total_iters; i++) {
            // preload col_buffer if necessary
            // need to avoid deadlock due to sync threads
            index_t current_step = i % COL_BUFFER_ITER;

            if (current_step == 0) {
                __syncthreads();
                for (int j = threadIdx.x; j < stride * COL_BUFFER_ITER; j += blockDim.x) {
                    index_t idx = row_start + i * stride + j;
                    if (idx >= row_end) {
                        break;
                    }
                    col_buffer[j] = col_indices[idx];
                }
                __syncthreads();
            }

            index_t idx = row_start + i * stride + my_offset;
            if (idx >= row_end) {
                continue;
            }

            //index_t col = col_indices[idx];
            index_t col = col_buffer[my_offset + current_step * stride]; //col_indices[current_step * stride + my_offset];
            sum += values[col * unit + my_element_id];
            // __syncthreads();
        }
        // Use atomicAdd to update the value
        atomicAdd(&values[row * unit + my_element_id], sum);
        
    }
}

template<class index_t, class data_t>
__global__ void cudaTraversalDownMultiElementKernel(
        index_t* row_offsets,
        index_t* col_indices,
        data_t* values,
        size_t start_row,
        size_t end_row,
        size_t rows_per_block,
        index_t unit=1
    ) {
    index_t my_start_row = start_row + (blockIdx.x * rows_per_block);
    index_t my_end_row = my_start_row + rows_per_block;
    if (my_end_row > end_row) {
        my_end_row = end_row;
    }
    int my_element_id = threadIdx.x % unit;
    int my_offset = threadIdx.x / unit;
    int stride = blockDim.x / unit;

    extern __shared__ char buf[];
    index_t* col_buffer = (index_t*)(buf);

    for (index_t row = my_start_row; row < my_end_row; row++) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        index_t row_nnz = row_end - row_start;
        index_t total_iters = (row_nnz + stride - 1) / stride;
        data_t sum = 0;

        auto my_val = values[row * unit + my_element_id];

        for (index_t i = 0; i < total_iters; i++) {
            // preload col_buffer if necessary
            // need to avoid deadlock due to sync threads
            index_t current_step = i % COL_BUFFER_ITER;

            if (current_step == 0) {
                __syncthreads();
                for (int j = threadIdx.x; j < stride * COL_BUFFER_ITER; j += blockDim.x) {
                    index_t idx = row_start + i * stride + j;
                    if (idx >= row_end) {
                        break;
                    }
                    col_buffer[j] = col_indices[idx];
                }
                __syncthreads();
            }

            index_t idx = row_start + i * stride + my_offset;
            if (idx >= row_end) {
                continue;
            }

            //index_t col = col_indices[idx];
            index_t col = col_buffer[my_offset + current_step * stride]; //col_indices[current_step * stride + my_offset];
            atomicAdd(&values[col * unit + my_element_id], my_val);
            // __syncthreads();
        }
        
    }
}

} // namespace grgl

#endif // GRGL_CUDA_ENABLED
#endif // GRGL_CUDA_KERNELS_CUH