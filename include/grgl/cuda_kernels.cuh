#ifndef GRGL_CUDA_KERNELS_CUH
#define GRGL_CUDA_KERNELS_CUH

#ifdef GRGL_CUDA_ENABLED
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <chrono>

namespace grgl {

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

} // namespace grgl

#endif // GRGL_CUDA_ENABLED
#endif // GRGL_CUDA_KERNELS_CUH