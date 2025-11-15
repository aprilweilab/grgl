// #define GRGL_CUDA_ENABLED
#ifndef GRGL_CUDA_MATMUL_CU
#define GRGL_CUDA_MATMUL_CU

#ifdef GRGL_CUDA_ENABLED


#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <chrono>
// #include "grgl/grg.h"
// #include "grgl/cuda_matmul.h"
#include "grgl/gpu_grg.h"


namespace grgl {

template<class index_t, class data_t, int THREADS_PER_ROW>
__global__ void cudaGRGTraversalUPSingleLevelKernel(
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
__global__ void cudaGRGTraversalDOWNSingleLevelKernel(
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

#define COL_BUFFER_ITER 16


template<class index_t, class data_t>
__global__ void cudaGRGTraversalUPMultiElementSingleLevelKernel(
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

// this kernel deals with the situation where each node has multiple elements, i.e. the value for a node is a vector
// One threadblock is responsible for one row at a time.
// Unit should be smaller or equal to the number of threads per threadblock for this kernel
// Number of threads per threadblock should be a multiplication of unit
template<class index_t, class data_t>
__global__ void cudaGRGTraversalDOWNMultiElementSingleLevelKernel(
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

    for (index_t row = my_start_row; row < my_end_row; row ++) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        data_t my_val = values[row * unit + my_element_id];

        for (index_t idx = row_start + my_offset; idx < row_end; idx += stride) {
            index_t col = col_indices[idx];
            atomicAdd(&values[col * unit + my_element_id], my_val);
            // printf("Thread %d adding %f to row %d. New value is %f. My id is %d\n", threadIdx.x, my_val, col, values[col], row);
        }
    }
}

/*
// this kernel reorders a single element
template <class index_t, class data_t, bool PERMUTATION_IS_SRC_TO_DST>
__global__ void cudaReorderSingleElementKernel(
    data_t* dst,
    const data_t* src,
    const index_t* permutation,
    size_t st,
    size_t ed,
    size_t element_per_thread=1
) {
    index_t my_st = (blockIdx.x * blockDim.x + threadIdx.x) * element_per_thread + st;
    if (my_st >= ed) {
        return;
    }
    index_t my_ed = my_st + element_per_thread * blockDim.x;
    if (my_ed > ed) {
        my_ed = ed;
    }

    for (index_t i = my_st; i < my_ed; i += element_per_thread) {
        if (PERMUTATION_IS_SRC_TO_DST) {
            index_t new_idx = permutation[i];
            dst[new_idx] = src[i];
        } else {
            index_t old_idx = permutation[i];
            dst[i] = src[old_idx];
        }
    }
}
*/


template <class index_t, class data_t, bool INPUT_NODE_MAJOR, bool PERMUTATION_IS_SRC_TO_DST, bool ATOMIC_ADD=false>
__global__ void cudaReorderKernel(
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



// Reorder data on CPU
// The i-th element in src is moved to permutation[i] in dst
// NOT optimized for performance
/*
template <class index_t, class data_t>
void reorder(data_t* dst, const data_t* src, const std::vector<index_t>& permutation, size_t st, size_t ed, size_t unit=1) {
    for (size_t i = st; i < ed; i++) {
        size_t new_idx = permutation[i];
        for (size_t u = 0; u < unit; u++) {
            dst[new_idx * unit + u] = src[i * unit + u];
        }
        // std::cout << "Reorder input: old idx " << i << " to new idx " << new_idx << " with value " << src[i] << std::endl;
    }
}
    */

template<class index_t, class data_t>
void cudaTraverseUPKernelLauncher(GPUGRG<index_t>& gpu_grg, data_t* values, size_t start_row, size_t end_row, double avg_edge_per_node, cudaStream_t stream) {

    int block_size = 64;
    int rows_per_block = 32;
    int num_rows = end_row - start_row;
    int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;

    if (avg_edge_per_node <= 8){
        constexpr int THREADS_PER_ROW = 4;
        cudaGRGTraversalUPSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            gpu_grg.row_offsets,
            gpu_grg.col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 16){
        constexpr int THREADS_PER_ROW = 8;
        cudaGRGTraversalUPSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            gpu_grg.row_offsets,
            gpu_grg.col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 32){
        constexpr int THREADS_PER_ROW = 16;
        cudaGRGTraversalUPSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            gpu_grg.row_offsets,
            gpu_grg.col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaGRGTraversalUPSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            gpu_grg.row_offsets,
            gpu_grg.col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    }
}


template<class index_t, class data_t>
void cudaTraverseUP(GPUGRG<index_t>& gpu_grg, data_t*d_inout_values, std::vector<cudaStream_t>& streams, index_t unit=1) {

    // print level and level cutoffs
    /*
    std::cout << "GPU CSR max height: " << gpu_csr.max_height <<", height cutoffs: ";
    for (size_t level = 0; level <= gpu_csr.max_height; level++) {
        std::cout << gpu_csr.host_height_cutoffs[level] << " ";
    }
    std::cout << std::endl;
    */

    // use a C++ timer
    // auto time_start = std::chrono::high_resolution_clock::now();
    std::cout << "Starting GPU traversal..." << std::endl;

    // 2. launch the kernels level by level
    for (size_t level = 1; level < gpu_grg.max_height; level++) {
        size_t start_row = gpu_grg.host_height_cutoffs[level];
        size_t end_row = gpu_grg.host_height_cutoffs[level+1];
        // size_t num_rows = end_row - start_row;
        index_t heavy_cutoff = gpu_grg.host_heavy_cutoffs[level];

        if (unit == 1) {
           cudaTraverseUPKernelLauncher<index_t, data_t>(gpu_grg, d_inout_values, start_row, end_row, gpu_grg.host_avg_child_counts[level], streams[0]);
        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;
            if (heavy_cutoff != 0) {
                int threads_per_block = HEAVY_THREADS_PER_UNIT * unit;
                int rows_per_block = 1;
                int num_blocks = (heavy_cutoff + rows_per_block - 1) / rows_per_block;
                
                cudaGRGTraversalUPMultiElementSingleLevelKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * HEAVY_THREADS_PER_UNIT, streams[1]>>>(
                    gpu_grg.row_offsets,
                    gpu_grg.col_indices,
                    d_inout_values,
                    start_row,
                    start_row + heavy_cutoff,
                    rows_per_block,
                    unit
                );
            }
            
            int threads_per_block = unit * THREADS_PER_UNIT;
            int rows_per_block = 1;
            int num_blocks = (end_row - start_row - heavy_cutoff + rows_per_block - 1) / rows_per_block;

            cudaGRGTraversalUPMultiElementSingleLevelKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * THREADS_PER_UNIT, streams[0]>>>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row + heavy_cutoff,
                end_row,
                rows_per_block,
                unit
            );
        }

        cudaDeviceSynchronize();
        std::cout<< "Warning: cudaTraverseUP is using cudaDeviceSynchronize for synchronization. This may cause performance issues if multiple streams are used." << std::endl;
        // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

    // auto time_end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    // std::cout << "GPU traversal finished in " << duration << " ms for " << 1 << " iteration" << std::endl;

    // gpu_csr.print();
    //std::vector<data_t> values_host(gpu_csr.num_rows);
    // std::cout << "Finished GPU traversal." << std::endl;
    
    // copy back the result to host
    // cudaMemcpy(h_output_values, values, gpu_csr.num_rows * sizeof(data_t) * unit, cudaMemcpyDeviceToHost);
    // cudaFree(values);
}

template<class index_t, class data_t>
void cudaTraverseDOWN(GPUGRG<index_t>& gpu_grg, data_t* d_inout_values, cudaStream_t stream, index_t unit=1) {
    // 1. set the initial values
    // assume the input values are reordered according to the GPU CSR format
    /*
    data_t* values;
    cudaMalloc((void**)&values, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemset(values, 0, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemcpy(values, h_mutation_values, (gpu_csr.num_rows) * sizeof(data_t) * unit, cudaMemcpyHostToDevice);
    */

    // print level and level cutoffs
    /*
    std::cout << "GPU CSR max height: " << gpu_csr.max_height <<", height cutoffs: ";
    for (size_t level = 0; level <= gpu_csr.max_height; level++) {
        std::cout << gpu_csr.host_height_cutoffs[level] << " ";
    }
    std::cout << std::endl;
    */

    // use a C++ timer
    // auto time_start = std::chrono::high_resolution_clock::now();
    std::cout << "Starting GPU traversal..." << std::endl;

    // 2. launch the kernels level by level
    for (int level = gpu_grg.max_height - 1; level > 0; level--) {
        size_t start_row = gpu_grg.host_height_cutoffs[level];
        size_t end_row = gpu_grg.host_height_cutoffs[level+1];
        size_t num_rows = end_row - start_row;

        if (unit == 1) {
            constexpr int THREADS_PER_ROW = 8;
            int block_size = 64;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalDOWNSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row,
                end_row,
                rows_per_block
            );
        } else {
            constexpr int THREADS_PER_UNIT = 10;
            int threads_per_block = unit * THREADS_PER_UNIT;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalDOWNMultiElementSingleLevelKernel<index_t, data_t><<<num_blocks, threads_per_block, 0, stream>>>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row,
                end_row,
                rows_per_block,
                unit
            );
            
        }
        // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

    /*
    cudaDeviceSynchronize();
    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "GPU traversal finished in " << duration << " ms." << std::endl;
    */

    // gpu_csr.print();
    //std::vector<data_t> values_host(gpu_csr.num_rows);
    std::cout << "Finished GPU traversal DOWN." << std::endl;
    
    // copy back the result to host
    /*
    cudaMemcpy(h_output_values, values, gpu_csr.num_rows * sizeof(data_t) * unit, cudaMemcpyDeviceToHost);
    cudaFree(values);
    */
}

template<class index_t>
template <class data_t>
void GPUGRG<index_t>::matrixMultiplication(const data_t* inputMatrix,
                               size_t inputCols,
                               size_t inputRows,
                               TraversalDirection direction,
                               data_t* outputMatrix,
                               size_t outputSize,
                               bool emitAllNodes,
                               data_t* buffer,
                               size_t buffer_size
                               ) {
    
    // make sure IOType and NodeValueType are the same type
    // static_assert(std::is_same<IOType, NodeValueType>::value, "IOType and NodeValueType must be the same type for GPU matmul.");
    // make data_t an alias for IOType/NodeValueType

    // using data_t = T;
    // using index_t = uint32_t;

    release_assert(inputCols > 0);
    release_assert(inputRows > 0);
    release_assert(buffer_size >= this->num_rows * inputRows * sizeof(data_t)); // buffer should be large enough to hold all node values
    cudaMemsetAsync(buffer, 0, this->num_rows * inputRows * sizeof(data_t), *(this->work_stream_ptr));
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
   auto feature_num = inputRows;   

    if (direction == DIRECTION_DOWN) {

        // Downward, we are calculating "how do the mutations impact the samples?"
        
        int node_per_block = 128 / feature_num;
        int threads_per_block = node_per_block * feature_num;
        int num_blocks = (this->num_mutations + node_per_block - 1) / node_per_block;

        cudaReorderPermutationPairKernel<index_t, data_t, false, true, true><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
            buffer,
            inputMatrix,
            this->mutation_and_new_mapping,
            0,
            this->num_mutations,
            feature_num,
            this->num_mutations,
            node_per_block
        );
        CHECK_CUDA_LAST_ERROR();

        cudaTraverseDOWN<index_t, data_t>(*this, buffer, *(this->work_stream_ptr), feature_num);
        CHECK_CUDA_LAST_ERROR();

        if (emitAllNodes) {
            num_blocks = (this->num_rows + node_per_block - 1) / node_per_block;
            cudaReorderKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
                outputMatrix,
                reinterpret_cast<data_t*>(buffer),
                this->old_to_new_mapping,
                0,
                this->num_rows,
                feature_num,
                this->num_rows,
                node_per_block
            );
        } else {
            cudaMemsetAsync(outputMatrix, 0, this->num_samples * feature_num * sizeof(data_t), *(this->work_stream_ptr));
            num_blocks = (this->num_samples + node_per_block - 1) / node_per_block;
            cudaReorderKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
                outputMatrix,
                reinterpret_cast<data_t*>(buffer),
                this->old_to_new_mapping,
                0,
                this->num_samples,
                feature_num,
                this->num_samples,
                node_per_block
            );
        }
        CHECK_CUDA_LAST_ERROR();

    } else {
        // Upward, we are calculating "how do the samples impact the mutations?"
        int node_per_block = 128 / feature_num;
        int threads_per_block = node_per_block * feature_num;
        int num_blocks = (this->num_samples + node_per_block - 1) / node_per_block;

        cudaReorderKernel<index_t, data_t, false, true><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
            reinterpret_cast<data_t*>(buffer),
            inputMatrix,
            this->old_to_new_mapping,
            0,
            this->num_samples,
            feature_num,
            this->num_samples,
            node_per_block
        );

        CHECK_CUDA_LAST_ERROR();

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

        cudaTraverseUP<index_t, data_t>(
            *this, 
            buffer, 
            streams, 
            feature_num
        );
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
            num_blocks = (this->num_rows + node_per_block - 1) / node_per_block;
            cudaReorderKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
                outputMatrix,
                reinterpret_cast<data_t*>(buffer),
                this->old_to_new_mapping,
                0,
                this->num_rows,
                feature_num,
                this->num_rows,
                node_per_block
            );
        } else {
            cudaMemsetAsync(outputMatrix, 0, this->num_mutations * feature_num * sizeof(data_t), *(this->work_stream_ptr));
            num_blocks = (this->num_mutations + node_per_block - 1) / node_per_block;
            cudaReorderPermutationPairKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
                outputMatrix,
                reinterpret_cast<data_t*>(buffer),
                this->mutation_and_new_mapping,
                0,
                this->num_mutations,
                feature_num,
                this->num_mutations,
                node_per_block
            );
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
template void GPUGRG<uint32_t>::matrixMultiplication<float>(
    const float*, size_t, size_t, TraversalDirection, 
    float*, size_t, bool, float*, size_t
);

template void GPUGRG<uint32_t>::matrixMultiplication<double>(
    const double*, size_t, size_t, TraversalDirection, 
    double*, size_t, bool, double*, size_t
);

template void GPUGRG<uint32_t>::matrixMultiplication<int>(
    const int*, size_t, size_t, TraversalDirection, 
    int*, size_t, bool, int*, size_t
);

/*
template void GPUGRG<uint32_t>::matrixMultiplication<long long>(
    const long long*, size_t, size_t, TraversalDirection, 
    long long*, size_t, bool, long long*, size_t
);
*/


template void GPUGRG<uint64_t>::matrixMultiplication<float>(
    const float*, size_t, size_t, TraversalDirection, 
    float*, size_t, bool, float*, size_t
);


template void GPUGRG<uint64_t>::matrixMultiplication<double>(
    const double*, size_t, size_t, TraversalDirection, 
    double*, size_t, bool, double*, size_t
);

template void GPUGRG<uint64_t>::matrixMultiplication<int>(
    const int*, size_t, size_t, TraversalDirection, 
    int*, size_t, bool, int*, size_t
);

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


/*
template void GRG::matrixMultiplicationGPU<float, float, true>(
    const float*, size_t, size_t, TraversalDirection, float*, size_t, 
    bool, bool, const float*, NodeInitEnum, float*);
*/
/*
template void GRG::matrixMultiplicationGPU<double, double, true>(
    const double*, size_t, size_t, TraversalDirection, double*, size_t, 
    bool, bool, const double*, NodeInitEnum, double*);
*/
} // namespace grgl

#endif // GRGL_CUDA_ENABLED
#endif // GRGL_CUDA_MATMUL_CU
