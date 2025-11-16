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
#include "grgl/cuda_kernels.cuh"


namespace grgl {


#define MIN_BLOCK_SIZE 64   // set to 64 for better occupancy
#define CONCURRENT_BLOCKS_PER_SM 32 // 2048 / MIN_BLOCK_SIZE
#define NUM_SMS 108  // for A100 GPU
#define MIN_LUNCHED_BLOCKS (CONCURRENT_BLOCKS_PER_SM * NUM_SMS)


template<class index_t, class data_t>
void cudaTraverseUpSingleElementLauncher(
    index_t* row_offsets,
    index_t* col_indices, 
    data_t* values, 
    size_t start_row, 
    size_t end_row, 
    size_t heavy_cutoff,
    double avg_edge_per_node, 
    cudaStream_t stream) {

    if (heavy_cutoff != 0) {
        // std::cout << "Warning: heavy_cutoff is not zero in cudaTraverseUpSingleElementLauncher. This may lead to performance issues." << std::endl;
        int block_size = MIN_BLOCK_SIZE;
        int rows_per_block = 2;
        int num_blocks = (heavy_cutoff + rows_per_block - 1) / rows_per_block;
        cudaTraversalUpSingleElementKernel<index_t, data_t, 32><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            start_row + heavy_cutoff,
            rows_per_block
        );
        start_row += heavy_cutoff;
    }

    int block_size = MIN_BLOCK_SIZE;
    int num_rows = end_row - start_row;
    int rows_per_block_max_occupancy = (num_rows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rows_per_block = 32;
    if (rows_per_block_max_occupancy < rows_per_block) {
        rows_per_block = rows_per_block_max_occupancy;
    }
    int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;

    /*
    std::cout << "Launching cudaTraverseUpSingleElementLauncher with block_size " << block_size 
              << ", rows_per_block " << rows_per_block 
              << ", num_blocks " << num_blocks 
              << ", avg_edge_per_node " << avg_edge_per_node 
              << std::endl;
    */

    if (avg_edge_per_node <= 8){
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalUpSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 16){
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalUpSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 32){
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalUpSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalUpSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    }
}

template<class index_t, class data_t>
void cudaTraverseDownSingleElementLauncher(
    index_t* row_offsets,
    index_t* col_indices, 
    data_t* values, 
    size_t start_row, 
    size_t end_row, 
    size_t heavy_cutoff, 
    double avg_edge_per_node, 
    cudaStream_t stream) {

    if (heavy_cutoff != 0) {
        // std::cout << "Warning: heavy_cutoff is not zero in cudaTraverseDownSingleElementLauncher. This may lead to performance issues." << std::endl;
        int block_size = MIN_BLOCK_SIZE;
        int rows_per_block = 1;
        int num_blocks = (heavy_cutoff + rows_per_block - 1) / rows_per_block;
        cudaTraversalDownSingleElementKernel<index_t, data_t, 32><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            start_row + heavy_cutoff,
            rows_per_block
        );
        start_row += heavy_cutoff;
    }

    int block_size = MIN_BLOCK_SIZE;
    int num_rows = end_row - start_row;
    int rows_per_block_max_occupancy = (num_rows + MIN_LUNCHED_BLOCKS - 1) / MIN_LUNCHED_BLOCKS;
    int rows_per_block = 16;
    if (rows_per_block_max_occupancy < rows_per_block) {
        rows_per_block = rows_per_block_max_occupancy;
    }

    int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;

    // std::cout << "Launching cudaTraverseDownSingleElementLauncher with " << num_blocks << " blocks, block size " << block_size << ", rows per block " << rows_per_block << ", avg_edge_per_node " << avg_edge_per_node << std::endl;

    if (avg_edge_per_node <= 8){
        constexpr int THREADS_PER_ROW = 4;
        cudaTraversalDownSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 16){
        constexpr int THREADS_PER_ROW = 8;
        cudaTraversalDownSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else if (avg_edge_per_node <= 32){
        constexpr int THREADS_PER_ROW = 16;
        cudaTraversalDownSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    } else {
        constexpr int THREADS_PER_ROW = 32;
        cudaTraversalDownSingleElementKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size, 0, stream>>>(
            row_offsets,
            col_indices,
            values,
            start_row,
            end_row,
            rows_per_block
        );
    }
}


template<class index_t, class data_t>
void cudaTraverseUP(GPUGRG<index_t>& gpu_grg, data_t*d_inout_values, cudaStream_t& stream, index_t unit=1) {
    std::cout << "Starting GPU traversal..." << std::endl;

    // launch the kernels level by level
    for (size_t level = 1; level < gpu_grg.max_height; level++) {
        size_t start_row = gpu_grg.host_height_cutoffs[level];
        size_t end_row = gpu_grg.host_height_cutoffs[level+1];
        // size_t num_rows = end_row - start_row;
        index_t heavy_cutoff = gpu_grg.host_heavy_cutoffs[level];

        if (unit == 1) {
           cudaTraverseUpSingleElementLauncher<index_t, data_t>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values, 
                start_row, 
                end_row, 
                heavy_cutoff,
                gpu_grg.host_avg_child_counts[level], 
                stream
            );
        
        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;
            if (heavy_cutoff != 0) {
                int threads_per_block = HEAVY_THREADS_PER_UNIT * unit;
                int rows_per_block = 1;
                int num_blocks = (heavy_cutoff + rows_per_block - 1) / rows_per_block;
                
                cudaTraversalUpMultiElementKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * HEAVY_THREADS_PER_UNIT, stream>>>(
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

            cudaTraversalUpMultiElementKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * THREADS_PER_UNIT, stream>>>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row + heavy_cutoff,
                end_row,
                rows_per_block,
                unit
            );
        }
    // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

}

template<class index_t, class data_t>
void cudaTraverseDOWN(GPUGRG<index_t>& gpu_grg, data_t* d_inout_values, cudaStream_t stream, index_t unit=1) {
    std::cout << "Starting GPU traversal..." << std::endl;

    // launch the kernels level by level
    for (int level = gpu_grg.max_height - 1; level > 0; level--) {
        size_t start_row = gpu_grg.host_height_cutoffs[level];
        size_t end_row = gpu_grg.host_height_cutoffs[level+1];
        double avg_edge_per_node = gpu_grg.host_avg_child_counts[level];
        // size_t num_rows = end_row - start_row;
        size_t heavy_cutoff = gpu_grg.host_heavy_cutoffs[level];

        if (unit == 1) {
            cudaTraverseDownSingleElementLauncher<index_t, data_t>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row,
                end_row,
                heavy_cutoff,
                avg_edge_per_node,
                stream
            );
        } else {
            constexpr int THREADS_PER_UNIT = 4;
            constexpr int HEAVY_THREADS_PER_UNIT = 48;

            if (heavy_cutoff != 0) {
                int threads_per_block = HEAVY_THREADS_PER_UNIT * unit;
                int rows_per_block = 1;
                int num_blocks = (heavy_cutoff + rows_per_block - 1) / rows_per_block;
                cudaTraversalDownMultiElementKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * HEAVY_THREADS_PER_UNIT, stream>>>(
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

            cudaTraversalDownMultiElementKernel<index_t, data_t><<<num_blocks, threads_per_block, COL_BUFFER_ITER * sizeof(index_t) * THREADS_PER_UNIT, stream>>>(
                gpu_grg.row_offsets,
                gpu_grg.col_indices,
                d_inout_values,
                start_row + heavy_cutoff,
                end_row,
                rows_per_block,
                unit
            );
            
        }
        // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

    std::cout << "Finished GPU traversal DOWN." << std::endl;
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
            cudaReorderMapKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
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
            cudaReorderMapKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
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

        cudaReorderMapKernel<index_t, data_t, false, true><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
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

        cudaTraverseUP<index_t, data_t>(
            *this, 
            buffer, 
            *(this->work_stream_ptr),
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
            cudaReorderMapKernel<index_t, data_t, true, false><<<num_blocks, threads_per_block, 0, *(this->work_stream_ptr)>>>(
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
