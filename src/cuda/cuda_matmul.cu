// #define GRGL_CUDA_ENABLED
#ifdef GRGL_CUDA_ENABLED
#define WARMUP_GPU
#ifndef GRGL_CUDA_MATMUL_CU
#define GRGL_CUDA_MATMUL_CU

#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <chrono>
#include "grgl/grg.h"
#include "grgl/cuda_matmul.h"

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
        double_t my_val = values[row];
        for (index_t idx = row_start + my_thread_id; idx < row_end; idx += THREADS_PER_ROW) {
            index_t col = col_indices[idx];
            atomicAdd(&values[col], my_val);
            // printf("Thread %d adding %f to row %d. New value is %f. My id is %d\n", threadIdx.x, my_val, col, values[col], row);
        }
    }
}

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

    for (index_t row = my_start_row; row < my_end_row; row++) {
        index_t row_start = row_offsets[row];
        index_t row_end = row_offsets[row + 1];
        data_t sum = 0;
        for (index_t idx = row_start + my_offset; idx < row_end; idx += stride) {
            index_t col = col_indices[idx];
            sum += values[col * unit + my_element_id];
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

// Reorder data on CPU
// The i-th element in src is moved to permutation[i] in dst
// NOT optimized for performance
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

template<class index_t, class data_t>
void cudaTraverseUP(CudaCSR<index_t, data_t>& gpu_csr, data_t* h_sample_values, data_t* h_output_values, index_t unit=1) {
    // 1. set the initial values
    // assume the input values are reordered according to the GPU CSR format
    data_t* values;
    cudaMalloc((void**)&values, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemset(values, 0, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemcpy(values, h_sample_values, gpu_csr.host_height_cutoffs[1] * sizeof(data_t) * unit, cudaMemcpyHostToDevice);

    // print level and level cutoffs
    std::cout << "GPU CSR max height: " << gpu_csr.max_height <<", height cutoffs: ";
    for (size_t level = 0; level <= gpu_csr.max_height; level++) {
        std::cout << gpu_csr.host_height_cutoffs[level] << " ";
    }
    std::cout << std::endl;

    // use a C++ timer
    auto time_start = std::chrono::high_resolution_clock::now();
    std::cout << "Starting GPU traversal..." << std::endl;
    // 2. launch the kernels level by level

    // multiple iters here may not help
    for (size_t level = 1; level < gpu_csr.max_height; level++) {
        size_t start_row = gpu_csr.host_height_cutoffs[level];
        size_t end_row = gpu_csr.host_height_cutoffs[level+1];
        size_t num_rows = end_row - start_row;

        if (unit == 1) {
            constexpr int THREADS_PER_ROW = 8;
            int block_size = 64;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalUPSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size>>>(
                gpu_csr.row_offsets,
                gpu_csr.col_indices,
                values,
                start_row,
                end_row,
                rows_per_block
            );
        } else {
            constexpr int THREADS_PER_UNIT = 10;
            int threads_per_block = unit * THREADS_PER_UNIT;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalUPMultiElementSingleLevelKernel<index_t, data_t><<<num_blocks, threads_per_block>>>(
                gpu_csr.row_offsets,
                gpu_csr.col_indices,
                values,
                start_row,
                end_row,
                rows_per_block,
                unit
            );
        }
        // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

    cudaDeviceSynchronize();
    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "GPU traversal finished in " << duration << " ms for " << 1 << " iteration" << std::endl;

    // gpu_csr.print();
    //std::vector<data_t> values_host(gpu_csr.num_rows);
    std::cout << "Finished GPU traversal." << std::endl;
    
    // copy back the result to host
    cudaMemcpy(h_output_values, values, gpu_csr.num_rows * sizeof(data_t) * unit, cudaMemcpyDeviceToHost);
    cudaFree(values);
}

template<class index_t, class data_t>
void cudaTraverseDOWN(CudaCSR<index_t, data_t>& gpu_csr, data_t* h_mutation_values, data_t* h_output_values, index_t unit=1) {
    // 1. set the initial values
    // assume the input values are reordered according to the GPU CSR format
    data_t* values;
    cudaMalloc((void**)&values, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemset(values, 0, gpu_csr.num_rows * sizeof(data_t) * unit);
    cudaMemcpy(values, h_mutation_values, (gpu_csr.num_rows) * sizeof(data_t) * unit, cudaMemcpyHostToDevice);

    // print level and level cutoffs
    std::cout << "GPU CSR max height: " << gpu_csr.max_height <<", height cutoffs: ";
    for (size_t level = 0; level <= gpu_csr.max_height; level++) {
        std::cout << gpu_csr.host_height_cutoffs[level] << " ";
    }
    std::cout << std::endl;

    // use a C++ timer
    auto time_start = std::chrono::high_resolution_clock::now();
    std::cout << "Starting GPU traversal..." << std::endl;
    // 2. launch the kernels level by level
    for (int level = gpu_csr.max_height - 1; level > 0; level--) {
        size_t start_row = gpu_csr.host_height_cutoffs[level];
        size_t end_row = gpu_csr.host_height_cutoffs[level+1];
        size_t num_rows = end_row - start_row;

        if (unit == 1) {
            constexpr int THREADS_PER_ROW = 8;
            int block_size = 64;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalDOWNSingleLevelKernel<index_t, data_t, THREADS_PER_ROW><<<num_blocks, block_size>>>(
                gpu_csr.row_offsets,
                gpu_csr.col_indices,
                values,
                start_row,
                end_row,
                rows_per_block
            );
        } else {
            constexpr int THREADS_PER_UNIT = 10;
            int threads_per_block = unit * THREADS_PER_UNIT;
            int rows_per_block = 32;
            int num_blocks = (num_rows + rows_per_block - 1) / rows_per_block;
            cudaGRGTraversalDOWNMultiElementSingleLevelKernel<index_t, data_t><<<num_blocks, threads_per_block>>>(
                gpu_csr.row_offsets,
                gpu_csr.col_indices,
                values,
                start_row,
                end_row,
                rows_per_block,
                unit
            );
            
        }
        // std::cout << "Launched kernel for level " << level << " with start_row " << start_row << " and end_row " << end_row << std::endl;
    }

    cudaDeviceSynchronize();
    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "GPU traversal finished in " << duration << " ms." << std::endl;

    // gpu_csr.print();
    //std::vector<data_t> values_host(gpu_csr.num_rows);
    std::cout << "Finished GPU traversal." << std::endl;
    
    // copy back the result to host
    cudaMemcpy(h_output_values, values, gpu_csr.num_rows * sizeof(data_t) * unit, cudaMemcpyDeviceToHost);
    cudaFree(values);
}

template <typename IOType, typename NodeValueType, bool useBitVector>
void GRG::matrixMultiplicationGPU(const IOType* inputMatrix,
                               size_t inputCols,
                               size_t inputRows,
                               TraversalDirection direction,
                               IOType* outputMatrix,
                               size_t outputSize,
                               bool emitAllNodes,
                               bool byIndividual,
                               const IOType* initMatrix,
                               NodeInitEnum nodeInit,
                               IOType* missMatrix) {
    
    // make sure IOType and NodeValueType are the same type
    static_assert(std::is_same<IOType, NodeValueType>::value, "IOType and NodeValueType must be the same type for GPU matmul.");
    // make data_t an alias for IOType/NodeValueType

    using data_t = IOType;
    using index_t = uint32_t;

    release_assert(inputCols > 0);
    release_assert(inputRows > 0);
    release_assert(useBitVector == false); // Bitvector not supported on GPU yet
    // cols are num of samples or mutations
    // rows are num of different input vectors per sample/mutation
    const size_t outputCols = outputSize / inputRows;
    validateMatMulInputs(this, inputCols, inputRows, direction, outputSize, emitAllNodes, byIndividual, outputCols);
    // When we do bitvector calculations, we must have the number of input rows be a multiple of the
    // element size the bitvector is using.
    const size_t effectiveInputRows =
        useBitVector ? roundUpToMultiple(inputRows, sizeof(NodeValueType) * 8) : inputRows;

    // The node value storage stores the different row values consecutively (so
    // similar to column-major order), in contrast to the input/output matrices
    // which are row-major.
    const size_t nodesPerElem = 1; //useBitVector ? sizeof(NodeValueType) * 8 : 1;
    release_assert(effectiveInputRows % nodesPerElem == 0);
    std::vector<NodeValueType> nodeValues(numNodes() * effectiveInputRows / nodesPerElem);

    // release_assert(effectiveInputRows / nodesPerElem == 1); // only support 1 row for GPU now

    switch (nodeInit) {
    case NIE_XTX:
        release_assert(false); // NIE_XTX not supported on GPU yet
        break;
    case NIE_VECTOR:
        release_assert(false); // NIE_VECTOR not supported on GPU yet
        break;
    case NIE_MATRIX:
        release_assert(false); // NIE_MATRIX not supported on GPU yet
        break;
    case NIE_ZERO:
    default: break;
    }

    if (!this->nodesAreOrdered()) {
        release_assert(false); // Non-ordered nodes not supported on GPU
    }

    std::vector<data_t> output_reordered(numNodes() * effectiveInputRows / nodesPerElem, 0);

    // convert the graph to GPU CSR format
    // currently this only supports ordered nodes, and is executed every time for matrixMultiplicationGPU
    // This can be saved onto disk later if needed
    GPUCsrVisitor<index_t, data_t> visitor(numNodes(), numEdges());
    auto seeds = this->getRootNodes();

    auto time_st = std::chrono::high_resolution_clock::now();
    this->visitDfs(
        visitor, 
        TraversalDirection::DIRECTION_DOWN,
        seeds
    );
    CudaCSR<index_t, data_t> gpu_csr;
    visitor.constructCSR(this, gpu_csr);
    auto time_ed = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_ed - time_st).count();
    std::cout << "Converted GRG to GPU CSR format with " << gpu_csr.num_rows << " rows and " << gpu_csr.nnz << " non-zeros." << std::endl;
    std::cout << "Time taken for conversion: " << duration << " ms" << std::endl;
    visitor.printChildCounts();

    if (direction == DIRECTION_DOWN) {
        // Downward, we are calculating "how do the mutations impact the samples?"
        for (const auto& tuple : this->getNodesAndMutations<GRG::NodeMutMiss>()) {
            const NodeID& nodeId = std::get<0>(tuple);
            const MutationId& mutId = std::get<1>(tuple);
            const NodeID& missingnessNode = std::get<2>(tuple);
            assert(mutId < inputCols);
            if (nodeId != INVALID_NODE_ID) {
                const size_t base = nodeId * effectiveInputRows;
                for (size_t row = 0; row < inputRows; row++) {
                    const size_t rowStart = row * inputCols;
                    matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                        nodeValues.data(), base + row, inputMatrix, rowStart + mutId);
                    // std::cout << "Initial input for node " << nodeId << " at row " << row << ": " << nodeValues[base + row] << " Mutation ID is " << mutId << std::endl;
                }
            }
            if (missMatrix != nullptr && missingnessNode != INVALID_NODE_ID) {
                const size_t base = missingnessNode * effectiveInputRows;
                // Apply the single missing-data value for this mutation to all the node values
                // associated with the missingness node.
                matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                    nodeValues.data(), base, missMatrix[mutId], inputRows);
            }
        }

        std::vector<data_t> input_reordered(numNodes() * effectiveInputRows, 0);

        std::vector<data_t> output(numNodes() * effectiveInputRows, 0);

        reorder<index_t, data_t>(input_reordered.data(), nodeValues.data(), visitor.getNewId(), 0, numNodes(), effectiveInputRows);

#ifdef WARMUP_GPU
        // warm up
        for (size_t i = 0; i < 3; i++) {
            cudaTraverseDOWN<index_t, data_t>(gpu_csr, input_reordered.data(), output.data(), effectiveInputRows);
        }
#endif

        cudaTraverseDOWN<index_t, data_t>(gpu_csr, input_reordered.data(), output.data(), effectiveInputRows);

        reorder<index_t, data_t>(output_reordered.data(), output.data(), visitor.getOldId(), 0, numNodes(), effectiveInputRows);

        if (!emitAllNodes) {
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * outputCols;
                for (NodeID sampleId = 0; sampleId < numSamples(); sampleId++) {
                    const size_t base = sampleId * effectiveInputRows;
                    assert(rowStart + sampleId < outputSize);
                    const size_t sampleIndex = byIndividual ? (sampleId / m_ploidy) : sampleId;
                    matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                        outputMatrix, rowStart + sampleIndex, output_reordered.data(), base + row);
                }
            }
        }

    } else {
        // Upward, we are calculating "how do the samples impact the mutations?"
        // data_t* raw_input = new data_t[numSamples()];
        // std::vector<data_t> raw_input(numSamples() * effectiveInputRows, 0);
        // This is for initial values
        // release_assert(effectiveInputRows == 1); // only support 1 row for GPU now
        for (NodeID sampleId = 0; sampleId < numSamples(); sampleId++) {
            assert(sampleId < inputCols);
            const size_t base = sampleId * effectiveInputRows;
            const size_t sampleIndex = byIndividual ? (sampleId / m_ploidy) : sampleId;
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * inputCols;
                matmulPerformIOAddition<NodeValueType, IOType, useBitVector>(
                   nodeValues.data(), base + row, inputMatrix, rowStart + sampleIndex);
                // std::cout << "Initial input for sample " << sampleId << " at row " << row << ": " << raw_input[base + row] << std::endl;
            }
        }
        // reorder the input according to the GPU CSR format
        std::vector<data_t> values(numSamples() * effectiveInputRows, 0);
        std::vector<data_t> output(numNodes() * effectiveInputRows, 0);

        reorder<index_t, data_t>(values.data(), nodeValues.data(), visitor.getNewId(), 0, numSamples(), effectiveInputRows);

// #define GRGL_GPU_CONTROLLED_EXP
#ifdef GRGL_GPU_CONTROLLED_EXP
        auto max_height = gpu_csr.max_height;
        for (int height = 1; height < max_height; height++) {
            std::cout << "Processing height " << height << std::endl;
            gpu_csr.max_height = height;
            cudaTraverseUP<index_t, data_t>(gpu_csr, values.data(), output.data(), effectiveInputRows);
        }   
        gpu_csr.max_height = max_height;
#endif

#ifdef WARMUP_GPU
        // warm up
        for (size_t i = 0; i < 3; i++) {
            cudaTraverseUP<index_t, data_t>(gpu_csr, values.data(), output.data(), effectiveInputRows);
        }
#endif

        cudaTraverseUP<index_t, data_t>(gpu_csr, values.data(), output.data(), effectiveInputRows);

        reorder<index_t, data_t>(output_reordered.data(), output.data(), visitor.getOldId(), 0, numNodes(), effectiveInputRows);

        if (!emitAllNodes) {            
            for (size_t row = 0; row < inputRows; row++) {
                const size_t rowStart = row * outputCols;
                for (const auto& triple : this->getNodesAndMutations<GRG::MutNodeMiss>()) {
                    const NodeID& nodeId = std::get<0>(triple);
                    const MutationId& mutId = std::get<1>(triple);
                    assert(rowStart + mutId < outputSize);
                    if (nodeId != INVALID_NODE_ID) {
                        const size_t base = nodeId * effectiveInputRows;
                        matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                            outputMatrix, rowStart + mutId, output_reordered.data(), base + row);
                    }
                    const NodeID& missingnessNode = std::get<2>(triple);
                    if (missMatrix != nullptr && missingnessNode != INVALID_NODE_ID) {
                        const size_t base = missingnessNode * effectiveInputRows;
                        // Add the missing-data value for this mutation to the output vector position
                        // associated with the mutation.
                        matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                            missMatrix[mutId], output_reordered.data(), base, inputRows);
                    }
                }
            }
                
        }
    }
    if (emitAllNodes) {
        // release_assert(false); // emitAllNodes not supported on GPU yet
        
        for (size_t row = 0; row < inputRows; row++) {
            const size_t rowStart = row * outputCols;
            for (NodeID nodeId = 0; nodeId < numNodes(); nodeId++) {
                const size_t base = nodeId * effectiveInputRows;
                matmulPerformIOAddition<IOType, NodeValueType, useBitVector>(
                    outputMatrix, rowStart + nodeId, output_reordered.data(), base + row);
            }
        }
        
    }
}

bool GRG::hasCudaSupport() const {
    int deviceCount;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    return (error == cudaSuccess && deviceCount > 0);
}

/*
template void GRG::matrixMultiplicationGPU<float, float, false>(
    const float*, size_t, size_t, TraversalDirection, float*, size_t, 
    bool, bool, const float*, NodeInitEnum, float*);
*/
template void GRG::matrixMultiplicationGPU<double, double, false>(
    const double*, size_t, size_t, TraversalDirection, double*, size_t, 
    bool, bool, const double*, NodeInitEnum, double*);

template void GRG::matrixMultiplicationGPU<int, int, false>(
    const int*, size_t, size_t, TraversalDirection, int*, size_t, 
    bool, bool, const int*, NodeInitEnum, int*);

template void GRG::matrixMultiplicationGPU<uint32_t, uint32_t, false>(
    const uint32_t*, size_t, size_t, TraversalDirection, uint32_t*, size_t, 
    bool, bool, const uint32_t*, NodeInitEnum, uint32_t*);

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

#endif // GRGL_CUDA_MATMUL_CU
#endif // GRGL_CUDA_ENABLED