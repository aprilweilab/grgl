#ifndef GRGL_CUDA_ENABLED
#error "This file should only be compiled when CUDA is enabled"
#endif

#include "grgl/cuda/gpu_grg.h"

namespace grgl {

void storeGPUGRGToDisk(const GPUGRG& gpu_grg, const std::string& filename) {
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
    uint64_t indexSize = sizeof(NodeIDSizeT);

    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    file.write(reinterpret_cast<const char*>(&numRows), sizeof(numRows));
    file.write(reinterpret_cast<const char*>(&numEdges), sizeof(numEdges));
    file.write(reinterpret_cast<const char*>(&maxHeight), sizeof(maxHeight));
    file.write(reinterpret_cast<const char*>(&numSamples), sizeof(numSamples));
    file.write(reinterpret_cast<const char*>(&numMutations), sizeof(numMutations));
    file.write(reinterpret_cast<const char*>(&indexSize), sizeof(indexSize));

    // Copy data from GPU to CPU for storage
    std::vector<NodeIDSizeT> h_RowOffsets(numRows + 1);
    // std::vector<NodeIDSizeT> host_col_indices(nnz);
    std::vector<NodeIDSizeT> h_oldToNewMapping(numRows);
    std::vector<NodeIDSizeT> h_newToOldMapping(numRows);
    std::vector<NodeIDSizeT> h_mutationNewPairs(numMutations * 2);

    cudaMemcpy(h_RowOffsets.data(), gpu_grg.rowOffsets, (numRows + 1) * sizeof(NodeIDSizeT), cudaMemcpyDeviceToHost);
    // cudaMemcpy(host_col_indices.data(), gpu_grg.col_indices,
    //            nnz * sizeof(NodeIDSizeT), cudaMemcpyDeviceToHost);
    cudaMemcpy(
        h_oldToNewMapping.data(), gpu_grg.oldToNewMapping, numRows * sizeof(NodeIDSizeT), cudaMemcpyDeviceToHost);
    cudaMemcpy(
        h_newToOldMapping.data(), gpu_grg.newToOldMapping, numRows * sizeof(NodeIDSizeT), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mutationNewPairs.data(),
               gpu_grg.mutationAndNewMapping,
               numMutations * 2 * sizeof(NodeIDSizeT),
               cudaMemcpyDeviceToHost);

    // Write data sections
    file.write(reinterpret_cast<const char*>(h_RowOffsets.data()), (numRows + 1) * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostColIndices.data()), numEdges * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostHeightCutoffs.data()), (maxHeight + 1) * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostHeavyCutoffs.data()), (maxHeight + 1) * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(gpu_grg.hostAvgChildCounts.data()), (maxHeight + 1) * sizeof(double));
    file.write(reinterpret_cast<const char*>(h_oldToNewMapping.data()), numRows * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(h_newToOldMapping.data()), numRows * sizeof(NodeIDSizeT));
    file.write(reinterpret_cast<const char*>(h_mutationNewPairs.data()), numMutations * 2 * sizeof(NodeIDSizeT));

    if (!file.good()) {
        throw std::runtime_error("Error writing to file: " + filename);
    }

    file.close();
}

GPUGRG loadGPUGRGFromDisk(const std::string& filename) {
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

    if (indexSize != sizeof(NodeIDSizeT)) {
        throw std::runtime_error("Incompatible index or data type size in file: " + filename);
    }

    // Create and allocate GPUGRG structure
    GPUGRG gpuGRG;
    gpuGRG.init(numRows, numSamples, numEdges, numMutations, maxHeight);

    // Read data sections
    std::vector<NodeIDSizeT> h_rowOffsets(numRows + 1);
    std::vector<NodeIDSizeT> h_colIndices(numEdges);
    std::vector<NodeIDSizeT> h_heightCutoffs(maxHeight + 1);
    std::vector<NodeIDSizeT> h_heavyCutoffs(maxHeight + 1);
    std::vector<double> h_avgChildCounts(maxHeight + 1);
    std::vector<NodeIDSizeT> h_oldToNewMapping(numRows);
    std::vector<NodeIDSizeT> h_newToOldMapping(numRows);
    std::vector<NodeIDSizeT> h_mutationNewPairs(numMutations * 2);

    file.read(reinterpret_cast<char*>(h_rowOffsets.data()), (numRows + 1) * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_colIndices.data()), numEdges * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_heightCutoffs.data()), (maxHeight + 1) * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_heavyCutoffs.data()), (maxHeight + 1) * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_avgChildCounts.data()), (maxHeight + 1) * sizeof(double));
    file.read(reinterpret_cast<char*>(h_oldToNewMapping.data()), numRows * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_newToOldMapping.data()), numRows * sizeof(NodeIDSizeT));
    file.read(reinterpret_cast<char*>(h_mutationNewPairs.data()), numMutations * 2 * sizeof(NodeIDSizeT));

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

GPUGRG convertGRGToGPUGRG(GRG* grg) {
    // Validate input requirements
    if (!grg->nodesAreOrdered()) {
        throw std::runtime_error("GRG nodes must be topologically ordered for GPU conversion");
    }

    if (grg->numNodes() == 0 || grg->numEdges() == 0) {
        throw std::runtime_error("GRG must have nodes and edges for GPU conversion");
    }

    // Check for potential overflow with NodeIDSizeT
    if (grg->numNodes() > std::numeric_limits<NodeIDSizeT>::max() ||
        grg->numEdges() > std::numeric_limits<NodeIDSizeT>::max()) {
        throw std::runtime_error("GRG too large for specified index type");
    }

    std::cout << "Converting GRG to GPUGRG format..." << std::endl;
    auto time_start = std::chrono::high_resolution_clock::now();

    // Create visitor to compute height-based ordering
    GPUCsrVisitor visitor(grg->numNodes(), grg->numEdges());
    auto seeds = grg->getRootNodes();

    // Perform DFS traversal to compute node heights
    grg->visitDfs(visitor, TraversalDirection::DIRECTION_DOWN, seeds);

    // Rearrange nodes within each height by child counts for better locality
    visitor.rearrange(grg);

    // Construct the GPUGRG structure
    GPUGRG gpu_grg;
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

} // namespace grgl