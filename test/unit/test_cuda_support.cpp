#ifdef GRGL_CUDA_ENABLED

#include <gtest/gtest.h>
#include "grgl/grg.h"
#include <memory>

using namespace grgl;

class CudaSupportTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a minimal GRG for testing
        grg = std::make_shared<MutableGRG>(4, 2, true);  // 4 samples, diploid, phased
    }
    
    std::shared_ptr<GRG> grg;
};

TEST_F(CudaSupportTest, CudaHardwareMustBeAvailable) {
    // This test FAILS if CUDA is compiled but no hardware is available
    bool cudaSupported = grg->hasCudaSupport();
    
    EXPECT_TRUE(cudaSupported) 
        << "CUDA support was compiled in but no CUDA hardware detected. "
        << "Either install CUDA drivers/hardware or build without -DENABLE_CUDA=ON";
    
    if (cudaSupported) {
        std::cout << "✓ CUDA hardware confirmed available" << std::endl;
    }
}

#endif // GRGL_CUDA_ENABLED