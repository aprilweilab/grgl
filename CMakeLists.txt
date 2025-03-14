cmake_minimum_required(VERSION 3.14)
project(grgl)
include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_AVX CACHE)
unset(COMPILER_SUPPORTS_AVX2 CACHE)
include(ExternalProject)

option(PYTHON_SUPPORT "Build Python loadable module" OFF)
option(ENABLE_CHECKERS "Enable external tools like clang-tidy" OFF)
option(ENABLE_TESTS "Enable automated test execution" ON)
option(ENABLE_BGEN "Enable BGEN support" OFF)
option(ENABLE_GSL "Enable GSL support for stat calculations" OFF)
option(ENABLE_CPU_PROF "Enable gperftools CPU profiler" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Check for optional POSIX functions
include(CheckSymbolExists)
check_symbol_exists(posix_fadvise "fcntl.h" HAVE_POSIX_FADVISE)

# Global flags that affect more than one target.
add_compile_options(-DCOMPACT_NODE_IDS -DVCF_GZ_SUPPORT -Werror)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory and run CMake from there.\n")
endif()

# Apple's clang does not properly detect the (lack of) these options.
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  message("Using AppleClang, so AVX manually disabled")
  set(COMPILER_SUPPORTS_AVX OFF)
  set(COMPILER_SUPPORTS_AVX2 OFF)
else()
  check_cxx_compiler_flag(-mavx COMPILER_SUPPORTS_AVX)
  check_cxx_compiler_flag(-mavx2 COMPILER_SUPPORTS_AVX2)
endif()

# The current potential use for AVX (beyond libvbyte, which is enabled elsewhere) is not
# sufficient performance improvement to enable. Enabling means we should support dynamic
# detection as well, since the Python bdist's are built with AVX available.
#if(${COMPILER_SUPPORTS_AVX} OR ${COMPILER_SUPPORTS_AVX2})
#  add_compile_options(-DUSE_AVX=1 -mavx)
#endif()

if (${ENABLE_CHECKERS})
  find_program(CLANGTIDY clang-tidy)
  if(CLANGTIDY)
    set(CMAKE_CXX_CLANG_TIDY ${CLANGTIDY} -extra-arg=-Wno-unknown-warning-option)
    message("Clang-Tidy finished setting up.")
  else()
    message(SEND_ERROR "Clang-Tidy requested but executable not found.")
  endif()

  function(add_clang_format_target)
    if(NOT ${PROJECT_NAME}_CLANG_FORMAT_BINARY)
            find_program(${PROJECT_NAME}_CLANG_FORMAT_BINARY clang-format)
    endif()

    if(${PROJECT_NAME}_CLANG_FORMAT_BINARY)
      if(${PROJECT_NAME}_BUILD_EXECUTABLE)
        add_custom_target(clang-format
                          COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY}
                          -i ${exe_sources} ${headers}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
      elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
        add_custom_target(clang-format
                          COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY}
                          -i ${headers}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
      else()
        add_custom_target(clang-format
                          COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY}
                          -i ${sources} ${headers}
                          WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
      endif()
      message(STATUS "Format the project using the `clang-format` target (i.e: cmake --build build --target clang-format).\n")
    endif()
  endfunction()

  add_clang_format_target()
endif()

set(GRGL_TEST_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_bloom_filter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_common.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_construct.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_csr_storage.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_grg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_map_muts.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_mult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_mutation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_serialize.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_util.cpp
    ${CMAKE_CURRENT_LIST_DIR}/test/unit/test_visitor.cpp
  )

include(third-party/CMakeLists.txt)

set(GRGL_PUBLIC_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/include/grgl/grg.h
    ${CMAKE_CURRENT_LIST_DIR}/include/grgl/grgnode.h
    ${CMAKE_CURRENT_LIST_DIR}/include/grgl/mutation.h
    ${CMAKE_CURRENT_LIST_DIR}/include/grgl/ts2grg.h
  )

set(GRGL_CORE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/src/build_shape.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/hap_index.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/grg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/grg_merge.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/map_mutations.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/mut_iterator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/serialize.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/ts2grg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/windowing.cpp
  )

include_directories(${CMAKE_CURRENT_LIST_DIR}/include/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/tskit/c/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/args/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/picohash/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/picovcf/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/pybind11/
                    ${CMAKE_CURRENT_LIST_DIR}/third-party/bgen/include/
                    ${CMAKE_BINARY_DIR}/third-party/bgen/
                    ${CMAKE_CURRENT_LIST_DIR}/src/
                    )

if (${ENABLE_BGEN})
  set(BGEN_LIBS bgen)
  add_compile_options(-DBGEN_ENABLED=1)
endif()

# In order to use this, grab https://github.com/gperftools/gperftools/releases/download/gperftools-2.15/gperftools-2.15.tar.gz
# and unzip it to the top-level grgl repo directory, and do "./configure && make"
if (${ENABLE_CPU_PROF})
  include_directories(${CMAKE_CURRENT_LIST_DIR}/gperftools-2.15/src/)
  add_compile_options(-g) # Turn on symbols (but not debug mode)
  set(GRGP_LIBS ${CMAKE_CURRENT_LIST_DIR}/gperftools-2.15/.libs/libprofiler.a unwind)
endif()

# We don't need AMPL bindings from gsl.
if (${ENABLE_GSL})
  ExternalProject_Add(
    gsl
    SOURCE_DIR ${CMAKE_BINARY_DIR}/third-party/gsl
    BINARY_DIR ${CMAKE_BINARY_DIR}/third-party/gsl
    URL "https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz"
    CMAKE_ARGS ""
    CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/third-party/gsl/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND ""
    TEST_COMMAND "")

  set(GRGP_LIBS ${CMAKE_BINARY_DIR}/third-party/gsl/.libs/libgsl.a)
  set(GRGP_DEPS gsl)
  include_directories(${CMAKE_BINARY_DIR}/third-party/gsl/)
  add_compile_options(-DGSL_ENABLED=1)
endif()

# Library libgrgl.a for C++ projects that want to use GRGL.
add_library(grgl STATIC ${GRGL_CORE_SOURCES})
target_compile_options(grgl PRIVATE
                       -fPIC)

# Tool for doing some basic format conversion tasks.
add_executable(grgl_tool
               ${CMAKE_CURRENT_LIST_DIR}/src/grgl.cpp
               )
set_target_properties(grgl_tool PROPERTIES OUTPUT_NAME "grgl")
target_link_libraries(grgl_tool tskit grgl ${BGEN_LIBS} z vbyte)

# Tool for processing GRG files (stats, etc).
add_executable(grgp_tool
               ${CMAKE_CURRENT_LIST_DIR}/src/grgp.cpp
               ${CMAKE_CURRENT_LIST_DIR}/src/calculations.cpp
               )
set_target_properties(grgp_tool PROPERTIES OUTPUT_NAME "grgp")
target_link_libraries(grgp_tool grgl ${GRGP_LIBS} vbyte)
if (${ENABLE_GSL})
  add_dependencies(grgp_tool ${GRGP_DEPS})
endif()


# Tool for merging GRGs.
add_executable(grgl_merge
               ${CMAKE_CURRENT_LIST_DIR}/src/merge.cpp
               )
set_target_properties(grgl_merge PROPERTIES OUTPUT_NAME "grg-merge")
target_link_libraries(grgl_merge grgl vbyte)

# Tool for converting file formats.
add_executable(gconverter
               ${CMAKE_CURRENT_LIST_DIR}/src/gconverter.cpp
               ${CMAKE_CURRENT_LIST_DIR}/src/mut_iterator.cpp
               )
target_link_libraries(gconverter ${BGEN_LIBS} z)

# Tool for converting file formats.
add_executable(gtime
               ${CMAKE_CURRENT_LIST_DIR}/extra/gtime.cpp
               ${CMAKE_CURRENT_LIST_DIR}/src/mut_iterator.cpp
               )
target_link_libraries(gtime ${BGEN_LIBS} z)

# Tool for indexing input files.
add_executable(gindexer
               ${CMAKE_CURRENT_LIST_DIR}/src/gindexer.cpp
               )
target_link_libraries(gindexer ${BGEN_LIBS} z)

# Tool for timing VCF scanning.
add_executable(vcftime
               ${CMAKE_CURRENT_LIST_DIR}/extra/vcftime.cpp
               )
target_link_libraries(vcftime ${BGEN_LIBS} z)

install(TARGETS grgl_tool grgl_merge gconverter gindexer vcftime)

# Build Python support as "_grgl"
if(PYTHON_SUPPORT)
    add_subdirectory(third-party/pybind11)
    pybind11_add_module(_grgl src/python/_grgl.cpp)
    target_link_libraries(_grgl PRIVATE grgl tskit vbyte)
endif()


install(
  TARGETS ${PROJECT_NAME}
  COMPONENT library
  INCLUDES
  DESTINATION include
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(
  FILES ${GRGL_PUBLIC_HEADERS}
  COMPONENT headers
  DESTINATION include/grgl
)

if(ENABLE_TESTS)
  # https://github.com/google/googletest/tree/master/googletest#incorporating-into-an-existing-cmake-project

  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  # Now simply link against gtest or gtest_main as needed. Eg
  add_executable(grgl_test ${GRGL_TEST_SOURCES})
  target_link_libraries(grgl_test gtest_main tskit grgl z ${BGEN_LIBS} vbyte)
  add_test(NAME grgl_test COMMAND grgl_test vbyte)
endif()