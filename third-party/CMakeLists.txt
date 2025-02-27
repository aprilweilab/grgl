cmake_minimum_required(VERSION 3.10)
project(grgl)
################ TSKIT converted to CMake ################
# Source files and compiler options copied from tskit/c/meson.build
add_library(tskit STATIC
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/core.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/tables.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/trees.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/genotypes.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/stats.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/convert.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/tskit/haplotype_matching.c
            ${CMAKE_CURRENT_LIST_DIR}/tskit/c/subprojects/kastore/kastore.c
            )
target_compile_options(tskit PRIVATE
                       -Wall
                       -Wextra
                       -Werror
                       -Wpedantic
                       -W
                       -Wmissing-prototypes
                       -Wstrict-prototypes
                       -Wconversion
                       -Wshadow
                       -Wpointer-arith
                       -Wcast-align
                       -Wcast-qual
                       -Wwrite-strings
                       -Wnested-externs
                       -fshort-enums
                       -fno-common
                       -fPIC)

include_directories(${CMAKE_CURRENT_LIST_DIR}/tskit/c/
                    ${CMAKE_CURRENT_LIST_DIR}/tskit/c/subprojects/kastore/
                    )
##########################################################

if (${ENABLE_BGEN})
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/bgen/)
endif()

############## libvbyte converted to CMake ###############
if (${COMPILER_SUPPORTS_AVX2})
  set(VBYTE_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/libvbyte/varintdecode.c
      ${CMAKE_CURRENT_LIST_DIR}/libvbyte/vbyte.cc)
  set(AVX_OPTION -mavx -mavx2)
elseif (${COMPILER_SUPPORTS_AVX})
  set(VBYTE_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/libvbyte/varintdecode.c
      ${CMAKE_CURRENT_LIST_DIR}/libvbyte/vbyte.cc)
  target_compile_options(vbyte PRIVATE -mavx)
  set(AVX_OPTION -mavx)
else()
  set(VBYTE_SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/libvbyte/vbyte.cc)
  set(AVX_OPTION)
endif()
add_library(vbyte STATIC ${VBYTE_SOURCES})
target_compile_options(vbyte PRIVATE
                       -Wall
                       -Wextra
                       -pedantic
                       -fPIC
                       ${AVX_OPTION}
                       )
include_directories(${CMAKE_CURRENT_LIST_DIR}/libvbyte)

##########################################################
