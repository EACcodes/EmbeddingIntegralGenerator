cmake_minimum_required(VERSION 2.6)
project(EmbeddingIntegralGenerator C CXX)


#**************************
# Compiler options
#**************************
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_AR "xiar")   # needed for libraries with -ipo
  if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo -xHost -opt-prefetch -unroll4")
  endif()
endif()

option(OMP "if OMP is on then the program is compiler with OpenMP support" OFF)
if (OMP)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp" )
  else()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
  endif()
endif()

# Make sure we user O0 for debug
if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

# Earlier CMAKE version don't have direct handling for C++11
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#*****************
# Math libraries
#*****************

# BLAS/LAPACK
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  link_directories(/usr/lib/atlas-base)
  link_directories(/usr/lib/atlas-base/atlas)
else()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")
endif()


#*******************
# Other libraries
#*******************

# Armadillo
include_directories(/home/EAC/TigerCI_libs/armadillo/include)

# FFTW
include_directories(/home/EAC/TigerCI_libs/fftw3/include)
link_directories(/home/EAC/TigerCI_libs/fftw3/lib/)

#********************
# Building EmbeddingIntegralGenerator
#********************

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/integral_generator)

# Define all the sub-parts and EmbeddingIntegralGenerator
add_subdirectory(integral_generator)
add_executable(emb_ints.exe ${INTEGRAL_GENERATOR_SOURCE})

# Link EmbeddingIntegralGenerator to FFTW3
target_link_libraries(emb_ints.exe fftw3)

# if we are using the GNU compiler we need to explicit ask to link to LAPACK/BLAS
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  target_link_libraries(emb_ints.exe lapack cblas f77blas atlas)
endif()

