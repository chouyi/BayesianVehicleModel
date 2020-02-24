cmake_minimum_required(VERSION 3.15)
project(BayesianVehicleModel)
# Find out the arch.
EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
message( STATUS "Architecture: ${ARCHITECTURE}" )
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-includes)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++17 -Wall -Dregister= ")
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
include_directories(.)
include_directories(src)
include_directories(/usr/local/include ${PROJECT_SOURCE_DIR}/include)
link_directories(/usr/local/lib /usr/lib ${PROJECT_SOURCE_DIR}/lib/mac64)

include(FindMPFI)
if (NOT MPFI_FOUND)
    message(FATAL_ERROR "ERROR: MPFI not found")
endif(NOT MPFI_FOUND)

if (NOT MPFR_FOUND)
    message(FATAL_ERROR "ERROR: MPFR not found.")
endif(NOT MPFR_FOUND)

set(POLYFORM_LIBRARIES libPolyForm.a)

find_package(MPFI REQUIRED)

find_package(GSL REQUIRED)
include_directories(GSL_INCLUDE_DIRS)

add_executable(UAV_CTModel src/blas0.hpp src/blas1_d.hpp src/linpack_d.hpp  src/blas0.cpp src/blas1_d.cpp
            src/dynamic.h src/dynamic.cpp  src/inference.h src/inference.cpp src/linpack_d.hpp src/linpack_d.cpp
            src/UAV_CTmodel.cpp
        )

target_link_libraries(UAV_CTModel
        ${MPFR_LIBRARIES}
        ${MPFI_LIBRARIES}
        ${GMP_LIBRARIES}
        ${GSL_LIBRARY}
        ${GSL_CBLAS_LIBRARY}
        m
        )