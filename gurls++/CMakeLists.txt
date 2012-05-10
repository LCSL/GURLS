#
# Copyright (C) 2012 Matteo Santoro
# CopyPolicy: Released under the terms of BSD
#

cmake_minimum_required(VERSION 2.6)

project(gurls++)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

enable_language(Fortran)
find_package(BLAS)
find_package(LAPACK)

set(GURLSLIBRARY gurls++)
set(SRCDIR src/)
set(HDRDIR include/)
set(TUTORIALDIR tutorial/)

set(gurls_headers
        ${HDRDIR}common.h
        ${HDRDIR}gurls.h
        ${HDRDIR}gurls.hpp
        ${HDRDIR}kernel.h
        ${HDRDIR}linearkernel.h
        ${HDRDIR}norm.h
        ${HDRDIR}exceptions.h
        ${HDRDIR}utils.h
        ${HDRDIR}basearray.h
        ${HDRDIR}basearray.hpp
        ${HDRDIR}gvec.h
        ${HDRDIR}gvec.hpp
        ${HDRDIR}gmat2d.h
        ${HDRDIR}gmat2d.hpp
        ${HDRDIR}gmath.h
        ${HDRDIR}options.h
        ${HDRDIR}optlist.h
        ${HDRDIR}paramsel.h
        ${HDRDIR}loocvprimal.h
        ${HDRDIR}fixlambda.h
        ${HDRDIR}optimization.h
        ${HDRDIR}rlsprimal.h
        ${HDRDIR}rlsdual.h
        ${HDRDIR}rlspegasos.h
        ${HDRDIR}pred.h
        ${HDRDIR}primal.h
        ${HDRDIR}dual.h
        ${HDRDIR}perf.h
        ${HDRDIR}precisionrecall.h
        ${HDRDIR}macroavg.h)

set(gurls_sources
        ${SRCDIR}kernel.cpp
        ${SRCDIR}norm.cpp
        ${SRCDIR}utils.cpp
        ${SRCDIR}exceptions.cpp
        ${SRCDIR}gvec.cpp
        ${SRCDIR}gmat2d.cpp
        ${SRCDIR}gmath.cpp
        ${SRCDIR}options.cpp
        ${SRCDIR}optlist.cpp
        ${SRCDIR}paramsel.cpp
        ${SRCDIR}loocvprimal.cpp
        ${SRCDIR}optimization.cpp
        ${SRCDIR}rlsprimal.cpp
        ${SRCDIR}pred.cpp
        ${SRCDIR}primal.cpp
        ${SRCDIR}dual.cpp
        ${SRCDIR}perf.cpp)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost 1.42 REQUIRED serialization)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories( ${Boost_LIBRARY_DIRS} )
endif(Boost_FOUND)

include_directories(${HDRDIR})

set(BUILD_SHARED_LIBS ON)

set(USE_BINARY_ARCHIVES OFF CACHE BOOL "If ON all the data structures in GURLS are saved/loaded using binary files.")

if(USE_BINARY_ARCHIVES)
        add_definitions(-DUSE_BINARY_ARCHIVES)
endif(USE_BINARY_ARCHIVES)

set(USE_MKL OFF CACHE BOOL "")

if(USE_MKL)

    SET(MKL_ARCH "em64t" CACHE STRING "System architecture (32, 64, em64t)")
    SET(MKLROOT "" CACHE PATH "Path to MKL")

    SET(MKL_LIB_PATH ${MKLROOT}/lib/${MKL_ARCH})

    include_directories(${MKLROOT}/include)
    include_directories(${MKL_LIB_PATH})

    if(variable STREQUAL "32")
        SET(ARCH_PREFIX "")
    else(variable STREQUAL "32")
        SET(ARCH_PREFIX "_lp64")
    endif(variable STREQUAL "32")

    SET(BLAS_LIBRARIES  "")

    SET(LAPACK_LIBRARIES
    ${MKL_LIB_PATH}/libmkl_scalapack${ARCH_PREFIX}.a
    ${MKL_LIB_PATH}/libmkl_solver${ARCH_PREFIX}.a
    -Wl,--start-group  ${MKL_LIB_PATH}/libmkl_intel${ARCH_PREFIX}.a ${MKL_LIB_PATH}/libmkl_gnu_thread.a ${MKL_LIB_PATH}/libmkl_core.a ${MKL_LIB_PATH}/libmkl_blacs_openmpi${ARCH_PREFIX}.a
    -Wl,--end-group -fopenmp -lpthread -lm
    )

endif(USE_MKL)

add_library(${GURLSLIBRARY} SHARED ${gurls_headers} ${gurls_sources})
target_link_libraries(${GURLSLIBRARY} gfortran pthread ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )


add_executable(gurlsTutorial ${TUTORIALDIR}examplegurls.cpp)
target_link_libraries(gurlsTutorial ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  pthread)

add_executable(gurlsIO ${TUTORIALDIR}IOexample.cpp)
target_link_libraries(gurlsIO ${GURLSLIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${Boost_LIBRARIES} pthread)

add_executable(gurlsSimpleExample ${TUTORIALDIR}example.cpp)
target_link_libraries(gurlsSimpleExample ${GURLSLIBRARY} pthread  ${Boost_LIBRARIES})

add_executable(gurlsHowToUseOPT ${TUTORIALDIR}howtouseopt.cpp)
target_link_libraries(gurlsHowToUseOPT ${GURLSLIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  pthread)

add_executable(gurlsBenchmark ${TUTORIALDIR}benchmark.cpp)
target_link_libraries(gurlsBenchmark ${GURLSLIBRARY} blas lapack pthread)
