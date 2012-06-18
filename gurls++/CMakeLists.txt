#
# Copyright (C) 2012 Matteo Santoro
# CopyPolicy: Released under the terms of GPL
#

project(gurls++)

set(GURLSLIBRARY gurls++)
set(SRCDIR src/)
set(HDRDIR include/)
set(TUTORIALDIR tutorial/)
set(TESTDIR test/)

set(gurls_headers "")

set(gurls_sources
        ${SRCDIR}gmath.cpp
        ${SRCDIR}options.cpp
        ${SRCDIR}optlist.cpp
)

set(USE_CUDA OFF CACHE BOOL "")

if(USE_CUDA)

    set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/cuda" ${CMAKE_MODULE_PATH})

    find_package(CUDA QUIET REQUIRED)

    if(NOT CUDA_FOUND)
        message("CUDA could not be located on your system.")
    endif(NOT CUDA_FOUND)

    set(gurls_sources
        ${gurls_sources}
        ${SRCDIR}gmath_gpu.cpp
    )

    set(GPULIBRARY ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

else(USE_CUDA)
    set(GPULIBRARY "")
endif(USE_CUDA)


set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost COMPONENTS serialization date_time REQUIRED)
OPTION(Boost_USE_STATIC_LIBS "Link statically to boost libs" OFF)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)

include_directories(
    ${HDRDIR}
)

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
    -Wl,--start-group
        ${MKL_LIB_PATH}/libmkl_intel${ARCH_PREFIX}.a
        ${MKL_LIB_PATH}/libmkl_gnu_thread.a
        ${MKL_LIB_PATH}/libmkl_core.a
        ${MKL_LIB_PATH}/libmkl_blacs_openmpi${ARCH_PREFIX}.a
    -Wl,--end-group -fopenmp -lpthread -lm
    )

#    SET(LAPACK_LIBRARIES
#    ${MKL_LIB_PATH}/libmkl_solver_lp64_sequential.a -Wl,--start-group  ${MKL_LIB_PATH}/libmkl_intel_lp64.a ${MKL_LIB_PATH}/libmkl_sequential.a ${MKL_LIB_PATH}/libmkl_core.a -Wl,--end-group -lpthread -lm
#    )

endif(USE_MKL)

set(USE_ACML OFF CACHE BOOL "")
if(USE_ACML)

    SET(ACMLROOT "" CACHE PATH "Path to ACML")

    SET(ACML_LIB_PATH ${ACMLROOT}/lib/)

    include_directories(${ACMLROOT}/include)

    SET(BLAS_LIBRARIES  "")

    SET(LAPACK_LIBRARIES
        ${ACML_LIB_PATH}/libacml_mp_dll.lib
    )

    add_definitions(-D_ACML)

endif(USE_ACML)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build a shared library (ON) or a static one (OFF)")
if(BUILD_SHARED_LIBS)
    set( LIB_TYPE SHARED)
else(BUILD_SHARED_LIBS)
    set(LIB_TYPE STATIC)
endif(BUILD_SHARED_LIBS)

SET(CBLAS_INCLUDE_PATH "" CACHE PATH "")
include_directories(${CBLAS_INCLUDE_PATH})

#ADD_DEFINITIONS( -Wall )
IF(WIN32)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)

    add_definitions(-DGOTOBLAS)
    set(FORTRAN "")
    set(THREADINGLIB "")

    SET(BLAS_LIBRARIES  "")
    SET(LAPACK_LIBRARIES "" CACHE FILEPATH "")

ELSEIF(WIN32)
    set(FORTRAN "gfortran")
    set(THREADINGLIB "pthread")
ENDIF(WIN32)

add_library(${GURLSLIBRARY} ${LIB_TYPE} ${gurls_headers} ${gurls_sources} )
target_link_libraries(${GURLSLIBRARY} ${FORTRAN} ${THREADINGLIB} ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GPULIBRARY})

IF(WIN32)
    set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")
    if(BUILD_SHARED_LIBS)
        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
    endif(BUILD_SHARED_LIBS)
ENDIF(WIN32)


add_executable(gurlsTutorial ${TUTORIALDIR}examplegurls.cpp )
target_link_libraries(gurlsTutorial ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

add_executable(gurlsSimpleExample ${TUTORIALDIR}example.cpp)
target_link_libraries(gurlsSimpleExample ${GURLSLIBRARY} ${THREADINGLIB}  ${Boost_LIBRARIES} ${GPULIBRARY})


set(BUILD_TEST OFF CACHE BOOL "")
if(BUILD_TEST)

    add_executable(GURLShoduallin ${TESTDIR}GURLShoduallin.cpp )
    target_link_libraries(GURLShoduallin ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLShoprimal ${TESTDIR}GURLShoprimal.cpp )
    target_link_libraries(GURLShoprimal ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLSloodualrbf ${TESTDIR}GURLSloodualrbf.cpp )
    target_link_libraries(GURLSloodualrbf ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLShodualrbf ${TESTDIR}GURLShodualrbf.cpp )
    target_link_libraries(GURLShodualrbf ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLSlooduallin ${TESTDIR}GURLSlooduallin.cpp )
    target_link_libraries(GURLSlooduallin ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLSlooprimal ${TESTDIR}GURLSlooprimal.cpp )
    target_link_libraries(GURLSlooprimal ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

    add_executable(GURLSpegasos ${TESTDIR}GURLSpegasos.cpp )
    target_link_libraries(GURLSpegasos ${GURLSLIBRARY} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}  ${THREADINGLIB} ${GPULIBRARY})

endif(BUILD_TEST)

