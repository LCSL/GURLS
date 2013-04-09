# - Find MKL library
# This module finds an installed MKL library
#
# This module sets the following variables:
#  MKL_INCLUDE_DIRS - uncached list of include directories
#  MKL_LIBRARY_DIRS - uncached list of library directories
#
#  MKL_DEFINITIONS - uncached list of required definitions
#
#  MKL_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use MKL

SET(MKL_ARCH "em64t" CACHE STRING "System architecture (32, 64, em64t)")
SET(MKL_ROOT "" CACHE PATH "Path to MKL")

SET(MKL_LIBRARY_DIRS ${MKL_ROOT}/lib/${MKL_ARCH})
IF(MKL_ARCH STREQUAL "em64t")
    SET(MKL_LIBRARY_DIRS ${MKL_LIBRARY_DIRS} ${MKL_ROOT}/lib/intel64)
ENDIF(MKL_ARCH STREQUAL "em64t")

SET(MKL_INCLUDE_DIRS ${MKL_ROOT}/include)


#set(MKL_INTEGER64 OFF CACHE BOOL "Enable 64 bit integers on MKL")
set(MKL_INTEGER64 OFF)

if(variable STREQUAL "32")
    SET(ARCH_PREFIX "")
else(variable STREQUAL "32")
    if(MKL_INTEGER64)
        SET(ARCH_PREFIX "_ilp64")
    else(MKL_INTEGER64)
        SET(ARCH_PREFIX "_lp64")
    endif(MKL_INTEGER64)
endif(variable STREQUAL "32")

set(MKL_MULTITHREADING ON CACHE BOOL "Enable multithreaded libraries")
#set(MKL_STATIC ON CACHE BOOL "Link statically against MKL libraries")

if(MSVC)
    set(MKL_IOMP ON)
else(MSVC)
    set(MKL_IOMP ON CACHE BOOL "Use intel OpenMP implementation")
endif(MSVC)


if(MSVC)
    SET(LIB_SUFFIX "_dll.lib")
    SET(OMP_SUFFIX "md.lib")
else(MSVC)
    SET(LIB_SUFFIX "")
    SET(OMP_SUFFIX "")
endif(MSVC)

if(MKL_MULTITHREADING)
    SET(SEQUENTIAL "")
    if(MKL_IOMP)
        SET(MULTITHREADING_LIBS
            mkl_intel_thread${LIB_SUFFIX}
            iomp5${OMP_SUFFIX}
            )
        SET(MKL_IOMP_LIB_DIR "" CACHE PATH "Path to the intel OpenMP library")
        SET(MKL_LIBRARY_DIRS ${MKL_LIBRARY_DIRS} ${IOMP_LIB_DIR})
    else(MKL_IOMP)
        SET(MULTITHREADING_LIBS
            mkl_gnu_thread${LIB_SUFFIX}
            -fopenmp
            )
    endif(MKL_IOMP)
else(MKL_MULTITHREADING)
    SET(SEQUENTIAL "_sequential")
    SET( MULTITHREADING_LIBS mkl_sequential${LIB_SUFFIX} )
endif(MKL_MULTITHREADING)


if(MSVC)
    SET(MKL_SOLVER mkl_solver${ARCH_PREFIX}${SEQUENTIAL}.lib)
else(MSVC)
    SET(MKL_SOLVER "")
endif(MSVC)

SET(MKL_LIBRARIES
    mkl_intel${ARCH_PREFIX}${LIB_SUFFIX}
    mkl_core${LIB_SUFFIX}
    ${MULTITHREADING_LIBS}
    ${MKL_SOLVER}
)

#SET(MKL_LIBRARIES
#    mkl_rt
#)

if(NOT MSVC)
    if(MKL_INTEGER64)
        ADD_DEFINITIONS(-DMKL_ILP64)
    endif(MKL_INTEGER64)

    SET(MKL_LIBRARIES
    ${MKL_LIBRARIES}
    -lpthread
    -lm
    )
endif(NOT MSVC)
