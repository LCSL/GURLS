# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>


#Check if BLAS and LAPACK are already available on the system

if(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

    find_package(MKL)

    set(BLAS_LAPACK_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
    set(BLAS_LAPACK_LIBRARY_DIRS ${MKL_LIBRARY_DIRS})
    set(BLAS_LAPACK_DEFINITIONS ${MKL_DEFINITIONS})
    set(BLAS_LAPACK_LIBRARIES ${MKL_LIBRARIES})
    set(BLAS_LAPACK_FOUND ${MKL_FOUND})

elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ACML")

    find_package(ACML)
    set(BLAS_LAPACK_INCLUDE_DIRS )
    set(BLAS_LAPACK_LIBRARY_DIRS )
    set(BLAS_LAPACK_DEFINITIONS ${ACML_DEFINITIONS})
    set(BLAS_LAPACK_LIBRARIES ${ACML_LIBRARIES})
    set(BLAS_LAPACK_FOUND ${ACML_FOUND})

elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ATLAS")

    find_package(ATLAS)
    set(BLAS_LAPACK_INCLUDE_DIRS )
    set(BLAS_LAPACK_LIBRARY_DIRS )
    set(BLAS_LAPACK_DEFINITIONS ${ATLAS_DEFINITIONS})
    set(BLAS_LAPACK_LIBRARIES ${ATLAS_LAPACK_LIBS} ${ATLAS_BLAS_LIBS} ${ATLAS_LIBS} gfortran)
    set(BLAS_LAPACK_FOUND ${ATLAS_FOUND})

elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "OPENBLAS")

    if(CMAKE_COMPILER_IS_GNUCC)
            enable_language(Fortran)
    endif()
    set (OPENBLAS_IGNORE_HEADERS ON)
    find_package(Openblas)

    set(BLAS_LAPACK_INCLUDE_DIRS ${Openblas_INCLUDE_DIRS})
    set(BLAS_LAPACK_LIBRARY_DIRS )
    set(BLAS_LAPACK_DEFINITIONS )
    set(BLAS_LAPACK_LIBRARIES ${Openblas_LIBRARIES})
    set(BLAS_LAPACK_FOUND ${Openblas_FOUND})

else(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")
# includes case "NETLIB"

    if(CMAKE_COMPILER_IS_GNUCC)
            enable_language(Fortran)
    endif()
    find_package(BLAS)
    find_package(LAPACK)

    set(BLAS_LAPACK_INCLUDE_DIRS )
    set(BLAS_LAPACK_LIBRARY_DIRS )
    set(BLAS_LAPACK_DEFINITIONS )
    set(BLAS_LAPACK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(${BLAS_FOUND} AND ${LAPACK_FOUND})
        set(BLAS_LAPACK_FOUND TRUE)
    else()
        set(BLAS_LAPACK_FOUND FALSE)
    endif()


endif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

if(NOT BLAS_LAPACK_FOUND)
    option(GURLS_USE_EXTERNAL_BLAS_LAPACK "build external project Openblas" OFF)
endif()

if(GURLS_USE_EXTERNAL_BLAS_LAPACK)
    include(BuildOpenblas)
endif()
#################### BOOST


set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
option(Boost_USE_STATIC_LIBS "Link statically against boost libs" ON)

#should not be needed #set(CMAKE_PREFIX_PATH $ENV{GURLSPP_ROOT} ${CMAKE_PREFIX_PATH})
find_package( Boost ${BOOST_MINIMUM_VERSION} COMPONENTS serialization date_time filesystem unit_test_framework system signals REQUIRED)
mark_as_advanced(Boost_DIR)

    if(MSVC)
        add_definitions(-DBOOST_ALL_NO_LIB)
    endif()

if(NOT Boost_FOUND)
    option(GURLS_USE_EXTERNAL_BOOST "build external project Boost" OFF)
endif()

if(GURLS_USE_EXTERNAL_BOOST) #if GURLS_USE_EXTERNAL_BOOST was enabled by the user...
    include(BuildBoost)
endif()


# HDF_5
find_package(HDF5 COMPONENTS C)
if (NOT HDF5_FOUND)
    option(GURLS_USE_EXTERNAL_HDF5 "build external project HDF5" OFF)
endif()

if(GURLS_USE_EXTERNAL_HDF5) #if GURLS_USE_EXTERNAL_HDF5 was enabled by the user...
    include(BuildHDF5)
endif()

