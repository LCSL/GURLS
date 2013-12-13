# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if( ${CMAKE_VERSION} VERSION_GREATER "2.8.5" OR ${CMAKE_VERSION} VERSION_EQUAL "2.8.5")
    include(ProcessorCount)
    ProcessorCount(N)
    if(${N} GREATER 1)
        set(MAKE_ARGS -j${N})
    endif()
endif()

include(ExternalProject)

if (MSVC)
    message(WARNING "Cannot build external project for OpenBLAS package on Windows")
    set(Openblas_FOUND FALSE)
else()
    ExternalProject_add(buildOpenblas
        GIT_REPOSITORY "git://github.com/xianyi/OpenBLAS.git"
        GIT_TAG "v0.2.8"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        SOURCE_DIR  ${EXTERNAL_PREFIX}/src/openblas
        INSTALL_DIR ${EXTERNAL_PREFIX}
        CONFIGURE_COMMAND ""
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_LIST_DIR}/../external/Makefile.install.openblas" "<SOURCE_DIR>/Makefile.install"
        BUILD_COMMAND NO_SHARED=1 ${CMAKE_MAKE_PROGRAM} libs netlib
        INSTALL_COMMAND NO_SHARED=1 ${CMAKE_MAKE_PROGRAM} install PREFIX=<INSTALL_DIR>
    )

    set(Openblas_FOUND TRUE)
    #set(Openblas_LIBRARIES ${Openblas_LIBRARY})
    set(BLAS_LAPACK_LIBRARIES ${EXTERNAL_PREFIX}/lib/libopenblas_sandybridgep-r0.2.8.a)
	if( CMAKE_COMPILER_IS_GNUCC)
		list(APPEND BLAS_LAPACK_LIBRARIES gfortran pthread)
	endif()
    set(BLAS_LAPACK_FOUND ${Openblas_FOUND})
    set(BLAS_LAPACK_INCLUDE_DIRS  ${EXTERNAL_PREFIX}/include)
    set_target_properties(buildOpenblas PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

