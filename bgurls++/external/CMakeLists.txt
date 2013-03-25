
cmake_minimum_required(VERSION 2.6)

if(NOT UNIX)
    message(FATAL_ERROR "Externals build currently works only on unix systems")
endif()

include(ExternalProject)


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ CACHE PATH "Install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(EXTERNAL_PREFIX ${CMAKE_INSTALL_PREFIX})


if( ${CMAKE_VERSION} VERSION_GREATER "2.8.5" OR ${CMAKE_VERSION} VERSION_EQUAL "2.8.5")
    include(ProcessorCount)
    ProcessorCount(N)
    if(${N} GREATER 1)
        set(MAKE_ARGS -j${N})
    endif()
endif()



############## Zlib

ExternalProject_add(zlib
    URL http://zlib.net/zlib-1.2.7.tar.gz
    URL_MD5 60df6a37c56e7c1366cca812414f7b85
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/zlib
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR> --static
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
)

#set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL 1)



############## MPICH2

ExternalProject_add(MPICH2
    URL http://www.mpich.org/static/tarballs/3.0.2/mpich-3.0.2.tar.gz
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/MPICH2
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND CFLAGS=-w CXXFLAGS=-w <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR> --disable-f77 --disable-fc --enable-g=none --enable-romio --enable-shared --disable-static
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
)

#set_target_properties(MPICH2 PROPERTIES EXCLUDE_FROM_ALL 1)



############## hdf5

ExternalProject_add(hdf5
    URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/hdf5
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND CC=${EXTERNAL_PREFIX}/bin/mpicc CFLAGS=-w <SOURCE_DIR>/configure
    --prefix=<INSTALL_DIR> --enable-production --enable-parallel --disable-shared --disable-fortran --disable-hl --with-zlib=${EXTERNAL_PREFIX}/
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
)
add_dependencies(hdf5 zlib MPICH2)

#set_target_properties(hdf5 PROPERTIES EXCLUDE_FROM_ALL 1)



############## Boost

ExternalProject_add(boost
    URL http://downloads.sourceforge.net/boost/boost_1_53_0.tar.gz
    URL_MD5 57a9e2047c0f511c4dfcf00eb5eb2fbb
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/boost
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.sh --prefix=${EXTERNAL_PREFIX} --with-libraries=serialization,date_time,filesystem,test,system,signals
    BUILD_COMMAND <SOURCE_DIR>/b2 ${MAKE_ARGS} -d0 --layout=tagged variant=release link=static threading=multi runtime-link=shared install
    INSTALL_COMMAND ""
)

#set_target_properties(boost PROPERTIES EXCLUDE_FROM_ALL 1)


#add_custom_target(External
#  DEPENDS netCDF hdf5 zlib MPICH2
#  COMMENT "Building external targets"
#)
