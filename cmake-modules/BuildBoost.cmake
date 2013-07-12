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

############## Boost
if(MSVC)
    set(EXT "bat")
#address-model=64
    if(MSVC_VERSION GREATER 1699)
            set(CONF_TOOLSET "vc11")
            set(BUILD_TOOLSET "msvc-11.0")
    elseif(MSVC_VERSION GREATER 1599)
            set(CONF_TOOLSET "vc10")
            set(BUILD_TOOLSET "msvc-10.0")
    elseif(MSVC_VERSION GREATER 1499)
            set(CONF_TOOLSET "vc9")
            set(BUILD_TOOLSET "msvc-9.0")
    elseif(MSVC_VERSION GREATER 1399)
            set(CONF_TOOLSET "vc8")
            set(BUILD_TOOLSET "msvc-8.0")
    elseif(MSVC_VERSION GREATER 1299)
            set(CONF_TOOLSET "vc7")
            set(BUILD_TOOLSET "msvc-7.0")
    elseif(MSVC_VERSION GREATER 1199)
            message(FATAL_ERROR "vc6 is not supported by Boost: Please upgrade your compiler.")
    endif()
else()
    set(EXT "sh")
endif()

ExternalProject_add(buildBoost
    URL http://downloads.sourceforge.net/boost/boost_1_53_0.tar.gz
    URL_MD5 57a9e2047c0f511c4dfcf00eb5eb2fbb
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/boost
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.${EXT} ${CONF_TOOLSET}
    BUILD_COMMAND <SOURCE_DIR>/b2 toolset=${BUILD_TOOLSET} ${MAKE_ARGS} -d0 --layout=tagged variant=release link=static threading=multi runtime-link=shared
                --prefix=${EXTERNAL_PREFIX}
                --with-serialization
                --with-date_time
                --with-filesystem
                --with-test
                --with-system
                --with-signals
                install
    INSTALL_COMMAND ""
)

set_target_properties(buildBoost PROPERTIES EXCLUDE_FROM_ALL 1)
