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
    set(LAYOUT "versioned")
    set(EXT "bat")
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
    set(BUILD_TOOLSET "toolset=${BUILD_TOOLSET}")
    
    if (CMAKE_CL_64)
        set(ADDRESS_MODEL "address-model=64")
    else()
        set(ADDRESS_MODEL "address-model=32")
    endif()
	set(Boost_VARIANT "debug" CACHE STRING "Possible values: debug, release")

else()
    set(LAYOUT "tagged")
    set(BUILD_TOOLSET "")
    set(EXT "sh")
    set(ADDRESS_MODEL "") # for linux??
	set(Boost_VARIANT "debug,release" CACHE STRING "Possible values: debug, release")

endif()


ExternalProject_add(buildBoost
    URL http://downloads.sourceforge.net/boost/boost_1_55_0.tar.gz
    URL_MD5 93780777cfbf999a600f62883bd54b17
    BUILD_IN_SOURCE 1
    SOURCE_DIR  ${EXTERNAL_PREFIX}/src/boost
    INSTALL_DIR ${EXTERNAL_PREFIX}/
    CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.${EXT} ${CONF_TOOLSET}
    BUILD_COMMAND <SOURCE_DIR>/b2 ${BUILD_TOOLSET} ${MAKE_ARGS} ${ADDRESS_MODEL} 
                variant=${Boost_VARIANT}
                link=static 
                threading=multi 
                runtime-link=shared
                -d0 
                --layout=${LAYOUT}
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

if(NOT MSVC)

    set(Boost_SERIALIZATION_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_serialization-mt.a)

    set(Boost_DATE_TIME_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_date_time-mt.a)

    set(Boost_FILESYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_filesystem-mt.a)

    set(Boost_UNIT_TEST_FRAMEWORK_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_unit_test_framework-mt.a)

    set(Boost_SYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_system-mt.a)

    set(Boost_SIGNALS_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_signals-mt.a)

set(Boost_INCLUDE_DIR ${EXTERNAL_PREFIX}/include)

else()

set(TAG_TOOLSET ${CONF_TOOLSET}0)
if(${Boost_VARIANT} STREQUAL release)
    set(Boost_SERIALIZATION_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_serialization-${TAG_TOOLSET}-mt-1_55.lib)

    set(Boost_DATE_TIME_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_date_time-${TAG_TOOLSET}-mt-1_55.lib)

    set(Boost_FILESYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_filesystem-${TAG_TOOLSET}-mt-1_55.lib)

    set(Boost_UNIT_TEST_FRAMEWORK_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_unit_test_framework-${TAG_TOOLSET}-mt-1_55.lib)

    set(Boost_SYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_system-${TAG_TOOLSET}-mt-1_55.lib)

    set(Boost_SIGNALS_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_signals-${TAG_TOOLSET}-mt-1_55.lib)
else()
    set(Boost_SERIALIZATION_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_serialization-${TAG_TOOLSET}-mt-gd-1_55.lib)

    set(Boost_DATE_TIME_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_date_time-${TAG_TOOLSET}-mt-gd-1_55.lib)

    set(Boost_FILESYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_filesystem-${TAG_TOOLSET}-mt-gd-1_55.lib)

    set(Boost_UNIT_TEST_FRAMEWORK_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_unit_test_framework-${TAG_TOOLSET}-mt-gd-1_55.lib)

    set(Boost_SYSTEM_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_system-${TAG_TOOLSET}-mt-gd-1_55.lib)

    set(Boost_SIGNALS_LIBRARY ${EXTERNAL_PREFIX}/lib/libboost_signals-${TAG_TOOLSET}-mt-gd-1_55.lib)
endif()
set(Boost_INCLUDE_DIRS ${EXTERNAL_PREFIX}/include/boost-1_55)

set(Boost_INCLUDE_DIR ${EXTERNAL_PREFIX}/include)


endif()



set(Boost_FOUND 1)
set_target_properties(buildBoost PROPERTIES EXCLUDE_FROM_ALL 1)

