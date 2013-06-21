cmake_minimum_required(VERSION 2.6)

if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
    message( FATAL_ERROR "You are trying to run CMake from the bgurls++ directory, instead of just from the top directory")
endif()

project(bgurls++)

get_filename_component(GURLS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../gurls++" ABSOLUTE)

set(GURLSLIBRARY bgurls++)
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/)
set(HDRDIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)
set(DEMODIR ${CMAKE_CURRENT_SOURCE_DIR}/demo/)
set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/test/)
set(MISCDIR ${CMAKE_CURRENT_SOURCE_DIR}/misc/)

set(GURLS_SRCDIR ${GURLS_ROOT}/src/)
set(GURLS_HDRDIR ${GURLS_ROOT}/include/)

add_definitions( -D_BGURLS )


file(GLOB bgurls_headers RELATIVE ${PROJECT_SOURCE_DIR} ${HDRDIR}*.h*)
file(GLOB gurls_headers RELATIVE ${PROJECT_SOURCE_DIR} ${GURLS_HDRDIR}*.h*)
set(gurls_headers ${gurls_headers} ${bgurls_headers})

file(GLOB bgurls_sources RELATIVE ${PROJECT_SOURCE_DIR} ${SRCDIR}*.cpp)
file(GLOB gurls_sources RELATIVE ${PROJECT_SOURCE_DIR} ${GURLS_SRCDIR}*.cpp)
set(gurls_sources ${gurls_sources} ${bgurls_sources})

include_directories( ${HDRDIR} ${GURLS_HDRDIR})


if(GURLS_USE_EXTERNAL_HDF5)

    if(WIN32)
      set(prefix "")
      set(suffix ".lib")
    else()
      set(prefix "lib")
      set(suffix ".a")
    endif()

    set(MPI_COMPILE_FLAGS "")
    set(MPI_LINK_FLAGS "-Wl,-rpath  -Wl,${GURLS_EXTERNALS_PREFIX}/lib")
    set(MPI_INCLUDE_PATH ${GURLS_EXTERNALS_PREFIX}/include)

    set(HDF5_INCLUDE_DIRS ${GURLS_EXTERNALS_PREFIX}/include)
    set(HDF5_LIBRARY_DIRS ${GURLS_EXTERNALS_PREFIX}/lib)
    set(HDF5_DEFINITIONS -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE)

    set(ZLIB_INCLUDE_DIRS ${GURLS_EXTERNALS_PREFIX}/include)
    set(ZLIB_LIBRARY_DIRS ${GURLS_EXTERNALS_PREFIX}/lib)
    set(ZLIB_DEFINITIONS "")

    set(MPI_LIBRARIES mpichcxx mpich opa mpl rt pthread)
    set(HDF5_LIBRARIES m ${prefix}hdf5${suffix})
    if(MSVC)
        set(ZLIB_LIBRARIES zlibstatic.lib)
    else()
        set(ZLIB_LIBRARIES libz.a)
    endif()

    set(MPI_FOUND ON)
    set(HDF5_FOUND ON)
    set(ZLIB_FOUND ON)

else(GURLS_USE_EXTERNAL_HDF5)

    find_package(MPI REQUIRED)

    find_package(HDF5 COMPONENTS C REQUIRED)
    mark_as_advanced(HDF5_DIR)

    find_package(ZLIB REQUIRED)

    #set ( ZLIB_LIBRARIES ${HDF5_z_LIBRARY})

endif(GURLS_USE_EXTERNAL_HDF5)


#message("MPI_COMPILE_FLAGS=${MPI_COMPILE_FLAGS}")
#message("MPI_LINK_FLAGS=${MPI_LINK_FLAGS}")
#message("MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}")

#message("HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
#message("HDF5_LIBRARY_DIRS=${HDF5_LIBRARY_DIRS}")
#message("HDF5_DEFINITIONS=${HDF5_DEFINITIONS}")

#message("ZLIB_INCLUDE_DIRS=${ZLIB_INCLUDE_DIRS}")
#message("ZLIB_LIBRARY_DIRS=${ZLIB_LIBRARY_DIRS}")
#message("ZLIB_DEFINITIONS=${ZLIB_DEFINITIONS}")

#message("MPI_LIBRARIES=${MPI_LIBRARIES}")
#message("HDF5_LIBRARIES=${HDF5_LIBRARIES}")
#message("ZLIB_LIBRARIES=${ZLIB_LIBRARIES}")


if(MPI_FOUND)
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    link_directories(${HDF5_LIBRARY_DIRS})
    add_definitions(${HDF5_DEFINITIONS})
endif(HDF5_FOUND)

if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARY_DIRS})
    add_definitions(${ZLIB_DEFINITIONS})
endif(ZLIB_FOUND)


option(BGURLSPP_USE_MPI_IO "Access HDF5 files using MPI_IO instead of MPI_POSIX" OFF)
if(BGURLSPP_USE_MPI_IO)
    add_definitions(-DUSE_MPIIO)
endif(BGURLSPP_USE_MPI_IO)


add_library(${GURLSLIBRARY} ${GURLS_LIB_LINK} ${gurls_headers} ${gurls_sources} )
target_link_libraries(${GURLSLIBRARY} ${GurlsDependencies_LIBRARIES} ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES})

if(GURLS_USE_EXTERNAL_BOOST)
    add_dependencies(${GURLSLIBRARY} boost)
endif(GURLS_USE_EXTERNAL_BOOST)
if(GURLS_USE_EXTERNAL_HDF5)
    add_dependencies(${GURLSLIBRARY} hdf5)
endif(GURLS_USE_EXTERNAL_HDF5)

INSTALL(TARGETS ${GURLSLIBRARY}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

if(MSVC)

    #set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")

    if(GURLS_BUILD_SHARED_LIBS)
        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
#    else()
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_STATIC")
    endif(GURLS_BUILD_SHARED_LIBS)

endif(MSVC)


set(GURLS_LINK_LIBRARIES
    ${GURLSLIBRARY}
    ${GurlsDependencies_LIBRARIES}
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

set(BGURLSPP_BUILD_DEMO ON CACHE BOOL "")
if(BGURLSPP_BUILD_DEMO)

    #SET(MATIO_INCLUDE_DIR "" CACHE PATH "Libmatio include dir")
    #SET(MATIO_LIBRARIES_DIR "" CACHE PATH "Libmatio libraries dir")
    #set(GURLS_LINK_LIBRARIES_DEMO matio
    #	${GURLS_LINK_LIBRARIES}
    #)
    #include_directories( ${MATIO_INCLUDE_DIR})
    #link_directories(${MATIO_LIBRARIES_DIR})

    add_all_executables(${DEMODIR} ${GURLS_LINK_LIBRARIES})

endif(BGURLSPP_BUILD_DEMO)

#set(BGURLSPP_BUILD_TEST OFF CACHE BOOL "")
#mark_as_advanced(FORCE BGURLSPP_BUILD_TEST)
#if(BGURLSPP_BUILD_TEST)
#    set(BGURLSPP_DATA_DIR "" CACHE PATH "Path to the bGURLS++ data directory")
#
#    if(BGURLSPP_DATA_DIR STREQUAL "")
#        message(WARNING, "BGURLSPP_DATA_DIR value is set to an empty string. Please be sure this is the right configuration for you")
#    endif(BGURLSPP_DATA_DIR STREQUAL "")
#
#    add_definitions(-DGURLS_DATA_DIR="${BGURLSPP_DATA_DIR}")
#
#    add_all_executables(${TESTDIR} ${GURLS_LINK_LIBRARIES})
#endif(BGURLSPP_BUILD_TEST)

set(BGURLSPP_BUILD_MISC OFF CACHE BOOL "")
mark_as_advanced(FORCE BGURLSPP_BUILD_MISC)
if(BGURLSPP_BUILD_MISC)
    add_all_executables(${MISCDIR} ${GURLS_LINK_LIBRARIES})
endif(BGURLSPP_BUILD_MISC)

# add a target to generate API documentation with Doxygen
set(BGURLSPP_BUILD_DOC OFF CACHE BOOL "Build Doxygen documentation")
if(BGURLSPP_BUILD_DOC)
    find_package(Doxygen REQUIRED)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(bgurlsppdoc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(BGURLSPP_BUILD_DOC)
