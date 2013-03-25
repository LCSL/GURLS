cmake_minimum_required(VERSION 2.6)

project(bgurls++)

set(GURLS_ROOT "" CACHE PATH "The root directory of gurls++")

option(USE_EXTERNALS "Link against external libraries" ON)
if(USE_EXTERNALS)
    set(EXTERNALS_PREFIX "" CACHE PATH "Directory where externals have been installed")
    if(EXTERNALS_PREFIX STREQUAL "")
        message(FATAL_ERROR "Please specify a valid path for EXTERNALS_PREFIX if you want to use externals")
    else()
        set(CMAKE_PREFIX_PATH ${EXTERNALS_PREFIX} ${CMAKE_PREFIX_PATH})
        set(Boost_NO_SYSTEM_PATHS ON)
    endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${GURLS_ROOT} ${CMAKE_MODULE_PATH})

if(CMAKE_COMPILER_IS_GNUCC OR MSVC)

    set(BLAS_LAPACK_IMPLEMENTATION "" CACHE STRING "Specify your blas/lapack implementation (MKL, ACML, ATLAS, NETLIB or an empty string)")

    if(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

        unset(BLAS_LIBRARIES CACHE)
        unset(LAPACK_LIBRARIES CACHE)

        find_package(MKL)

        include_directories(${MKL_INCLUDE_DIRS})
        link_directories(${MKL_LIBRARY_DIRS})
        add_definitions( ${MKL_DEFINITIONS})

        set(BLAS_LIBRARIES "")
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})

    elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ACML")

        unset(BLAS_LIBRARIES CACHE)
        unset(LAPACK_LIBRARIES CACHE)

        find_package(ACML)
        set(BLAS_LIBRARIES "")
        add_definitions( ${ACML_DEFINITIONS})
        set(LAPACK_LIBRARIES ${ACML_LIBRARIES})


    elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ATLAS")

        unset(BLAS_LIBRARIES CACHE)
        unset(LAPACK_LIBRARIES CACHE)

        find_package(ATLAS)
        set(BLAS_LIBRARIES ${ATLAS_BLAS_LIBS} ${ATLAS_LIBS} gfortran)
        add_definitions( ${ATLAS_DEFINITIONS})
        set(LAPACK_LIBRARIES ${ATLAS_LAPACK_LIBS})

    elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "NETLIB")

        unset(BLAS_LIBRARIES CACHE)
        unset(LAPACK_LIBRARIES CACHE)

        enable_language(Fortran)
        find_package(BLAS)
        find_package(LAPACK)

    else(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

        set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
        set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

    endif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

else(CMAKE_COMPILER_IS_GNUCC OR MSVC)

    set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
    set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

endif(CMAKE_COMPILER_IS_GNUCC OR MSVC)


set(GURLSLIBRARY bgurls++)
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/)
set(HDRDIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)
set(DEMODIR ${CMAKE_CURRENT_SOURCE_DIR}/demo/)
set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/test/)
set(MISCDIR ${CMAKE_CURRENT_SOURCE_DIR}/misc/)

set(GURLS_SRCDIR ${GURLS_ROOT}/src/)
set(GURLS_HDRDIR ${GURLS_ROOT}/include/)
#set(GURLS_DEMODIR ${GURLS_ROOT}/demo/)
#set(GURLS_TESTDIR ${GURLS_ROOT}/test/)
#set(GURLS_MISCDIR ${GURLS_ROOT}/misc/)

add_definitions( -D_BGURLS )

function(add_all_executables SRCDIR LINK_LIBS)

    file(GLOB sources ${SRCDIR}*.cpp)

    set(executables "")

    foreach(source_file ${sources})

        get_filename_component(executable "${source_file}" NAME_WE)

        set(executables ${executables} ${executable})

        add_executable(${executable} ${source_file} )
        target_link_libraries(${executable} ${LINK_LIBS})

        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH ${INSTALL_RPATH} ${CMAKE_INSTALL_PREFIX}/lib)

    endforeach(source_file)


    INSTALL(TARGETS ${executables}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
    )

endfunction(add_all_executables)


file(GLOB bgurls_headers RELATIVE ${PROJECT_SOURCE_DIR} ${HDRDIR}*.h*)
file(GLOB gurls_headers RELATIVE ${PROJECT_SOURCE_DIR} ${GURLS_HDRDIR}*.h*)
set(gurls_headers ${gurls_headers} ${bgurls_headers})

file(GLOB bgurls_sources RELATIVE ${PROJECT_SOURCE_DIR} ${SRCDIR}*.cpp)
file(GLOB gurls_sources RELATIVE ${PROJECT_SOURCE_DIR} ${GURLS_SRCDIR}*.cpp)
set(gurls_sources ${gurls_sources} ${bgurls_sources})

set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost 1.46.0 COMPONENTS serialization date_time  filesystem unit_test_framework system signals REQUIRED)

if(USE_EXTERNALS)
    unset(Boost_USE_STATIC_LIBS CACHE)
    set(Boost_USE_STATIC_LIBS ON)
else()
    option(Boost_USE_STATIC_LIBS "Link statically against boost libs" ON)
endif()

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)
mark_as_advanced(Boost_DIR)
include_directories( ${HDRDIR} ${GURLS_HDRDIR})



find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})


find_package(HDF5 COMPONENTS C REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})
add_definitions(${HDF5_DEFINITIONS})
mark_as_advanced(HDF5_DIR)


set(USE_BINARY_ARCHIVES ON CACHE BOOL "If ON all the data structures in bGURLS will be saved/loaded using binary files instead of text files")

if(USE_BINARY_ARCHIVES)
    add_definitions(-DUSE_BINARY_ARCHIVES)
endif(USE_BINARY_ARCHIVES)


if(USE_EXTERNALS)
    unset(BUILD_SHARED_LIBS CACHE)
    set(BUILD_SHARED_LIBS OFF)
else()
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build a shared library (ON) or a static one (OFF)")
endif()

if(BUILD_SHARED_LIBS)
    set( LIB_TYPE SHARED)
else(BUILD_SHARED_LIBS)
    set(LIB_TYPE STATIC)
endif(BUILD_SHARED_LIBS)


set ( ZLIB_LIBRARIES ${HDF5_z_LIBRARY})

add_definitions( -Wall )

add_library(${GURLSLIBRARY} ${LIB_TYPE} ${gurls_headers} ${gurls_sources} )
target_link_libraries(${GURLSLIBRARY} ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES})

INSTALL(TARGETS ${GURLSLIBRARY}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

if(MSVC)

    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")

    if(BUILD_SHARED_LIBS)
        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
    endif(BUILD_SHARED_LIBS)

endif(MSVC)


set(GURLS_LINK_LIBRARIES
    ${GURLSLIBRARY}
    ${Boost_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

set(BUILD_DEMO ON CACHE BOOL "")
if(BUILD_DEMO)

    #SET(MATIO_INCLUDE_DIR "" CACHE PATH "Libmatio include dir")
    #SET(MATIO_LIBRARIES_DIR "" CACHE PATH "Libmatio libraries dir")
    #set(GURLS_LINK_LIBRARIES_DEMO matio
    #	${GURLS_LINK_LIBRARIES}
    #)
    #include_directories( ${MATIO_INCLUDE_DIR})
    #link_directories(${MATIO_LIBRARIES_DIR})

    add_all_executables(${DEMODIR} ${GURLS_LINK_LIBRARIES})

endif(BUILD_DEMO)

#set(BUILD_TEST OFF CACHE BOOL "")
#if(BUILD_TEST)
#    set(BGURLS_DATA_DIR "" CACHE PATH "Path to the bGURLS++ data directory")
#
#    if(BGURLS_DATA_DIR STREQUAL "")
#        message(WARNING, "BGURLS_DATA_DIR value is set to an empty string. Please be sure this is the right configuration for you")
#    endif(BGURLS_DATA_DIR STREQUAL "")
#
#    add_definitions(-DGURLS_DATA_DIR="${BGURLS_DATA_DIR}")
#
#    add_all_executables(${TESTDIR} ${GURLS_LINK_LIBRARIES})
#endif(BUILD_TEST)

set(BUILD_MISC OFF CACHE BOOL "")
if(BUILD_MISC)
    add_all_executables(${MISCDIR} ${GURLS_LINK_LIBRARIES})
endif(BUILD_MISC)

# add a target to generate API documentation with Doxygen
set(BUILD_DOC OFF CACHE BOOL "Build Doxygen documentation")
if(BUILD_DOC)
    find_package(Doxygen REQUIRED)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(BUILD_DOC)
