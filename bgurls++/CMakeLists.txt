# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
    message( FATAL_ERROR "You are trying to run CMake from the bgurls++ directory, instead of just from the top directory")
endif()

set(GURLSLIBRARY bgurls++)
project(${GURLSLIBRARY})

add_definitions( -D_BGURLS )
set(export_definitions ${export_definitions} -D_BGURLS)
set(bgurls_headers       include/bgurls++/bgurls.h
                        include/bgurls++/bigarray.h
                        include/bgurls++/bigarray.hpp
                        include/bgurls++/bigmath.h
                        include/bgurls++/bigoptimization.h
                        include/bgurls++/bigoptimizer_rlspegasos.h
                        include/bgurls++/bigoptimizer_rlsprimal.h
                        include/bgurls++/bigoptlist.h
                        include/bgurls++/bigparamsel_calibratesgd.h
                        include/bgurls++/bigparamsel.h
                        include/bgurls++/bigparamsel_hoprimal.h
                        include/bgurls++/bigperf.h
                        include/bgurls++/bigperf_macroavg.h
                        include/bgurls++/bigpred.h
                        include/bgurls++/bigpred_primal.h
                        include/bgurls++/bigsplit.h
                        include/bgurls++/bigsplit_ho.h
                        include/bgurls++/mpi_utils.h
    )

set(bgurls_sources       src/bigarray.cpp
                        src/bigmath.cpp
                        src/bigoptlist.cpp
                        src/mpi_utils.cpp
    )

set(gurls_sources   ../gurls++/src/blas_lapack.cpp
                ../gurls++/src/gmath.cpp
                ../gurls++/src/optarray.cpp
                ../gurls++/src/optfunction.cpp
                ../gurls++/src/options.cpp
                ../gurls++/src/optlist.cpp
                ../gurls++/src/optmatrix.cpp
    )

set(BGurls++_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/ CACHE INTERNAL "" )

include_directories( ${BGurls++_INCLUDE_DIRS} ${Gurls++_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${BLAS_LAPACK_INCLUDE_DIRS})
add_definitions(${BLAS_LAPACK_DEFINITIONS})
link_directories(${BLAS_LAPACK_LIBRARY_DIRS})
add_definitions(${HDF5_DEFINITIONS})
set(export_definitions ${export_definitions} ${HDF5_DEFINITIONS}) #TODO check if needed here

#TODO install headers

# if(MPI_FOUND)
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
# endif(MPI_FOUND)


option(BGURLSPP_USE_MPI_IO "Access HDF5 files using MPI_IO instead of MPI_POSIX" OFF)
if(BGURLSPP_USE_MPI_IO)
    add_definitions(-DUSE_MPIIO)
    set(export_definitions ${export_definitions} -DUSE_MPIIO)
endif(BGURLSPP_USE_MPI_IO)

add_library(${GURLSLIBRARY} ${GURLS_LIB_LINK} ${bgurls_headers} ${bgurls_sources} ${gurls_sources})
target_link_libraries(${GURLSLIBRARY} ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_SIGNALS_LIBRARY})
set(BGurls++_LIBRARY ${GURLSLIBRARY}  CACHE INTERNAL "")
set(BGurls++_LIBRARIES ${GURLSLIBRARY} ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_SIGNALS_LIBRARY}) #to compile executables - maybe missing dependencies though

if(GURLS_USE_EXTERNAL_BOOST)
    add_dependencies(${GURLSLIBRARY} buildBoost)
endif(GURLS_USE_EXTERNAL_BOOST)
if(GURLS_USE_EXTERNAL_HDF5)
    add_dependencies(${GURLSLIBRARY} buildHdf5)
endif(GURLS_USE_EXTERNAL_HDF5)

INSTALL(TARGETS ${GURLSLIBRARY}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

#if(MSVC)

    #set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")

#    if(GURLS_BUILD_SHARED_LIBS)
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
#    else()
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_STATIC")
#    endif(GURLS_BUILD_SHARED_LIBS)

#endif(MSVC)


# set(GURLS_LINK_LIBRARIES
#     ${GURLSLIBRARY}
#     gurls++ #mancano le dipendenze di gurls++
#     ${MPI_LIBRARIES}
#     ${HDF5_LIBRARIES}
#     ${ZLIB_LIBRARIES}
# )

option(BGURLSPP_BUILD_DEMO "Build bGurls++ demos" ON)
if(BGURLSPP_BUILD_DEMO)

    #SET(MATIO_INCLUDE_DIR "" CACHE PATH "Libmatio include dir")
    #SET(MATIO_LIBRARIES_DIR "" CACHE PATH "Libmatio libraries dir")
    #set(GURLS_LINK_LIBRARIES_DEMO matio
    #	${GURLS_LINK_LIBRARIES}
    #)
    #include_directories( ${MATIO_INCLUDE_DIR})
    #link_directories(${MATIO_LIBRARIES_DIR})

    #add_all_executables(${DEMODIR} ${GURLS_LINK_LIBRARIES})
    add_subdirectory(demo)
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

option(BGURLSPP_BUILD_MISC "" OFF)
mark_as_advanced(FORCE BGURLSPP_BUILD_MISC)
if(BGURLSPP_BUILD_MISC)
    add_subdirectory(misc)
endif(BGURLSPP_BUILD_MISC)

# add a target to generate API documentation with Doxygen
option(BGURLSPP_BUILD_DOC "Build Doxygen documentation" OFF)
if(BGURLSPP_BUILD_DOC)
    find_package(Doxygen REQUIRED)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(bgurlsppdoc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(BGURLSPP_BUILD_DOC)
