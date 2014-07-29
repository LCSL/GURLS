# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
    message( FATAL_ERROR "You are trying to run CMake from the gurls++ directory, instead of just from the top directory")
endif()

set(GURLSLIBRARY gurls++)
project(${GURLSLIBRARY})

file(GLOB gurls_headers RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "include/gurls++/*.h*")
file(GLOB gurls_sources RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "src/*.cpp")

set(Gurls++_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "")
include_directories(${Gurls++_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${BLAS_LAPACK_INCLUDE_DIRS})
add_definitions(${BLAS_LAPACK_DEFINITIONS})
link_directories(${BLAS_LAPACK_LIBRARY_DIRS})

set (GurlsDependencies_LIBRARIES ${BLAS_LAPACK_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_DATE_TIME_LIBRARY})

add_library(${GURLSLIBRARY} ${GURLS_LIB_LINK} ${gurls_headers} ${gurls_sources} )

target_link_libraries(${GURLSLIBRARY} ${GurlsDependencies_LIBRARIES})

if(GURLS_USE_EXTERNAL_BLAS_LAPACK)
    add_dependencies(${GURLSLIBRARY} buildOpenblas)
endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)

if(GURLS_USE_EXTERNAL_BOOST)
    add_dependencies(${GURLSLIBRARY} buildBoost)
endif(GURLS_USE_EXTERNAL_BOOST)


install(TARGETS ${GURLSLIBRARY} EXPORT Gurlstargets-install
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install (FILES ${gurls_headers} DESTINATION include/gurls++)
install(EXPORT Gurlstargets-install DESTINATION lib/Gurls)

export(TARGETS gurls++ FILE ${gurls_BINARY_DIR}/Gurlstargets-buildtree.cmake)

if(MSVC)

    set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")

#    if(GURLS_BUILD_SHARED_LIBS)
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
#    else()
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_STATIC")
#    endif(GURLS_BUILD_SHARED_LIBS)

endif(MSVC)

set(Gurls++_LIBRARY ${GURLSLIBRARY} CACHE INTERNAL "")
set(Gurls++_LIBRARIES ${GURLSLIBRARY} ${GurlsDependencies_LIBRARIES} ) #to compile test executables

option(GURLSPP_BUILD_DEMO "" ON)
if(GURLSPP_BUILD_DEMO)
    add_subdirectory(demo)
endif(GURLSPP_BUILD_DEMO)

#option(GURLSPP_BUILD_TEST "" OFF)
#mark_as_advanced(FORCE GURLSPP_BUILD_TEST)
#if(GURLSPP_BUILD_TEST)
#     set(GURLSPP_DATA_DIR "" CACHE PATH "Path to the Gurls++ data directory")
#     mark_as_advanced(FORCE GURLSPP_DATA_DIR)
#
#     if(GURLSPP_DATA_DIR STREQUAL "")
#         message(WARNING, "GURLSPP_DATA_DIR value is set to an empty string. Please be sure this is the right configuration for you")
#     endif(GURLSPP_DATA_DIR STREQUAL "")
#
#     add_definitions(-DGURLS_DATA_DIR="${GURLSPP_DATA_DIR}")
#   add_subdirectory(test)
#   add_all_executables(${TESTDIR} ${GURLS_LINK_LIBRARIES})
#endif(GURLSPP_BUILD_TEST)

option(GURLSPP_BUILD_MISC "" OFF)
mark_as_advanced(FORCE GURLSPP_BUILD_MISC)
if(GURLSPP_BUILD_MISC)
    add_subdirectory(misc)
#    add_all_executables(${MISCDIR} ${GURLS_LINK_LIBRARIES})
endif(GURLSPP_BUILD_MISC)

# add a target to generate API documentation with Doxygen
option(GURLSPP_BUILD_DOC "Build Doxygen documentation" OFF)
if(GURLSPP_BUILD_DOC)
    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(gurlsppdoc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
    else(DOXYGEN_FOUND)
        message(WARNING "Doxygen documentation was enabled, but the Doxygen package was not found.")
    endif(DOXYGEN_FOUND)
endif(GURLSPP_BUILD_DOC)

