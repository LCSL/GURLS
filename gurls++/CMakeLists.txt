cmake_minimum_required(VERSION 2.6)

project(gurls++)

set(GURLSLIBRARY gurls++)
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/)
set(HDRDIR ${CMAKE_CURRENT_SOURCE_DIR}/include/)
set(DEMODIR ${CMAKE_CURRENT_SOURCE_DIR}/demo/)
set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/test/)
set(MISCDIR ${CMAKE_CURRENT_SOURCE_DIR}/misc/)


file(GLOB gurls_headers ${HDRDIR}*.h*)
file(GLOB gurls_sources ${SRCDIR}*.cpp)

include_directories( BEFORE ${HDRDIR} )


add_library(${GURLSLIBRARY} ${GURLS_LIB_LINK} ${gurls_headers} ${gurls_sources} )
target_link_libraries(${GURLSLIBRARY} ${GurlsDependencies_LIBRARIES})

if(GURLS_USE_EXTERNAL_BLAS_LAPACK)
    add_dependencies(${GURLSLIBRARY} openblas)
endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)

if(GURLS_USE_EXTERNAL_BOOST)
    add_dependencies(${GURLSLIBRARY} boost)
endif(GURLS_USE_EXTERNAL_BOOST)

INSTALL(TARGETS ${GURLSLIBRARY}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
INSTALL (FILES ${gurls_headers} DESTINATION include/gurls++)

if(MSVC)

    set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_FLAGS "/EHa")

    if(GURLS_BUILD_SHARED_LIBS)
        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_EXPORTS")
#    else()
#        set_target_properties(${GURLSLIBRARY} PROPERTIES COMPILE_DEFINITIONS "_GURLS_STATIC")
    endif(GURLS_BUILD_SHARED_LIBS)

endif(MSVC)

set(GURLS_LINK_LIBRARIES
    ${GURLSLIBRARY}
    ${GurlsDependencies_LIBRARIES}
)

set(GURLSPP_BUILD_DEMO OFF CACHE BOOL "")
if(GURLSPP_BUILD_DEMO)
    add_all_executables(${DEMODIR} ${GURLS_LINK_LIBRARIES})
endif(GURLSPP_BUILD_DEMO)

set(GURLSPP_BUILD_TEST OFF CACHE BOOL "")
if(GURLSPP_BUILD_TEST)
    set(GURLS_DATA_DIR "" CACHE PATH "Path to the Gurls++ data directory")

    if(GURLSPP_DATA_DIR STREQUAL "")
        message(WARNING, "GURLSPP_DATA_DIR value is set to an empty string. Please be sure this is the right configuration for you")
    endif(GURLSPP_DATA_DIR STREQUAL "")

    add_definitions(-DGURLS_DATA_DIR="${GURLSPP_DATA_DIR}")

    add_all_executables(${TESTDIR} ${GURLS_LINK_LIBRARIES})
endif(GURLSPP_BUILD_TEST)

set(GURLSPP_BUILD_MISC OFF CACHE BOOL "")
if(GURLSPP_BUILD_MISC)
    add_all_executables(${MISCDIR} ${GURLS_LINK_LIBRARIES})
endif(GURLSPP_BUILD_MISC)

# add a target to generate API documentation with Doxygen
set(GURLSPP_BUILD_DOC OFF CACHE BOOL "Build Doxygen documentation")
if(GURLSPP_BUILD_DOC)
    find_package(Doxygen REQUIRED)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(GURLSPP_BUILD_DOC)

