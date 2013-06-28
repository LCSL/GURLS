# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
    message( FATAL_ERROR "You are trying to run CMake from the gurls++ directory, instead of just from the top directory")
endif()

set(GURLSLIBRARY gurls++)
project(${GURLSLIBRARY})

#set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/)

#set(DEMODIR ${CMAKE_CURRENT_SOURCE_DIR}/demo/)
#set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR}/test/)
#set(MISCDIR ${CMAKE_CURRENT_SOURCE_DIR}/misc/)


set(gurls_headers   include/basearray.h
                    include/basearray.hpp
                    include/blas_lapack.h
                    include/boltzmangap.h
                    include/boltzman.h
                    include/calibratesgd.h
                    include/chisquaredkernel.h
                    include/common.h
                    include/confidence.h
                    include/dual.h
                    include/exceptions.h
                    include/exports.h
                    include/fixlambda.h
                    include/fixsiglam.h
                    include/gap.h
                    include/gmat2d.h
                    include/gmat2d.hpp
                    include/gmath.h
                    include/gurls.h
                    include/gurls.hpp
                    include/gvec.h
                    include/gvec.hpp
                    include/hodual.h
                    include/hogpregr.h
                    include/hoprimal.h
                    include/icholwrapper.h
                    include/icholwrapper.hpp
                    include/kernel.h
                    include/linearkernel.h
                    include/loocvdual.h
                    include/loocvprimal.h
                    include/loogpregr.h
                    include/macroavg.h
                    include/maxscore.h
                    include/norm.h
                    include/norml2.h
                    include/normtestzscore.h
                    include/normzscore.h
                    include/nystromwrapper.h
                    include/nystromwrapper.hpp
                    include/optarray.h
                    include/optfunction.h
                    include/optimization.h
                    include/options.h
                    include/optlist.h
                    include/optmatrix.h
                    include/paramsel.h
                    include/perf.h
                    include/precisionrecall.h
                    include/predgp.h
                    include/pred.h
                    include/predkernel.h
                    include/predkerneltraintest.h
                    include/primal.h
                    include/rbfkernel.h
                    include/recrlswrapper.h
                    include/recrlswrapper.hpp
                    include/rlsauto.h
                    include/rlsdual.h
                    include/rlsdualr.h
                    include/rlsgp.h
                    include/rlspegasos.h
                    include/rlsprimal.h
                    include/rlsprimalrecinit.h
                    include/rlsprimalrecupdate.h
                    include/rlsprimalr.h
                    include/rmse.h
                    include/serialization.h
                    include/siglam.h
                    include/siglamhogpregr.h
                    include/siglamho.h
                    include/siglamloogpregr.h
                    include/split.h
                    include/splitho.h
                    include/utils.h
                    include/wrapper.h
                    include/wrapper.hpp
    )

set(gurls_sources   src/gmath.cpp
                    src/optarray.cpp
                    src/optfunction.cpp
                    src/options.cpp
                    src/optlist.cpp
                    src/optmatrix.cpp
    )
set(GURLS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include/ CACHE INTERNAL "")
#if (GURLS_USE_BOOST) #should change something in the code as well...
    list(APPEND GURLS_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
    set (GurlsDependencies_LIBRARIES ${BLAS_LAPACK_LIBRARIES} ${Boost_LIBRARIES}) #should use only the serialization one?
#else(GURLS_USE_BOOST)
#    set (GurlsDependencies_LIBRARIES ${BLAS_LAPACK_LIBRARIES})
#endif(GURLS_USE_BOOST)

# if(GURLS_BUILD_BGURLSPP)
#     add_definitions( -D_BGURLS )
# endif()

include_directories( BEFORE ${GURLS_INCLUDE_DIRS} )
add_library(${GURLSLIBRARY} ${GURLS_LIB_LINK} ${gurls_headers} ${gurls_sources} )

target_link_libraries(${GURLSLIBRARY} ${GurlsDependencies_LIBRARIES})

if(GURLS_USE_EXTERNAL_BLAS_LAPACK)
    add_dependencies(${GURLSLIBRARY} buildOpenblas)
endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)

if(GURLS_USE_EXTERNAL_BOOST)
    add_dependencies(${GURLSLIBRARY} buildBoost)
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

set(GURLS_LIBRARIES
    ${GURLSLIBRARY}
    ${GurlsDependencies_LIBRARIES} CACHE INTERNAL ""
)

option(GURLSPP_BUILD_DEMO "" ON)
if(GURLSPP_BUILD_DEMO)
    add_subdirectory(demo)
endif(GURLSPP_BUILD_DEMO)

option(GURLSPP_BUILD_TEST "" OFF)
mark_as_advanced(FORCE GURLSPP_BUILD_TEST)
if(GURLSPP_BUILD_TEST)
#     set(GURLSPP_DATA_DIR "" CACHE PATH "Path to the Gurls++ data directory")
#     mark_as_advanced(FORCE GURLSPP_DATA_DIR)
# 
#     if(GURLSPP_DATA_DIR STREQUAL "")
#         message(WARNING, "GURLSPP_DATA_DIR value is set to an empty string. Please be sure this is the right configuration for you")
#     endif(GURLSPP_DATA_DIR STREQUAL "")
# 
#     add_definitions(-DGURLS_DATA_DIR="${GURLSPP_DATA_DIR}")
   add_subdirectory(test)
#   add_all_executables(${TESTDIR} ${GURLS_LINK_LIBRARIES})
endif(GURLSPP_BUILD_TEST)

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

    if (Doxygen_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(gurlsppdoc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
    else(Doxygen_FOUND)
        message(WARNING "Doxygen documentation was enabled, but the Doxygen package was not found.")
    endif(Doxygen_FOUND)
endif(GURLSPP_BUILD_DOC)

