# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
    message( FATAL_ERROR "You are trying to run CMake from the gurls++ directory, instead of just from the top directory")
endif()

set(GURLSLIBRARY gurls++)
project(${GURLSLIBRARY})

set(gurls_headers   include/gurls++/basearray.h
                    include/gurls++/basearray.hpp
                    include/gurls++/blas_lapack.h
                    include/gurls++/boltzmangap.h
                    include/gurls++/boltzman.h
                    include/gurls++/calibratesgd.h
                    include/gurls++/chisquaredkernel.h
                    include/gurls++/common.h
                    include/gurls++/confidence.h
                    include/gurls++/dual.h
                    include/gurls++/exceptions.h
                    include/gurls++/exports.h
                    include/gurls++/fixlambda.h
                    include/gurls++/fixsiglam.h
                    include/gurls++/gap.h
                    include/gurls++/gmat2d.h
                    include/gurls++/gmat2d.hpp
                    include/gurls++/gmath.h
                    include/gurls++/gurls.h
                    include/gurls++/gurls.hpp
                    include/gurls++/gvec.h
                    include/gurls++/gvec.hpp
                    include/gurls++/hodual.h
                    include/gurls++/hogpregr.h
                    include/gurls++/hoprimal.h
                    include/gurls++/icholwrapper.h
                    include/gurls++/icholwrapper.hpp
                    include/gurls++/kernel.h
                    include/gurls++/kernelrlswrapper.h
                    include/gurls++/kernelrlswrapper.hpp
                    include/gurls++/linearkernel.h
                    include/gurls++/loocvdual.h
                    include/gurls++/loocvprimal.h
                    include/gurls++/loogpregr.h
                    include/gurls++/macroavg.h
                    include/gurls++/maxscore.h
                    include/gurls++/norm.h
                    include/gurls++/norml2.h
                    include/gurls++/normtestzscore.h
                    include/gurls++/normzscore.h
                    include/gurls++/nystromwrapper.h
                    include/gurls++/nystromwrapper.hpp
                    include/gurls++/optarray.h
                    include/gurls++/optfunction.h
                    include/gurls++/optimization.h
                    include/gurls++/options.h
                    include/gurls++/optlist.h
                    include/gurls++/optmatrix.h
                    include/gurls++/paramsel.h
                    include/gurls++/perf.h
                    include/gurls++/precisionrecall.h
                    include/gurls++/predgp.h
                    include/gurls++/pred.h
                    include/gurls++/predkernel.h
                    include/gurls++/predkerneltraintest.h
                    include/gurls++/primal.h
                    include/gurls++/rbfkernel.h
                    include/gurls++/recrlswrapper.h
                    include/gurls++/recrlswrapper.hpp
                    include/gurls++/rlsauto.h
                    include/gurls++/rlsdual.h
                    include/gurls++/rlsdualr.h
                    include/gurls++/rlsgp.h
                    include/gurls++/rlspegasos.h
                    include/gurls++/rlsprimal.h
                    include/gurls++/rlsprimalrecinit.h
                    include/gurls++/rlsprimalrecupdate.h
                    include/gurls++/rlsprimalr.h
                    include/gurls++/rlswrapper.h
                    include/gurls++/rlswrapper.hpp
                    include/gurls++/rmse.h
                    include/gurls++/serialization.h
                    include/gurls++/siglam.h
                    include/gurls++/siglamhogpregr.h
                    include/gurls++/siglamho.h
                    include/gurls++/siglamloogpregr.h
                    include/gurls++/split.h
                    include/gurls++/splitho.h
                    include/gurls++/utils.h
                    include/gurls++/wrapper.h
                    include/gurls++/wrapper.hpp
    )

set(gurls_sources   src/gmath.cpp
                    src/optarray.cpp
                    src/optfunction.cpp
                    src/options.cpp
                    src/optlist.cpp
                    src/optmatrix.cpp
    )

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
install(EXPORT Gurlstargets-install DESTINATION lib)

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

