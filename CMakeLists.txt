# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

cmake_minimum_required(VERSION 2.8.1)

project(gurls)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules ${CMAKE_MODULE_PATH})

#function(add_all_executables SRCDIR LINK_LIBS)
#
#    file(GLOB sources ${SRCDIR}*.cpp)
#
#    set(executables "")
#
#    foreach(source_file ${sources})
#
#        get_filename_component(executable "${source_file}" NAME_WE)
#
#        set(executables ${executables} ${executable})
#
#        add_executable(${executable} ${source_file} )
#        target_link_libraries(${executable} ${LINK_LIBS})
#
#        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH ${INSTALL_RPATH} ${CMAKE_INSTALL_PREFIX}/lib)
#
#    endforeach(source_file)
#
#
#    INSTALL(TARGETS ${executables}
#      RUNTIME DESTINATION bin
#      LIBRARY DESTINATION lib
#      ARCHIVE DESTINATION lib
#    )
#
#endfunction(add_all_executables)

set(EXTERNAL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external)
include(GurlsFindDependencies)
if (BLAS_LAPACK_FOUND AND Boost_FOUND) #boost should be optional at some point though
    option(GURLS_BUILD_GURLSPP "Build gurls++" ON)
else()
    option(GURLS_BUILD_GURLSPP "Build gurls++" OFF)
endif()

if(GURLS_BUILD_GURLSPP)

    if(HDF5_FOUND)
        option(GURLS_BUILD_BGURLSPP "Build Bgurls++" OFF)
    else()
        option(GURLS_BUILD_BGURLSPP "Build Bgurls++" OFF)
    endif()
#   option(GURLS_USE_EXTERNALS "" ON)
#    option(GURLS_USE_BOOST "" ON)

    # Find available libraries in the system

    # Allow downloading and compilation of external libraries, if needed/wanted
    #    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

#     if(GURLS_USE_EXTERNALS)
#
#         if(MSVC)
#             set(GURLS_USE_EXTERNAL_BLAS_LAPACK OFF)
#         else()
#             set(GURLS_USE_EXTERNAL_BLAS_LAPACK ON CACHE BOOL "")
#         endif()
# 
#         set(GURLS_USE_EXTERNAL_BOOST ON CACHE BOOL "")
# 
#         if(GURLS_BUILD_BGURLSPP)
#             set(GURLS_USE_EXTERNAL_HDF5 ON CACHE BOOL "")
#         else()
#             set(GURLS_USE_EXTERNAL_HDF5 OFF CACHE BOOL "")
#         endif()
# 
#         set(GURLS_EXTERNALS_PREFIX ${CMAKE_INSTALL_PREFIX})
# 
#         add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)
# 
#     else(GURLS_USE_EXTERNALS)
# 
#         unset(GURLS_USE_EXTERNAL_BLAS_LAPACK CACHE)
#         unset(GURLS_USE_EXTERNAL_BOOST CACHE)
#         unset(GURLS_USE_EXTERNAL_HDF5 CACHE)
# 
#     endif(GURLS_USE_EXTERNALS)
# 
#     #################### BLAS/LAPACK
# 
#     if(GURLS_USE_EXTERNAL_BLAS_LAPACK)
# 
#         enable_language(Fortran)
#         link_directories(${GURLS_EXTERNALS_PREFIX}/lib)
#         set(GurlsExtDependencies_LIBRARIES libopenblas.a gfortran pthread)
# 
#     else(GURLS_USE_EXTERNAL_BLAS_LAPACK)
# 
#             set(DEPENDENCIES_COMPONENTS ${DEPENDENCIES_COMPONENTS} BLAS_LAPACK)
# 
#     endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)
# 
#     #################### BOOST
# 
#     if(GURLS_USE_EXTERNAL_BOOST)
# 
#         if(WIN32)
#           set(prefix "lib")
#           set(suffix ".lib")
#         else()
#           set(prefix "lib")
#           set(suffix ".a")
#         endif()
# 
#         set( Boost_LIBRARY_DIRS ${GURLS_EXTERNALS_PREFIX}/lib )
#         set( Boost_INCLUDE_DIRS ${GURLS_EXTERNALS_PREFIX}/include )
#         set( Boost_LIBRARIES
#             ${prefix}boost_serialization-mt${suffix}
#             ${prefix}boost_date_time-mt${suffix}
#             ${prefix}boost_filesystem-mt${suffix}
#             ${prefix}boost_unit_test_framework-mt${suffix}
#             ${prefix}boost_system-mt${suffix}
#             ${prefix}boost_signals-mt${suffix}
#         )
#
#        if(MSVC)
#            add_definitions(-DBOOST_ALL_NO_LIB)
#        endif()
#
#        include_directories(${Boost_INCLUDE_DIRS})
#        link_directories(${Boost_LIBRARY_DIRS})
#        set(GurlsExtDependencies_LIBRARIES ${GurlsExtDependencies_LIBRARIES} ${Boost_LIBRARIES})
#
#    else(GURLS_USE_EXTERNAL_BOOST)
#
#       set(DEPENDENCIES_COMPONENTS ${DEPENDENCIES_COMPONENTS} BOOST)
#
#    endif(GURLS_USE_EXTERNAL_BOOST)
#
#    if(DEPENDENCIES_COMPONENTS)
#        find_package(GurlsDependencies COMPONENTS ${DEPENDENCIES_COMPONENTS} REQUIRED)
#
#        include_directories(${GurlsDependencies_INCLUDE_DIRS})
#        link_directories(${GurlsDependencies_LIBRARY_DIRS})
#        add_definitions(${GurlsDependencies_DEFINITIONS})
#    endif()
#
#    set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${GurlsExtDependencies_LIBRARIES})

    option(GURLS_USE_BINARY_ARCHIVES "If ON all the data structures in GURLS are saved/loaded using binary files." ON)

    if(GURLS_USE_BINARY_ARCHIVES)
        add_definitions(-DUSE_BINARY_ARCHIVES)
    endif(GURLS_USE_BINARY_ARCHIVES)

    if(GURLS_USE_EXTERNALS)
        unset (GURLS_BUILD_SHARED_LIBS CACHE )
        set (GURLS_BUILD_SHARED_LIBS OFF) #why?
    else()
        set(GURLS_BUILD_SHARED_LIBS OFF CACHE BOOL "Build a shared library (ON) or a static one (OFF)")
    endif()

    if(GURLS_BUILD_SHARED_LIBS)
        set(GURLS_LIB_LINK SHARED)
    else(GURLS_BUILD_SHARED_LIBS)
        set(GURLS_LIB_LINK STATIC)
        add_definitions( -D_GURLS_STATIC)
    endif(GURLS_BUILD_SHARED_LIBS)

    #add_definitions( -Wall )
    if(MSVC)
        add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    endif(MSVC)


    if(GURLS_BUILD_GURLSPP)
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gurls++)
    endif(GURLS_BUILD_GURLSPP)

    if(GURLS_BUILD_BGURLSPP)
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bgurls++)
    endif(GURLS_BUILD_BGURLSPP)

#TODO: create GURLSconfig.cmake etc.

else(GURLS_BUILD_GURLSPP)
	message(WARNING "You have disabled compilation!")
endif(GURLS_BUILD_GURLSPP)
