# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

cmake_minimum_required(VERSION 2.8.1)

project(gurls)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")


set(EXTERNAL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external)
set(BLAS_LAPACK_IMPLEMENTATION "OPENBLAS" CACHE STRING "Possible values: ACML, ATLAS, MKL, NETLIB, OPENBLAS")
if(WIN32)
	if (NOT CMAKE_CL64)
		set (Openblas_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/win32")
	else()
		message(WARNING "Openblas library for 64-bit compilers is not shipped yet")
	endif()
endif()
set(BOOST_MINIMUM_VERSION "1.46.0")
include(GurlsFindDependencies)
option(GURLS_BUILD_GURLSPP "Build gurls++" ON)

if (NOT BLAS_LAPACK_FOUND OR NOT Boost_FOUND) #boost should be optional at some point though
    set(GURLS_BUILD_GURLSPP OFF)
endif()

if(GURLS_BUILD_GURLSPP)

    if(HDF5_FOUND)
        option(GURLS_BUILD_BGURLSPP "Build Bgurls++" OFF)
    else()
        option(GURLS_BUILD_BGURLSPP "Build Bgurls++" OFF)
    endif()

    option(GURLS_USE_BINARY_ARCHIVES "If ON all the data structures in GURLS are saved/loaded using binary files." ON)

    if(GURLS_USE_BINARY_ARCHIVES)
        add_definitions(-DUSE_BINARY_ARCHIVES)
    endif(GURLS_USE_BINARY_ARCHIVES)

    if(GURLS_USE_EXTERNAL_BLAS_LAPACK OR GURLS_USE_EXTERNAL_BOOST OR GURLS_USE_EXTERNAL_HDF5)
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

	#Actually build something
    if(GURLS_BUILD_GURLSPP)
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gurls++)
    endif(GURLS_BUILD_GURLSPP)

    if(GURLS_BUILD_BGURLSPP)
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bgurls++)
    endif(GURLS_BUILD_BGURLSPP)

	### WRITE CMAKE CONFIG FILES
    #Gurls++Config.cmake for build tree
	set (export_include_dirs ${Gurls++_INCLUDE_DIRS})
	set (generated_export_file "Gurls++targets-buildtree.cmake")
	set (export_libraries)
	set (export_definitions)
	
	if (GURLS_BUILD_GURLSPP)
		set (export_libraries ${export_libraries} gurls++)
		get_property(gurlspp_export_definitions DIRECTORY gurls++ PROPERTY COMPILE_DEFINITIONS) 
		set (export_definitions ${export_definitions} ${gurlspp_export_definitions})
	endif()
	if (GURLS_BUILD_BGURLSPP)
		set (export_libraries ${export_libraries} bgurls++)
		get_property(bgurlspp_export_definitions DIRECTORY bgurls++ PROPERTY COMPILE_DEFINITIONS) 
		set (export_definitions ${export_definitions} ${bgurlspp_export_definitions}) #TODO there may be duplicates...
	endif()
	configure_file(cmake-modules/Gurls++Config.cmake.in ${PROJECT_BINARY_DIR}/Gurls++Config.cmake @ONLY)
	
	#Gurls++Config.cmake for installation
	set (export_include_dirs "${CMAKE_INSTALL_PREFIX}/include")
	set (generated_export_file "Gurls++targets-install.cmake")
	configure_file(cmake-modules/Gurls++Config.cmake.in ${PROJECT_BINARY_DIR}/Gurls++ConfigForInstall.cmake @ONLY)
	install(FILES ${PROJECT_BINARY_DIR}/Gurls++ConfigForInstall.cmake DESTINATION lib RENAME Gurls++Config.cmake)
	
	# uninstall target
	configure_file(
    "cmake-modules/gurls-uninstall.cmake.in"
    "${PROJECT_BINARY_DIR}/gurls-uninstall.cmake"
    IMMEDIATE @ONLY)

	add_custom_target(uninstall
		COMMAND ${CMAKE_COMMAND} -P "${PROJECT_BINARY_DIR}/gurls-uninstall.cmake")

else(GURLS_BUILD_GURLSPP)
    message(WARNING "You have disabled compilation!")
endif(GURLS_BUILD_GURLSPP)
