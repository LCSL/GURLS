cmake_minimum_required(VERSION 2.8)
project(gurls++)

set(GURLSLIBRARY gurls++)
set(SRCDIR src/)
set(HDRDIR include/)

set(gurls_headers
	${HDRDIR}common.h
	${HDRDIR}gurls.h
	${HDRDIR}gurls.hpp
	${HDRDIR}kernel.h
	${HDRDIR}norm.h
	${HDRDIR}exceptions.h
	${HDRDIR}utils.h
	${HDRDIR}basearray.h
	${HDRDIR}basearray.hpp
	${HDRDIR}gvec.h
	${HDRDIR}gvec.hpp
	${HDRDIR}gmat2d.h
	${HDRDIR}gmat2d.hpp
	${HDRDIR}gmath.h
	${HDRDIR}options.h
	${HDRDIR}optlist.h
	${HDRDIR}paramsel.h
	${HDRDIR}loocvprimal.h
	${HDRDIR}optimization.h
	${HDRDIR}rlsprimal.h
	${HDRDIR}pred.h)

set(gurls_sources
	${SRCDIR}kernel.cpp
	${SRCDIR}norm.cpp
	${SRCDIR}utils.cpp
	${SRCDIR}exceptions.cpp
	${SRCDIR}gvec.cpp
	${SRCDIR}gmat2d.cpp
	${SRCDIR}gmath.cpp
	${SRCDIR}options.cpp
	${SRCDIR}optlist.cpp
	${SRCDIR}paramsel.cpp
	${SRCDIR}loocvprimal.cpp
	${SRCDIR}optimization.cpp
	${SRCDIR}rlsprimal.cpp
	${SRCDIR}pred.cpp)

enable_language( Fortran )
INCLUDE( FindBLAS.cmake )
INCLUDE( FindLAPACK.cmake )
INCLUDE( FindBoost.cmake )

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost 1.40 REQUIRED serialization)

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories( ${Boost_LIBRARY_DIRS} )
endif(Boost_FOUND)

include_directories(${HDRDIR})

set(BUILD_SHARED_LIBS ON)

set(USE_BINARY_ARCHIVES OFF CACHE BOOL "If ON all the data structures in GURLS are saved/loaded using binary files.")

if(USE_BINARY_ARCHIVES)
	add_definitions(-DUSE_BINARY_ARCHIVES)
endif(USE_BINARY_ARCHIVES)

add_library(${GURLSLIBRARY} ${gurls_headers} ${gurls_sources})
target_link_libraries(${GURLSLIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} pthread gfortran)

add_executable(gurlsSimpleExample example.cpp)
target_link_libraries(gurlsSimpleExample ${GURLSLIBRARY} pthread)

add_executable(gurlsIO IOexample.cpp)
target_link_libraries(gurlsIO ${GURLSLIBRARY} ${Boost_LIBRARIES} pthread)

add_executable(example examplegurls.cpp)
target_link_libraries(example ${GURLSLIBRARY} ${Boost_LIBRARIES} pthread)

add_executable(gurlsHowToUseOPT howtouseopt.cpp)
target_link_libraries(gurlsHowToUseOPT ${GURLSLIBRARY} pthread)

add_executable(gurlsBenchmark benchmark.cpp)
target_link_libraries(gurlsBenchmark ${GURLSLIBRARY} blas lapack pthread)

# Set the following variable to ON if you want to enbale the use of YARP/ICUB
set(USE_YARP_ICUB OFF CACHE BOOL "Enables the compilation of the modules and examples that rely on YARP and ICUB libraries.")
# Set the following variable to ON if you want to enable the use of OpenCV
set(USE_OPENCV OFF CACHE BOOL "Enables the compilation of the modules and examples that rely on the OpenCV library")

# Set the following variable to ON if you want to enable the use of uBLAS
set(USE_UBLAS OFF CACHE BOOL "Enables the compilation of the modules and examples that rely on the uBLAS library")


IF(USE_YARP_ICUB)
	find_package(YARP)
	link_libraries(${YARP_LIBRARIES})
	include_directories(${YARP_INCLUDE_DIRS})
	set(CMAKE_MODULE_PATH ${YARP_MODULE_PATH} ${CMAKE_MODULE_PATH})

	find_package(ICUB)
	set(CMAKE_MODULE_PATH ${ICUB_MODULE_PATH} ${CMAKE_MODULE_PATH})
	include(iCubOptions)
	include(iCubHelpers)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ICUB_LINK_FLAGS}")

	add_executable(forChris forChris.cpp)
	target_link_libraries(forChris ${GURLSLIBRARY} ${YARP_LIBRARIES})
ENDIF(USE_YARP_ICUB)


IF(USE_OPENCV)
	find_package(OpenCV)
	include_directories(${OpenCV_INCLUDE_DIRS})

	add_executable(gurlsNopencv gurlsNopencv.cpp)
	target_link_libraries(gurlsNopencv  ${GURLSLIBRARY} ${OpenCV_LIBRARIES})
ENDIF(USE_OPENCV)


IF(USE_UBLAS)
	find_package(OpenCV)
	include_directories(${OpenCV_INCLUDE_DIRS})

	add_executable(ublas exampleublas.cpp)
	target_link_libraries(ublas  ${GURLSLIBRARY} ${Boost_LIBRARIES})
ENDIF(USE_UBLAS)
