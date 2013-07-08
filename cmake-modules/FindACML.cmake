# - Find ACML library
# This module finds an installed ACML library
#
# This module sets the following variables:
#  ACML_INCLUDE_DIRS - uncached list of include directories
#
#  ACML_DEFINITIONS - uncached list of required definitions
#
#  ACML_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use ACML

set(ACML_ROOT "$ENV{ACML_DIR}" CACHE PATH "Path to ACML")

set(ACML_LIB_PATH ${ACML_ROOT}/lib/)
set(ACML_STATIC OFF CACHE BOOL "Link statically against ACML libraries")

find_path(ACML_INCLUDE_DIRS acml.h ${ACML_ROOT}/include )

set(ACML_DEFINITIONS -D_ACML)

if(${ACML_ROOT} MATCHES "^.*mp.*$")
    set(ACML_MULTITHREADING "_mp")
else(${ACML_ROOT} MATCHES "^.*mp.*$")
    set(ACML_MULTITHREADING "")
endif(${ACML_ROOT} MATCHES "^.*mp.*$")

if(ACML_STATIC)
        set(SUFFIX "")
else(ACML_STATIC)
    if(MSVC)
        set(SUFFIX "_dll")
    else(MSVC)
        set(SUFFIX "")
    endif(MSVC)
endif(ACML_STATIC)
  
find_library(ACML_LIBRARIES libacml${ACML_MULTITHREADING}${SUFFIX} ${ACML_LIB_PATH})

if (ACML_INCLUDE_DIRS)
  set(ACML_FOUND TRUE)
endif (ACML_INCLUDE_DIRS)

if (ACML_FOUND)
  if (NOT ACML_FIND_QUIETLY)
	 message(STATUS "Found ACML: ${ACML_INCLUDE_DIRS}")
  endif (NOT ACML_FIND_QUIETLY)
else(ACML_FOUND)
  if (ACML_FIND_REQUIRED)
	 message(FATAL_ERROR "Could not find ACML")
  endif (ACML_FIND_REQUIRED)
endif (ACML_FOUND)