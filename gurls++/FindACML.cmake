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

set(ACML_ROOT "" CACHE PATH "Path to ACML")

set(ACML_LIB_PATH ${ACML_ROOT}/lib/)
set(ACML_STATIC OFF CACHE BOOL "Link statically against ACML libraries")

set(ACML_INCLUDE_DIRS ${ACML_ROOT}/include)

set(ACML_DEFINITIONS -D_ACML)

if(${ACML_ROOT} MATCHES "^.*mp.*$")
    set(ACML_MULTITHREADING "_mp")
else(${ACML_ROOT} MATCHES "^.*mp.*$")
    set(ACML_MULTITHREADING "")
endif(${ACML_ROOT} MATCHES "^.*mp.*$")

if(ACML_STATIC)
    if(MSVC)
        set(PREFIX ".lib")
    else(MSVC)
        set(PREFIX ".a")
    endif(MSVC)
else(ACML_STATIC)
    if(MSVC)
        set(PREFIX "_dll.lib")
    else(MSVC)
        set(PREFIX ".so")
    endif(MSVC)
endif(ACML_STATIC)

set(ACML_LIBRARIES ${ACML_LIB_PATH}libacml${ACML_MULTITHREADING}${PREFIX})
