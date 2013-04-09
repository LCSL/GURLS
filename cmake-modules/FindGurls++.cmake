# - Finds Gurls++ and all dependencies
# Once done this will define
#  Gurls++_FOUND - System has Gurls++
#  Gurls++_INCLUDE_DIRS - The Gurls++ include directories
#  Gurls++_LIBRARY_DIRS - The Gurls++ library directories
#  Gurls++_LIBRARIES - The libraries needed to use Gurls++
#  Gurls++_DEFINITIONS - Compiler switches required for using Gurls++

find_package(GurlsDependencies)


set(Gurls++_USE_BINARY_ARCHIVES ON CACHE BOOL "If ON all the data structures in GURLS are saved/loaded using binary files.")

find_library(Gurls++_LIBRARY
        NAMES libgurls++.a libgurls+.so libgurls++.lib gurls++.lib
        PATHS ${Gurls++_ROOT}/lib
)

find_path(Gurls++_INCLUDE_DIR
        NAMES gurls.h
        PATHS ${Gurls++_ROOT}/include/gurls++
)

if( ( Gurls++_LIBRARY STREQUAL "") OR ( Gurls++_INCLUDE_DIR STREQUAL "") )
    set(Gurls++_ROOT "" CACHE PATH "Path to the root of a Gurls++ installation")
    set(Gurls++_FOUND 0)
    message(WARNING "Gurls++ not found. Please try specifying Gurls++_ROOT")
else()
    set(Gurls++_FOUND 1)
    set(Gurls++_INCLUDE_DIRS ${Gurls++_INCLUDE_DIR} ${GurlsDependencies_INCLUDE_DIRS})
    set(Gurls++_LIBRARY_DIRS ${GurlsDependencies_LIBRARY_DIRS})
    set(Gurls++_LIBRARIES ${Gurls++_LIBRARY} ${GurlsDependencies_LIBRARIES})
    set(Gurls++_DEFINITIONS ${GurlsDependencies_DEFINITIONS})

    if(Gurls++_USE_BINARY_ARCHIVES)
        set(Gurls++_DEFINITIONS ${Gurls++_DEFINITIONS} -DUSE_BINARY_ARCHIVES)
    endif(Gurls++_USE_BINARY_ARCHIVES)

endif()

