cmake_minimum_required(VERSION 2.6)

project(gurls)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules ${CMAKE_MODULE_PATH})

function(add_all_executables SRCDIR LINK_LIBS)

    file(GLOB sources ${SRCDIR}*.cpp)

    set(executables "")

    foreach(source_file ${sources})

        get_filename_component(executable "${source_file}" NAME_WE)

        set(executables ${executables} ${executable})

        add_executable(${executable} ${source_file} )
        target_link_libraries(${executable} ${LINK_LIBS})

        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
        set_property(TARGET ${executable} PROPERTY INSTALL_RPATH ${INSTALL_RPATH} ${CMAKE_INSTALL_PREFIX}/lib)

    endforeach(source_file)


    INSTALL(TARGETS ${executables}
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
    )

endfunction(add_all_executables)

set(GURLS_BUILD_GURLSPP ON CACHE BOOL "Build gurls++")
set(GURLS_BUILD_BGURLSPP OFF CACHE BOOL "Build Bgurls++")

if(GURLS_BUILD_GURLSPP OR GURLS_BUILD_BGURLSPP)

    set(GURLS_USE_EXTERNALS ON CACHE BOOL "")

    if(GURLS_USE_EXTERNALS)

        #set(GURLS_USE_EXTERNAL_BLAS_LAPACK OFF CACHE BOOL "")
        set(GURLS_USE_EXTERNAL_BOOST ON CACHE BOOL "")

        if(GURLS_BUILD_BGURLSPP)
            set(GURLS_USE_EXTERNAL_HDF5 ON CACHE BOOL "")
        else()
            set(GURLS_USE_EXTERNAL_HDF5 OFF CACHE BOOL "")
        endif()

        set(GURLS_EXTERNALS_PREFIX ${CMAKE_INSTALL_PREFIX})

        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

    endif(GURLS_USE_EXTERNALS)

    #################### BLAS/LAPACK

    if(GURLS_USE_EXTERNAL_BLAS_LAPACK)

    else(GURLS_USE_EXTERNAL_BLAS_LAPACK)

        if(CMAKE_COMPILER_IS_GNUCC OR MSVC)

            set(BLAS_LAPACK_IMPLEMENTATION "" CACHE STRING "Specify your blas/lapack implementation (MKL, ACML, ATLAS, NETLIB or an empty string)")

            if(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

                unset(BLAS_LIBRARIES CACHE)
                unset(LAPACK_LIBRARIES CACHE)

                find_package(MKL)

                include_directories(${MKL_INCLUDE_DIRS})
                link_directories(${MKL_LIBRARY_DIRS})
                add_definitions( ${MKL_DEFINITIONS})

                set(BLAS_LIBRARIES "")
                set(LAPACK_LIBRARIES ${MKL_LIBRARIES})

            elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ACML")

                unset(BLAS_LIBRARIES CACHE)
                unset(LAPACK_LIBRARIES CACHE)

                find_package(ACML)
                set(BLAS_LIBRARIES "")
                add_definitions( ${ACML_DEFINITIONS})
                set(LAPACK_LIBRARIES ${ACML_LIBRARIES})


            elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ATLAS")

                unset(BLAS_LIBRARIES CACHE)
                unset(LAPACK_LIBRARIES CACHE)

                find_package(ATLAS)
                set(BLAS_LIBRARIES ${ATLAS_BLAS_LIBS} ${ATLAS_LIBS} gfortran)
                add_definitions( ${ATLAS_DEFINITIONS})
                set(LAPACK_LIBRARIES ${ATLAS_LAPACK_LIBS})

            elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "NETLIB")

                unset(BLAS_LIBRARIES CACHE)
                unset(LAPACK_LIBRARIES CACHE)

                enable_language(Fortran)
                find_package(BLAS)
                find_package(LAPACK)

            else(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

                set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
                set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

            endif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

        else(CMAKE_COMPILER_IS_GNUCC OR MSVC)

            set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
            set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

        endif(CMAKE_COMPILER_IS_GNUCC OR MSVC)

    endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)

    #################### BOOST

    if(GURLS_USE_EXTERNAL_BOOST)

        if(WIN32)
          set(prefix "")
          set(suffix ".lib")
        else()
          set(prefix "lib")
          set(suffix ".a")
        endif()

        set( Boost_LIBRARY_DIRS ${GURLS_EXTERNALS_PREFIX}/lib )
        set( Boost_INCLUDE_DIRS ${GURLS_EXTERNALS_PREFIX}/include )
        set( Boost_LIBRARIES
            ${prefix}boost_serialization-mt${suffix}
            ${prefix}boost_date_time-mt${suffix}
            ${prefix}boost_filesystem-mt${suffix}
            ${prefix}boost_unit_test_framework-mt${suffix}
            ${prefix}boost_system-mt${suffix}
            ${prefix}boost_signals-mt${suffix}
        )

        if(MSVC)
            add_definitions(-DBOOST_ALL_NO_LIB)
        endif()

        set( Boost_FOUND ON )

    else(GURLS_USE_EXTERNAL_BOOST)

        set(Boost_USE_MULTITHREADED      ON)
        set(Boost_USE_STATIC_RUNTIME    OFF)
        option(Boost_USE_STATIC_LIBS "Link statically against boost libs" OFF)
        find_package( Boost 1.46.0 COMPONENTS serialization date_time  filesystem unit_test_framework system signals REQUIRED)
        mark_as_advanced(Boost_DIR)

    endif(GURLS_USE_EXTERNAL_BOOST)

    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
        link_directories(${Boost_LIBRARY_DIRS})
    endif(Boost_FOUND)


    set(GURLS_USE_BINARY_ARCHIVES ON CACHE BOOL "If ON all the data structures in GURLS are saved/loaded using binary files.")

    if(GURLS_USE_BINARY_ARCHIVES)
        add_definitions(-DUSE_BINARY_ARCHIVES)
    endif(GURLS_USE_BINARY_ARCHIVES)

    if(GURLS_USE_EXTERNALS)
        unset (GURLS_BUILD_SHARED_LIBS CACHE )
        set (GURLS_BUILD_SHARED_LIBS OFF)
    else()
        set(GURLS_BUILD_SHARED_LIBS ON CACHE BOOL "Build a shared library (ON) or a static one (OFF)")
    endif()

    if(GURLS_BUILD_SHARED_LIBS)
        set( GURLS_LIB_LINK SHARED)
    else(GURLS_BUILD_SHARED_LIBS)
        set(GURLS_LIB_LINK STATIC)
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

endif(GURLS_BUILD_GURLSPP OR GURLS_BUILD_BGURLSPP)
