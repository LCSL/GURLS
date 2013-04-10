if( GurlsDependencies_FIND_COMPONENTS )
  foreach( component ${GurlsDependencies_FIND_COMPONENTS} )
    STRING( TOUPPER ${component} _COMPONENT )
    SET( GURLS_USE_${_COMPONENT} 1 )
  endforeach( component )
else()
    SET( GURLS_USE_BLAS_LAPACK 1 )
    SET( GURLS_USE_BOOST 1 )
endif()

set(GurlsDependencies_INCLUDE_DIRS "")
set(GurlsDependencies_LIBRARY_DIRS "")
set(GurlsDependencies_LIBRARIES "")
set(GurlsDependencies_DEFINITIONS "")

#################### BLAS/LAPACK

if(GURLS_USE_BLAS_LAPACK)

    if(CMAKE_COMPILER_IS_GNUCC OR MSVC)

        if(NOT MSVC)
            set(BLAS_LAPACK_IMPLEMENTATION "OPENBLAS" CACHE STRING "Specify your blas/lapack implementation (MKL, ACML, ATLAS, NETLIB, OPENBLAS or an empty string)")
        else()
            set(BLAS_LAPACK_IMPLEMENTATION "" CACHE STRING "Specify your blas/lapack implementation (MKL, ACML, ATLAS, NETLIB or an empty string)")
        endif()

        if(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

            unset(BLAS_LIBRARIES CACHE)
            unset(LAPACK_LIBRARIES CACHE)

            find_package(MKL)

            set(GurlsDependencies_INCLUDE_DIRS ${GurlsDependencies_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})
            set(GurlsDependencies_LIBRARY_DIRS ${GurlsDependencies_LIBRARY_DIRS} ${MKL_LIBRARY_DIRS})
            set(GurlsDependencies_DEFINITIONS ${GurlsDependencies_DEFINITIONS} ${MKL_DEFINITIONS})

            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${MKL_LIBRARIES})

        elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ACML")

            unset(BLAS_LIBRARIES CACHE)
            unset(LAPACK_LIBRARIES CACHE)

            find_package(ACML)
            set(GurlsDependencies_DEFINITIONS ${GurlsDependencies_DEFINITIONS} ${ACML_DEFINITIONS})
            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${ACML_LIBRARIES})


        elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "ATLAS")

            unset(BLAS_LIBRARIES CACHE)
            unset(LAPACK_LIBRARIES CACHE)

            find_package(ATLAS)
            set(GurlsDependencies_DEFINITIONS ${GurlsDependencies_DEFINITIONS} ${${ATLAS_DEFINITIONS})
            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${ATLAS_LAPACK_LIBS} ${ATLAS_BLAS_LIBS} ${ATLAS_LIBS} gfortran)

        elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "NETLIB")

            unset(BLAS_LIBRARIES CACHE)
            unset(LAPACK_LIBRARIES CACHE)

            enable_language(Fortran)
            find_package(BLAS)
            find_package(LAPACK)

            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

        elseif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "OPENBLAS")

            unset(BLAS_LIBRARIES CACHE)
            unset(LAPACK_LIBRARIES CACHE)

            set (Openblas_ROOT $ENV{GURLSPP_ROOT})
            enable_language(Fortran)
            find_package(Openblas)

            set(GurlsDependencies_INCLUDE_DIRS ${GurlsDependencies_INCLUDE_DIRS} ${Openblas_INCLUDE_DIRS})
            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${Openblas_LIBRARIES})


        else(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

            set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
            set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

            set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

        endif(BLAS_LAPACK_IMPLEMENTATION STREQUAL "MKL")

    else(CMAKE_COMPILER_IS_GNUCC OR MSVC)

        set(BLAS_LIBRARIES  "" CACHE FILEPATH "")
        set(LAPACK_LIBRARIES "" CACHE FILEPATH "")

        set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif(CMAKE_COMPILER_IS_GNUCC OR MSVC)

endif(GURLS_USE_BLAS_LAPACK)


#################### BOOST

if(GURLS_USE_BOOST)
    set(Boost_USE_MULTITHREADED      ON)
    set(Boost_USE_STATIC_RUNTIME    OFF)
    option(Boost_USE_STATIC_LIBS "Link statically against boost libs" ON)

    set(CMAKE_PREFIX_PATH $ENV{GURLSPP_ROOT} ${CMAKE_PREFIX_PATH})
    find_package( Boost 1.46.0 COMPONENTS serialization date_time filesystem unit_test_framework system signals REQUIRED)
    mark_as_advanced(Boost_DIR)

    if(Boost_FOUND)
        set(GurlsDependencies_INCLUDE_DIRS ${GurlsDependencies_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
        set(GurlsDependencies_LIBRARY_DIRS ${GurlsDependencies_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
        set(GurlsDependencies_LIBRARIES ${GurlsDependencies_LIBRARIES} ${Boost_LIBRARIES})
    endif(Boost_FOUND)

endif(GURLS_USE_BOOST)

SET(GurlsDependencies_FOUND 1)
