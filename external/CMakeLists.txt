
cmake_minimum_required(VERSION 2.6)

include(ExternalProject)


set(EXTERNAL_PREFIX ${GURLS_EXTERNALS_PREFIX})

if( ${CMAKE_VERSION} VERSION_GREATER "2.8.5" OR ${CMAKE_VERSION} VERSION_EQUAL "2.8.5")
    include(ProcessorCount)
    ProcessorCount(N)
    if(${N} GREATER 1)
        set(MAKE_ARGS -j${N})
    endif()
endif()


if(GURLS_USE_EXTERNAL_HDF5)

    ############## Zlib

    if(MSVC)
        ExternalProject_add(zlib
            URL http://zlib.net/zlib-1.2.7.tar.gz
            URL_MD5 60df6a37c56e7c1366cca812414f7b85
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/zlib
            #CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PREFIX}
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PREFIX}
        )
    else()
        ExternalProject_add(zlib
            URL http://zlib.net/zlib-1.2.7.tar.gz
            URL_MD5 60df6a37c56e7c1366cca812414f7b85
            BUILD_IN_SOURCE 1
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/zlib
            INSTALL_DIR ${EXTERNAL_PREFIX}/
            CONFIGURE_COMMAND <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR> --static
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
        )
    endif()

    set_target_properties(zlib PROPERTIES EXCLUDE_FROM_ALL 1)



    ############## MPICH2

    if(MSVC)
        set(DOWNLOAD_PATH ${EXTERNAL_PREFIX}/src/MPICH2)
        set(DOWNLOAD_URL http://www.mpich.org/static/tarballs/1.4.1p1)

        IF(CMAKE_CL_64)
            set(DOWNLOAD_FILE mpich2-1.4.1p1-win-x86-64.msi)
        else()
            set(DOWNLOAD_FILE mpich2-1.4.1p1-win-ia32.msi)
        endif()

        if (NOT EXISTS ${DOWNLOAD_PATH}/${DOWNLOAD_FILE})
            message(STATUS "Downloading MPICH2...")
            file(DOWNLOAD ${DOWNLOAD_URL}/${DOWNLOAD_FILE} ${DOWNLOAD_PATH}/${DOWNLOAD_FILE} SHOW_PROGRESS)
            message( STATUS "Done")
        endif()

        string(REPLACE "/" "\\" TARGET_PREFIX ${EXTERNAL_PREFIX})
        string(REPLACE "/" "\\" DOWNLOAD_PATH ${DOWNLOAD_PATH})

        ExternalProject_add(MPICH2
            DOWNLOAD_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            #INSTALL_COMMAND msiexec /passive TARGETDIR=${TARGET_PREFIX} /i "${DOWNLOAD_PATH}\\${DOWNLOAD_FILE}"
            INSTALL_COMMAND msiexec TARGETDIR=${TARGET_PREFIX} /i "${DOWNLOAD_PATH}\\${DOWNLOAD_FILE}"
        )
    else()
        ExternalProject_add(MPICH2
            URL http://www.mpich.org/static/tarballs/3.0.2/mpich-3.0.2.tar.gz
            BUILD_IN_SOURCE 1
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/MPICH2
            INSTALL_DIR ${EXTERNAL_PREFIX}/
            CONFIGURE_COMMAND CFLAGS=-w CXXFLAGS=-w <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR> --disable-f77 --disable-fc --enable-g=none --enable-romio --enable-shared --disable-static
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
        )
    endif()

    set_target_properties(MPICH2 PROPERTIES EXCLUDE_FROM_ALL 1)



    ############## hdf5
    if(MSVC)

        ExternalProject_add(hdf5
            URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/hdf5
            CMAKE_ARGS
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_INSTALL_PREFIX=${EXTERNAL_PREFIX}
                -DCMAKE_PREFIX_PATH=${EXTERNAL_PREFIX}
                -DBUILD_SHARED_LIBS=OFF
                -DBUILD_STATIC_EXECS=OFF
                -DBUILD_TESTING=OFF
                -DHDF5_BUILD_EXAMPLES=OFF
                -DHDF5_BUILD_FORTRAN=OFF
                -DHDF5_BUILD_HL_LIB=OFF
                -DHDF5_BUILD_TOOLS=OFF
                -DHDF5_DISABLE_COMPILER_WARNINGS=ON
                -DHDF5_ENABLE_PARALLEL=ON
                -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
                -DZLIB_INCLUDE_DIR=${EXTERNAL_PREFIX}/include
                -DZLIB_LIBRARY=${EXTERNAL_PREFIX}/lib/zlibstatic.lib
                -DZLIB_USE_EXTERNAL=OFF
        )
    else()
        ExternalProject_add(hdf5
            URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
            BUILD_IN_SOURCE 1
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/hdf5
            INSTALL_DIR ${EXTERNAL_PREFIX}/
            CONFIGURE_COMMAND CC=${EXTERNAL_PREFIX}/bin/mpicc CFLAGS=-w <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR> --enable-production --enable-parallel --disable-shared --disable-fortran --disable-hl --with-zlib=${EXTERNAL_PREFIX}/
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
        )
    endif()
    add_dependencies(hdf5 zlib MPICH2)

    set_target_properties(hdf5 PROPERTIES EXCLUDE_FROM_ALL 1)

endif(GURLS_USE_EXTERNAL_HDF5)

############## Boost

if(GURLS_USE_EXTERNAL_BOOST)

    if(MSVC)
        set(EXT "bat")
    else()
        set(EXT "sh")
    endif()

    #address-model=64

    ExternalProject_add(boost
        URL http://downloads.sourceforge.net/boost/boost_1_53_0.tar.gz
        URL_MD5 57a9e2047c0f511c4dfcf00eb5eb2fbb
        BUILD_IN_SOURCE 1
        SOURCE_DIR  ${EXTERNAL_PREFIX}/src/boost
        INSTALL_DIR ${EXTERNAL_PREFIX}/
        CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.${EXT}
        BUILD_COMMAND <SOURCE_DIR>/b2 ${MAKE_ARGS} -d0 --layout=tagged variant=release link=static threading=multi runtime-link=shared
                    --prefix=${EXTERNAL_PREFIX}
                    --with-serialization
                    --with-date_time
                    --with-filesystem
                    --with-test
                    --with-system
                    --with-signals
                    install
        INSTALL_COMMAND ""
    )

    set_target_properties(boost PROPERTIES EXCLUDE_FROM_ALL 1)

endif(GURLS_USE_EXTERNAL_BOOST)


if(GURLS_USE_EXTERNAL_BLAS_LAPACK)

    ExternalProject_add(openblas
        GIT_REPOSITORY "git://github.com/xianyi/OpenBLAS.git"
        GIT_TAG "v0.2.6"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        SOURCE_DIR  ${EXTERNAL_PREFIX}/src/openblas
        INSTALL_DIR ${EXTERNAL_PREFIX}/
        CONFIGURE_COMMAND ""
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.install.openblas" "<SOURCE_DIR>/Makefile.install"
        BUILD_COMMAND NO_SHARED=1 ${CMAKE_MAKE_PROGRAM} libs netlib
        INSTALL_COMMAND NO_SHARED=1 ${CMAKE_MAKE_PROGRAM} install PREFIX=<INSTALL_DIR>
    )

    set_target_properties(openblas PROPERTIES EXCLUDE_FROM_ALL 1)

endif(GURLS_USE_EXTERNAL_BLAS_LAPACK)

#add_custom_target(External
#  DEPENDS netCDF hdf5 zlib MPICH2
#  COMMENT "Building external targets"
#)
