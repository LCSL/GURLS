# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

if( ${CMAKE_VERSION} VERSION_GREATER "2.8.5" OR ${CMAKE_VERSION} VERSION_EQUAL "2.8.5")
    include(ProcessorCount)
    ProcessorCount(N)
    if(${N} GREATER 1)
        set(MAKE_ARGS -j${N})
    endif()
endif()

include(ExternalProject)

 ############## Zlib

    set(hdf5_dependencies)

    find_package(ZLIB)
    if( NOT ZLIB_FOUND)
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
        set(ZLIB_FOUND TRUE)
        set(ZLIB_INCLUDE_DIRS ${EXTERNAL_PREFIX}/include) #cached?
        #set(ZLIB_LIBRARY_DIRS ${EXTERNALS_PREFIX}/lib) #cached?
        if(MSVC)
            set(ZLIB_LIBRARIES ${EXTERNAL_PREFIX}/lib/zlibstatic.lib)
        else()
            set(ZLIB_LIBRARIES ${EXTERNAL_PREFIX}/lib/libz.a)
        endif()
        list(APPEND hdf5_dependencies zlib)

    endif()


    ############## MPICH2

    find_package(MPI)
    if(NOT MPI_FOUND)
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
        set(MPI_CXX_FOUND TRUE)
        set(MPI_CXX_LIBRARIES mpichcxx mpich opa mpl rt pthread)
        set(MPI_CXX_COMPILE_FLAGS "")
        set(MPI_CXX_LINK_FLAGS "-Wl,-rpath  -Wl,${EXTERNAL_PREFIX}/lib")
        set(MPI_CXX_INCLUDE_PATH ${EXTERNAL_PREFIX}/include)
        set(MPI_CXX_COMPILER ${EXTERNAL_PREFIX}/bin/mpicxx)
        set(MPI_C_COMPILER ${EXTERNAL_PREFIX}/bin/mpicc)

        list(APPEND hdf5_dependencies MPICH2)

    endif()

    ############## hdf5
    if(MSVC)

        ExternalProject_add(buildHdf5
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
                 -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIR}
                 -DZLIB_LIBRARY=${ZLIB_LIBRARIES}
                 -DZLIB_USE_EXTERNAL=OFF
        )
        set(HDF5_LIBRARIES optimized ${EXTERNAL_PREFIX}/lib/hdf5.lib debug ${EXTERNAL_PREFIX}/lib/hdf5d.lib)
    else()
        ExternalProject_add(buildHdf5
            URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
            BUILD_IN_SOURCE 1
            SOURCE_DIR  ${EXTERNAL_PREFIX}/src/hdf5
            INSTALL_DIR ${EXTERNAL_PREFIX}/
            CONFIGURE_COMMAND CC=${MPI_C_COMPILER} CFLAGS=-w <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR> --enable-production --enable-parallel --disable-shared --disable-fortran --disable-hl --with-zlib=${EXTERNAL_PREFIX}/
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} ${MAKE_ARGS}
        )
            set(HDF5_LIBRARIES ${EXTERNAL_PREFIX}/lib/libhdf5.a)
    endif()
    if(hdf5_dependencies)
        add_dependencies(buildHdf5 ${hdf5_dependencies})
    endif()

    set_target_properties(buildHdf5 PROPERTIES EXCLUDE_FROM_ALL 1)
	set(HDF5_FOUND TRUE)
    set(HDF5_INCLUDE_DIRS ${EXTERNAL_PREFIX}/include CACHE INTERNAL "")
    set(HDF5_LIBRARY_DIRS ${EXTERNAL_PREFIX}/lib CACHE INTERNAL "")
    set(HDF5_DEFINITIONS -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE)
