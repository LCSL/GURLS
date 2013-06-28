# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>
include_directories(${GURLS_INCLUDE_DIRS})

add_executable(GURLSloocvprimal GURLSloocvprimal.cpp)
target_link_libraries(GURLSloocvprimal ${GURLS_LIBRARIES})

add_executable(helloworld helloworld.cpp)
target_link_libraries(helloworld ${GURLS_LIBRARIES})

add_executable(recursiveRLS recursiveRLS.cpp)
target_link_libraries(recursiveRLS ${GURLS_LIBRARIES})