# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

include_directories(${BGURLS_INCLUDE_DIRS})
add_executable(bigdemo bigdemo.cpp)
target_link_libraries(bigdemo ${BGURLS_LIBRARIES})