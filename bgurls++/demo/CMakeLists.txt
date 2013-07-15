# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

include_directories(${BGurls++_INCLUDE_DIRS} ${Gurls++_INCLUDE_DIRS})
add_definitions( -D_BGURLS )
add_executable(bigdemo bigdemo.cpp)
target_link_libraries(bigdemo ${BGurls++_LIBRARIES})