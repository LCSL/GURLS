# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

include_directories(${Gurls++_INCLUDE_DIRS})

add_executable(examplegurls examplegurls.cpp)
target_link_libraries(examplegurls ${Gurls++_LIBRARIES})
