# Copyright (C) 2011-2013  Istituto Italiano di Tecnologia, Massachussets Institute of Techology
# Authors: Elena Ceseracciu <elena.ceseracciu@iit.it>, Matteo Santoro <msantoro@mit.edu>

add_executable(examplebgurls examplebgurls.cpp)
target_link_libraries(examplebgurls.cpp ${BGURLS_LIBRARIES})