############################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
SET( BASIC_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/libbasic )
SET( IO_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libio )
SET( BLOCKMATCHING_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/blockmatching ) 

# Independants
add_subdirectory(libio)
add_subdirectory(libbasic)

# Depend on BASIC and IO
add_subdirectory(blockmatching)








