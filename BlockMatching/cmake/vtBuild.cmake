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

## #################################################################
## 
## #################################################################
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

## #################################################################
## Configure path
## #################################################################
set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})


## #################################################################
## Build option
## #################################################################
# Compilation Option

# to add flags for the pre-compiler
# SET( CPPFLAGS "" )
# add_definitions( ${CPPFLAGS} )

# SET(CFLAGS "-ansi -Wall -fsigned-char -fsigned-bitfields -O")
# -Wall -Wextra -ansi -fsigned-char -fsigned-bitfields -pthread -fPIC
SET( CFLAGS "-ansi -Wall -Wextra -fsigned-char -fsigned-bitfields" )


# 
#    None (CMAKE_C_FLAGS or CMAKE_CXX_FLAGS used)
#    Debug (CMAKE_C_FLAGS_DEBUG or CMAKE_CXX_FLAGS_DEBUG)
#    Release (CMAKE_C_FLAGS_RELEASE or CMAKE_CXX_FLAGS_RELEASE)
#    RelWithDebInfo (CMAKE_C_FLAGS_RELWITHDEBINFO or CMAKE_CXX_FLAGS_RELWITHDEBINFO
#    MinSizeRel (CMAKE_C_FLAGS_MINSIZEREL or CMAKE_CXX_FLAGS_MINSIZEREL)
#

if ( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" )
  message( "GNU C compiler" )
  SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS}")
#  SET( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CFLAGS}")
#  SET( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAG_RELEASES} ${CFLAGS}")
#  SET( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CFLAGS}")
#  SET( CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} ${CFLAGS}")
#  message( "${CMAKE_C_FLAGS}" )
endif( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" )

if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  message( "GNU CXX compiler" )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CFLAGS}")
#  SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CFLAGS}")
#  SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAG_RELEASES} ${CFLAGS}")
#  SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${CFLAGS}")
#  SET( CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${CFLAGS}")
endif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )

