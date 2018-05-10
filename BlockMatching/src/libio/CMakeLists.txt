####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(libio)
set(LIB_NAME io)

## #################################################################
## Input Files 
## #################################################################
# Create file list in "${SRC_FILES}" var with "${CMAKE_CURRENT_SOURCE_DIR}/*.c" regex
# file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

SET(SRC_FILES  ImageIO.c
	       analyze.c
	       bmp.c
	       bmpendian.c
	       bmpread.c
	       gif.c
	       inr.c
	       iris.c
               metaImage.c
	       nachos.c
	       pnm.c
	       raw.c )

if (  KLB_FOUND ) 
  SET( SRC_FILES ${SRC_FILES}
          klb.c
  )
else( KLB_FOUND )
  message( WARNING, " klb was not found" )
endif( KLB_FOUND )

if (  NIFTI_FOUND ) 
  SET( SRC_FILES ${SRC_FILES}
          nifti.c
  )
else( NIFTI_FOUND )
  message( WARNING, " nifti was not found" )
endif( NIFTI_FOUND )

if (  LIBTIFF_FOUND ) 
  SET( SRC_FILES ${SRC_FILES}
          tif.c
  )
else( LIBTIFF_FOUND )
  message( WARNING, " libtiff was not found (from src/libio)" )
endif( LIBTIFF_FOUND )


# Exe names
SET(EXE_NAMES test-libio
)
  

## #################################################################
## Build
## #################################################################
# Build library
add_library( ${LIB_NAME} ${SRC_FILES} )

# Add dependency to zlib
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} )

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} io)
endforeach(E)

