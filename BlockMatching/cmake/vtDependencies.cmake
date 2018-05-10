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
## ZLIB
## #################################################################
# Look for Z lib. Automatically defines following values : 
# ZLIB_FOUND         : Bool : true if find, else false
# ZLIB_INCLUDE_DIRS  : ZLIB HEADERS
# ZLIB_LIBRARIES     : FILE FOR LINKING WITH ZLIB
find_package( ZLIB REQUIRED )
include_directories( ${ZLIB_INCLUDE_DIRS} )

if(vt_USE_OPENMP)
  ## #################################################################
  # Look for OpenMP. Automatically definecs following values : 
  # OPENMP_FOUND     : Bool : true if find, else false
  # OpenMP_C_FLAGS   : OpenMP flags for C compiler
  # OpenMP_CXX_FLAGS : OpenMP flags for CXX compiler
  find_package( OpenMP )
  if (OPENMP_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

## #################################################################
## M Library
## #################################################################
# Unix specific : link the m library (automated on windows)
if(UNIX)
    link_libraries(m)
    link_libraries(pthread)
endif(UNIX)


## #################################################################
## NIFTI
## #################################################################

find_package (NIFTI QUIET)
if ( NIFTI_FOUND )
  message(STATUS "nifti was found")
  include(${NIFTI_USE_FILE}) 
  include_directories(${NIFTI_INCLUDE_DIRS})
  link_directories(${NIFTI_LIBRARY_DIRS})
  link_libraries(niftiio znz)
#  message( "" )
#  message( "NIFTI_LIBRARIES=${NIFTI_LIBRARIES}" )
#  message( "NIFTI_LIBRARY_DIRS=${NIFTI_LIBRARY_DIRS}" )
#  message( "NIFTI_INCLUDE_DIRS=${NIFTI_INCLUDE_DIRS}" )
#  message( "NIFTI_USE_FILE=${NIFTI_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_NIFTICLIB_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( NIFTI_FOUND )
  message(STATUS "nifti was NOT found")
endif( NIFTI_FOUND )



## #################################################################
## KLB
## #################################################################

find_package (KLB QUIET)
if ( KLB_FOUND )
  message(STATUS "klb was found")
  include(${KLB_USE_FILE}) 
  include_directories(${KLB_INCLUDE_DIRS})
  link_directories(${KLB_LIBRARY_DIRS})
  link_libraries(klb)
#  link_libraries(klb z bzip2)
#  message( "" )
#  message( "KLB_LIBRARIES=${KLB_LIBRARIES}" )
#  message( "KLB_LIBRARY_DIRS=${KLB_LIBRARY_DIRS}" )
#  message( "KLB_INCLUDE_DIRS=${KLB_INCLUDE_DIRS}" )
#  message( "KLB_USE_FILE=${KLB_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_KLB_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( KLB_FOUND )
  message(STATUS "klb was NOT found")
endif( KLB_FOUND )



## #################################################################
## LIBTIFF
## #################################################################

find_package (LIBTIFF QUIET)
if ( LIBTIFF_FOUND )
  message(STATUS "libtiff was found")
  include(${LIBTIFF_USE_FILE}) 
  include_directories(${LIBTIFF_INCLUDE_DIRS})
  link_directories(${LIBTIFF_LIBRARY_DIRS})
  link_libraries(tiff)
#  message( "" )
#  message( "LIBTIFF_LIBRARIES=${LIBTIFF_LIBRARIES}" )
#  message( "LIBTIFF_LIBRARY_DIRS=${LIBTIFF_LIBRARY_DIRS}" )
#  message( "LIBTIFF_INCLUDE_DIRS=${LIBTIFF_INCLUDE_DIRS}" )
#  message( "LIBTIFF_USE_FILE=${LIBTIFF_USE_FILE}" )
#  message( "" )
  SET( CPPFLAGS "${CPPFLAGS} -D_LIBLIBTIFF_" )
  add_definitions( ${CPPFLAGS} )
#  MESSAGE( "CPPFLAGS=${CPPFLAGS}" )
else( LIBTIFF_FOUND )
  message(STATUS "libtiff was NOT found (from cmake/vtDependencies)")
endif( LIBTIFF_FOUND )









