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
## Header
## #################################################################
# Specify the version being used
cmake_minimum_required(VERSION 2.8)


# from VTK-6.2.0 CMakeLists.txt
# Set a consistent MACOSX_RPATH default across all CMake versions.
# When CMake 2.8.12 is required, change this default to 1.
# When CMake 3.0.0 is required, remove this block (see CMP0042).
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

# Project name
project(vt)

## #################################################################
## Options
## #################################################################
# Specify the version being used
option(BUILD_SHARED_LIBS      " Build shared libraries"  true)
option(BASIC_BUILD_TEMPLATES  " Regenerate templated files in basic lib" false)
option(BUILD_TESTING	       " Build testing" false)
option(vt_USE_OPENMP	       " Use OpenMP" true) 

SET(KLB_DIR              /media/simview/data/Leo/latest-comp/Python-Scripts/BlockMatching/external/KLBFILE/build)
SET(LIBTIFF_DIR             /media/simview/data/Leo/latest-comp/Python-Scripts/BlockMatching/external/TIFF/build)
# SET(FILTERS_BUILD_TEST CACHE BOOL ON)
# SET(FILTERS_BUILD_INSIGHT CACHE BOOL ON)
# SET(FILTERS_BUILD_PURIFY CACHE BOOL ON)
# SET(FILTERS_BUILD_QUANTIFY CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_COMPARE CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_TIME CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_STATS CACHE BOOL ON)


## #################################################################
## Config file
## #################################################################
# Cmake config file directory path
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Call include vtBuild.cmake file. This file configure builds ouput (build/bin, build/lib)
include(vtBuild)

# Call include vtDependencies.cmake file. This file configure the externals dependencies of vt (zlib, openmp, m)
include(vtDependencies)

include(vtPack)

# To improve
SET(vt_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src/blockmatching
  ${PROJECT_SOURCE_DIR}/src/libbasic
  ${PROJECT_SOURCE_DIR}/src/libio
  )


set(vt_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set(vt_CMAKE_DIR    ${CMAKE_MODULE_PATH})
set(vt_USE_FILE     ${PROJECT_BINARY_DIR}/vtUse.cmake)

## Setup Config File
if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/vtConfig.cmake.in)
  configure_file( ## Build tree configure file
    ${PROJECT_SOURCE_DIR}/cmake/vtConfig.cmake.in
    ${PROJECT_BINARY_DIR}/vtConfig.cmake
    @ONLY IMMEDIATE)
endif()

## Setup use file
if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/vtUse.cmake.in)
configure_file( ## Common use file
  ${PROJECT_SOURCE_DIR}/cmake/vtUse.cmake.in
  ${PROJECT_BINARY_DIR}/vtUse.cmake
  @ONLY IMMEDIATE)
endif()


## #################################################################
## ADD_SUBIDIRECTORY
## #################################################################

# Source directory
add_subdirectory(src)
