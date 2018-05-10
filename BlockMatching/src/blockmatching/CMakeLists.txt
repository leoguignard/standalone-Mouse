####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2012
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(blockmatching)
set(LIB_NAME blockmatchingLIB)

message(STATUS "Blockmatching : LIB suffix added to the lib") 
# message(WARNING "Blockmatching : LIB suffix added to the lib") 


## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  
#        api-execTemplate.c
        api-applyTrsfToPoints.c
        api-applyTrsf.c
        api-blockmatching.c
	api-buildPyramidImage.c
        api-composeTrsf.c
        api-copyTrsf.c
        api-createGrid.c
        api-createTrsf.c
        api-cropImage.c
        api-intermediaryTrsf.c
        api-interpolateImages.c
        api-invTrsf.c
        api-pointmatching.c
        bal-behavior.c
	bal-block-tools.c
	bal-block.c
 	bal-blockmatching-param-tools.c
	bal-blockmatching-param.c
	bal-blockmatching.c
	bal-estimator.c
	bal-field-tools.c 
	bal-field.c
	bal-image-tools.c
	bal-image.c
        bal-interpolation.c
	bal-lineartrsf-tools.c
	bal-lineartrsf.c
	bal-matrix.c
        bal-point.c
	bal-pyramid.c
	bal-stddef.c
	bal-tests.c
	bal-transformation-compose.c
	bal-transformation-copy.c
	bal-transformation-inversion.c
	bal-transformation-tools.c
	bal-transformation.c
	bal-vectorfield.c
	)
	

# Exe names
SET(EXE_NAMES 
#        execTemplate
        applyTrsf
        applyTrsfToPoints
        blockmatching
	buildPyramidImage
        composeTrsf
        copyTrsf
        createGrid
        createTrsf
        cropImage
        intermediaryTrsf
        interpolateImages
        invTrsf
        pointmatching
        printImage
        printTrsf
#	test-copy
)
  
## #################################################################
## Build
## #################################################################

## Add dependencies on object
#set_source_files_properties(
#  bal-lineartrsf.c 
#  PROPERTY 
#  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bal-lineartrsf-wls2D.c
#	OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bal-lineartrsf-wls3D.c
#	OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bal-lineartrsf-ls2D.c
#	OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bal-lineartrsf-ls3D.c
#)

# Add dependency on object
#set_source_files_properties(
#  bal-vectorfield.c 
#  PROPERTY
#  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bal-inverseVectorField.c
#)

# build library and link
add_library(${LIB_NAME} ${LIB_STATUS} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io)

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io)
endforeach(E)


## #################################################################
# Build valgrind. history.
## #################################################################
#
#SET(VALGRIND_NAMES blockmatching)
#
#SET(ALL_VALGRIND_NAMES test-copy
#   applyTrsf
#   copyTrsf
#   createGrid
#   createRandomTrsf
#   createVectorTrsf
#   invTrsf
#   printTrsf
#   composeTrsf
#   blockmatching
#  )
#
#if( BLOCKMATCHING_BUILD_VALGRIND )
#  foreach(V ${VALGRIND_NAMES})
#    add_executable(${NAME}_valgrind-${V} ${CMAKE_CURRENT_SOURCE_DIR}/${V}.c)
#    target_link_libraries(${NAME}_valgrind-${V} ${NAME} ${ZLIB_LIBRARIES} basic io)
#    set_target_properties(${NAME}_valgrind-${V} PROPERTIES OUTPUT_NAME ${V})
#  endforeach(V)
#endif( BLOCKMATCHING_BUILD_VALGRIND )
#
#if( BLOCKMATCHING_BUILD_ALL_VALGRIND )
#  foreach(AV ${ALL_VALGRIND_NAMES})
#    add_executable(${NAME}_allvalgrind-${AV} ${CMAKE_CURRENT_SOURCE_DIR}/${AV}.c)
#    target_link_libraries(${NAME}_allvalgrind-${AV} ${NAME} ${ZLIB_LIBRARIES} basic io)
#    set_target_properties(${NAME}_allvalgrind-${AV} PROPERTIES OUTPUT_NAME ${AV})
#  endforeach(AV)
#endif( BLOCKMATCHING_BUILD_ALL_VALGRIND )

