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
#

# Name
project(libbasic)
set(LIB_NAME basic)

## #################################################################
## LIB, BIN : headers, sources
## #################################################################
# Src files
SET(SRC_FILES basictransform.c
	chamfercoefficients.c
	chamferdistance.c
	chamferdistance-mask.c
	chunks.c
	connexe.c
	convert.c
	convolution.c
	convolution1D.c
	cspline-coeff.c
	cspline.c
	cspline4x4.c
	drawShapes.c
	eigens.c
	extract.c
	file-tools.c
	histogram.c
	iopnm.c
	issimple3D.c
	levenberg.c
	linearFiltering-common.c
	linearFiltering-contours.c
	linearFiltering.c
	local-operation.c
	morphotools.c
	parcelling.c
	pixel-operation.c
	recbuffer.c
	recline.c
	reech-def.c
        reech4x4-coeff.c
        reech4x4.c
	regionalext.c
	sour-filter.c
	string-tools.c
	systlin.c
	t04t08.c
	t06t26.c
        threshold.c
	topological-component-image.c
	topological-generic-tree.c
        topological-growing.c
        topological-operations-common.c
        topological-thinning.c
        topological-voxel-tree.c
	topological-voxel-tree-main.c
	topological-voxel-tree-order.c
	topological-voxel-tree-pruning.c
	transfo.c
        vtmalloc.c
	watershed.c
 	zcross.c  
 	morpho.c
  )

SET(EXE_NAMES
#    test-chamferdistance-mask
#    test-chunks-loop
#    test-chunks
#    test-connexe-pnm
#    test-convol-pnm
#    test-edges-pnm
#    test-edges
#    test-eigens
#    test-file-tools
#    test-hyster-pnm
#    test-iopnm
#    test-morpho
#    test-openmp
#    test-quater
#    test-reech-pnm
#    test-rot
#    test-schedule-openmp
#    test-systlin
#    test-thinning-pnm
#    test-transfo-spherical
#    test-vtmalloc
)



## #################################################################
## Build 
## #################################################################
# Add dependency on object
set_source_files_properties(
  linearFiltering.c 
  PROPERTY
  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linearFiltering-gradient.c
  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linearFiltering-hessian.c
  OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linearFiltering-laplacian.c
)

# Build lib
add_library(${LIB_NAME} ${SRC_FILES})

# Build exe
foreach(E ${EXE_NAMES})
  
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME})
  
endforeach(E)


## #################################################################
## ADD BUILD_TEMPLATES INSTRUCTIONS - not cross-platform
## #################################################################
# BASIC_BUILD_TEMPLATES : fixed via ccmake 
if(BASIC_BUILD_TEMPLATES)
  # execute non cross-platform process generateTemplates.sh in templates-gneration dir
  execute_process( COMMAND chmod +x ${CMAKE_CURRENT_SOURCE_DIR}/templates-generation/generateTemplates.sh )
  execute_process( COMMAND ./generateTemplates.sh
		   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/templates-generation )
endif(BASIC_BUILD_TEMPLATES)


  
## #################################################################
## Test : history.
## #################################################################
# Build test (cached var : determine via ccmake)
if(${BUILD_TESTING})
  foreach(T ${TEST_NAMES})
  
    add_executable(${T} ${CMAKE_CURRENT_SOURCE_DIR}/${T}.c)
    target_link_libraries(${T} ${LIB_NAME})
    add_test(${T} ${T})
    
  endforeach(T)
endif(${BUILD_TESTING})


####################################################################
## IOC Files (unused)
## #################################################################
# Unused.. Hystory. 
# SET(IOC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# SET(IOC_H_DIR   ${CMAKE_CURRENT_SOURCE_DIR})

# Src files
#SET(IOC_FILES ${IOC_SRC_DIR}/bmpendian.c
#	${IOC_SRC_DIR}/iobmp.c
#	${IOC_SRC_DIR}/iogif.c
#	${IOC_SRC_DIR}/readbmp.c)


