/****************************************************
 * topological-voxel-tree-main.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mar 20 jui 2017 11:38:23 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_voxel_tree_main_h_
#define _topological_voxel_tree_main_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <typedefs.h>
#include <topological-component-image.h>
#include <topological-voxel-tree.h>

extern void setVerboseInTopologicalVoxelTreeMain( int v );
extern void incrementVerboseInTopologicalVoxelTreeMain( );
extern void decrementVerboseInTopologicalVoxelTreeMain( );

extern void setDebugInTopologicalVoxelTreeMain( int d );
extern void incrementDebugInTopologicalVoxelTreeMain(  );
extern void decrementDebugInTopologicalVoxelTreeMain(  );

extern int relabelVoxelTree( typeComponentImage *treeImage,
                             void* ancBuf, bufferType ancType,
                             typeVoxelTree *iTree );

extern int buildVoxelTree( void* theBuf, bufferType theType,
                           void* ancBuf, bufferType ancType,
                           int *theDim, float *theSize,
                           typeVoxelTree *iTree,
                           float realBranchLength, int voxelBranchLength, int nmaxendedges,
                           int maxorder );

#ifdef __cplusplus
}
#endif

#endif
