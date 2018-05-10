/****************************************************
 * topological-voxel-tree-order.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 15 jui 2017 09:29:15 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_voxel_tree_order_h_
#define _topological_voxel_tree_order_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>
#include <topological-component-image.h>
#include <topological-voxel-tree.h>


extern void setVerboseInTopologicalVoxelTreeOrder( int v );
extern void incrementVerboseInTopologicalVoxelTreeOrder( );
extern void decrementVerboseInTopologicalVoxelTreeOrder( );

extern void setDebugInTopologicalVoxelTreeOrder( int d );
extern void incrementDebugInTopologicalVoxelTreeOrder(  );
extern void decrementDebugInTopologicalVoxelTreeOrder(  );


extern int branchesVoxelTree( typeComponentImage *treeImage, typeVoxelTree *iTree );

extern int maxOrderVoxelTree( typeVoxelTree *iTree );

extern int branchesTreeImage( void* theBuf, bufferType theType,
                       void* ancBuf, bufferType ancType,
                       void* resBuf, bufferType resType,
                       int *theDim, float *theSize );

#ifdef __cplusplus
}
#endif

#endif
