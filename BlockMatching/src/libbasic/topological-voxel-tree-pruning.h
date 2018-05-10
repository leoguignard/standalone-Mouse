/****************************************************
 * topological-voxel-tree-pruning.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 11:30:52 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_voxel_tree_pruning_h_
#define _topological_voxel_tree_pruning_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>

extern void setVerboseInTopologicalVoxelTreePruning( int v );
extern void incrementVerboseInTopologicalVoxelTreePruning( );
extern void decrementVerboseInTopologicalVoxelTreePruning( );

extern void setDebugInTopologicalVoxelTreePruning( int d );
extern void incrementDebugInTopologicalVoxelTreePruning(  );
extern void decrementDebugInTopologicalVoxelTreePruning(  );

extern int pruneVoxelTreeWithOrder( typeComponentImage *treeImage, typeVoxelTree *iTree,
                                    int maxorder );

extern int pruneVoxelTreeWithEdgeLength( typeComponentImage *treeImage, typeVoxelTree *iTree,
                           float realBranchLength, int voxelBranchLength,
                           int nmaxendedges );

extern int pruneTreeImage( void* theBuf, bufferType theType,
                           void* ancBuf, bufferType ancType,
                           void* resBuf, bufferType resType,
                           int *theDim, float *theSize,
                           float realBranchLength, int voxelBranchLength, int nmaxendedges );



#ifdef __cplusplus
}
#endif

#endif
