/****************************************************
 * topological-voxel-tree.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 11:10:49 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_voxel_tree_h_
#define _topological_voxel_tree_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>
#include <topological-component-image.h>

extern void setVerboseInTopologicalVoxelTree( int v );
extern void incrementVerboseInTopologicalVoxelTree( );
extern void decrementVerboseInTopologicalVoxelTree( );

extern void setDebugInTopologicalVoxelTree( int d );
extern void incrementDebugInTopologicalVoxelTree(  );
extern void decrementDebugInTopologicalVoxelTree(  );





typedef enum enumVoxelTreeComponent {
  _INTERN_TREE_UNKNOWN_COMPONENT_,
  _INTERN_TREE_REMOVED_COMPONENT_,
  _INTERN_TREE_MERGED_COMPONENT_,
  _INTERN_TREE_INNER_EDGE_,
  _INTERN_TREE_END_EDGE_,
  _INTERN_TREE_EDGE_,
  _INTERN_TREE_JUNCTION_
} enumVoxelTreeComponent;



typedef struct voxelTreeIntPoint {
  int x;
  int y;
  int z;
  int anchor;
  int order_branch;
  int label_branch;
} voxelTreeIntPoint;

typedef struct voxelTreeFloatPoint {
  float x;
  float y;
  float z;
  int anchor;
  int order_branch;
  int label_branch;
} voxelTreeFloatPoint;



typedef struct infoBranchingPoint {
  int label_branch;
  int order;
  float length;
  int inqueue;
} infoBranchingPoint;



typedef struct voxelTreeComponent {

  enumVoxelTreeComponent type;

  int n_data;
  int n_allocated_data;
  /* branches: points are ordered from one endpoint to the other
   * junctions: all points are gathered but without any order
   */
  voxelTreeIntPoint *data;
  /* junctions: contains the barycenter of junctions points
   * branches: not applicable
   */
  voxelTreeFloatPoint center;

  /* bounding box
   */
  int minCorner[3];
  int maxCorner[3];

  /* for branches:
   * these two fields contain label of neighboring structures
   * (should be junctions), or -1 if there is no neighbors
   */
  int indexFirstJunction;
  int indexLastJunction;
  float length;

  /* to compute branch order
   * the two first are for edges
   * the last one for junction
   */
  infoBranchingPoint infoFirstPoint;
  infoBranchingPoint infoLastPoint;
  infoBranchingPoint infoJunction;

  /* this field is used for pruning
   */
  int canbepruned;

  /* this field is used for further computation
   * ie symbolic tree building
   * see treeImageToTree()
   */
  int pointIndexOffset;
} voxelTreeComponent;



typedef struct typeVoxelTree {
  voxelTreeComponent *data;
  int n_data;
  int n_allocated_data;
  float voxelSize[3];
} typeVoxelTree;



extern void initVoxelTreeComponent( voxelTreeComponent *c );
extern void freeVoxelTreeComponent( voxelTreeComponent *c );

extern int addPointToVoxelTreeComponent( voxelTreeComponent *l, voxelTreeIntPoint *p );



extern void initVoxelTree( typeVoxelTree *t );
extern void freeVoxelTree( typeVoxelTree *t );

extern int getNewComponentFromVoxelTree( typeVoxelTree *t );



typedef struct typeVoxelTreeList {
  typeVoxelTree *data;
  int n_data;
  int n_allocated_data;
} typeVoxelTreeList;

extern void initVoxelTreeList( typeVoxelTreeList *l );
extern void freeVoxelTreeList( typeVoxelTreeList *l );
extern int allocVoxelTreeList( typeVoxelTreeList *l, int n );

extern void fprintEdgeVoxelTree( FILE *f, typeVoxelTree *tree, int i, int writebranches );

extern void fprintfVoxelTreeComponent( FILE *f, voxelTreeComponent *c );
extern void fprintVoxelTree( FILE *f, typeVoxelTree *tree, int writebranches );
extern void fprintfSummaryVoxelTree( FILE *f, typeVoxelTree *c, char *desc );

extern int extractJunctionEdgeVoxelTree( typeComponentImage *treeImage,
                                 voxelTreeComponent *c );
extern int extractEdgeVoxelTree( typeComponentImage *treeImage,
                                 voxelTreeComponent *c, voxelTreeIntPoint *firstPoint );
extern int voxelTreeJunctionCenter( voxelTreeComponent *c );
extern void lengthEdgeVoxelTree( typeVoxelTree *tree,
                    int index );


extern int componentImageToVoxelTree( typeComponentImage *treeImage,
                                   void *theMark,
                                   bufferType typeMark,
                                   typeVoxelTree *tree );

extern int labelBranchImageToBranches( typeComponentImage *treeImage,
                                       typeVoxelTree *tree );

extern int countEndEdgeInImage( void* theBuf, bufferType theType, int *theDim );


typedef enum enumVoxelTreeAttribute {
    _NO_ATTRIBUTE_,
    _POINT_NEIGHBORS_,
    _COMPONENT_LABEL_,
    _BRANCH_LABEL_,
    _BRANCH_ORDER_
} enumVoxelTreeAttribute;

extern int voxelTreeToImage( typeVoxelTree *iTree,
                             void *resBuf, bufferType resType, int *resDim,
                             enumVoxelTreeAttribute attribute );

extern int writeComponentLengthFromVoxelTree( char *name, typeVoxelTree *iTree, int fullwrite );

extern int writeBranchLengthFromVoxelTree( char *name, typeVoxelTree *iTree, int fullwrite );

#ifdef __cplusplus
}
#endif

#endif
