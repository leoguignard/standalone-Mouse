/****************************************************
 * topological-generic-tree.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 10:28:22 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_generic_tree_h_
#define _topological_generic_tree_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>
#include <topological-component-image.h>
#include <topological-voxel-tree.h>



extern void setVerboseInTopologicalGenericTree( int v );
extern void incrementVerboseInTopologicalGenericTree( );
extern void decrementVerboseInTopologicalGenericTree( );

extern void setDebugInTopologicalGenericTree( int d );
extern void incrementDebugInTopologicalGenericTree(  );
extern void decrementDebugInTopologicalGenericTree(  );





typedef enum _enumTreePointType {
  _TREE_UNKNOWN_POINT_,
  _TREE_END_POINT_,
  _TREE_EDGE_POINT_,
  _TREE_JUNCTION_POINT_
} _enumTreePointType;

typedef struct typeTreeFloatPoint {
  _enumTreePointType type;
  float x;
  float y;
  float z;
  int anchor; /* to mark point from main axon */
  int order_branch;
  int label_branch;
  int label_component;
  int color;  /* for output */
} typeTreeFloatPoint;



typedef enum enumTreeEdgeType {
  _UNDEFINED_EDGE_,
  _INNER_EDGE_,
  _END_EDGE_
} enumTreeEdgeType;

typedef struct typeTreeEdge {
  enumTreeEdgeType type;
  int label;
  int n_data;
  int n_allocated_data;
  int *data;
} typeTreeEdge;



typedef struct typeTree {
  int n_points;
  int n_allocated_points;
  typeTreeFloatPoint *points;
  int n_edges;
  int n_allocated_edges;
  typeTreeEdge *edges;
} typeTree;



extern void initTree( typeTree *t );
extern void freeTree( typeTree *t );


extern int voxelTreeToTree( typeVoxelTree *iTree, typeTree *tree );

extern int treeImageToTree( typeComponentImage *treeImage,
                            void* ancBuf, bufferType ancType,
                            typeTree *t );



typedef enum enumFileType {
  _MULTIPLE_FILES_,
  _SINGLE_FILE_
} enumFileType;

typedef enum enumTreePointColor {
  _NO_COLOR_,
  _COLOR_ANCHOR_,
  _COLOR_BRANCH_ORDER_,
  _COLOR_BRANCH_LABEL_,
  _COLOR_COMPONENT_LABEL_
} enumTreePointColor;

typedef enum enumColorTable {
  _DEFAULT_COLORTABLE_,
  _USER_COLORTABLE_
} enumColorTable;

typedef struct typeTreeOutputParam {
  enumTreePointColor color;
  enumFileType filetype;  /* one file per branch (multiples files)
                             or one file for the whole tree */
  enumColorTable colortable;
  float table[2][3];
} typeTreeOutputParam;



extern void initTreeOutputParam( typeTreeOutputParam *p );
extern void fprintfTreeOutputParam( FILE *f, typeTreeOutputParam *p, char *str );

extern int writeVTKLegacyFile( char *name, typeTree *tree,
                               typeTreeOutputParam *par,
                               char *desc );


#ifdef __cplusplus
}
#endif

#endif
