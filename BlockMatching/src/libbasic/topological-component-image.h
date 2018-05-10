/****************************************************
 * topological-component-image.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Ven  9 jui 2017 16:05:53 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_component_image_
#define _topological_component_image_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>


extern void setVerboseInTopologicalComponentImage( int v );
extern void incrementVerboseInTopologicalComponentImage( );
extern void decrementVerboseInTopologicalComponentImage( );

extern void setDebugInTopologicalComponentImage( int d );
extern void incrementDebugInTopologicalComponentImage(  );
extern void decrementDebugInTopologicalComponentImage(  );


extern int countNeighborsInImage( void *theBuf, bufferType theType,
                                      void *resBuf, bufferType resType, int *theDim );


/* describes an image containing a tree
 * components (either junction or branch)
 * are labeled
 */

typedef enum enumComponentType {
  _UNKNOWN_COMPONENT_,
  _EDGE_COMPONENT_,
  _JUNCTION_COMPONENT_
} enumComponentType;

typedef struct typeComponent {
  int label;
  int size;
  enumComponentType type;
} typeComponent;

typedef struct typeComponentImage {
  /*
  int firstComponentLabel;
  int lastComponentLabel;
  int firstJunctionLabel;
  int lastJunctionLabel;
  */
  typeComponent *component;
  int n_component;
  int n_allocated_component;

  void *componentBuf;
  bufferType componentBufferType;
  int theDim[3];
  float voxelSize[3];
} typeComponentImage;

extern void initComponentImage( typeComponentImage *c );
extern int allocComponentImage( typeComponentImage *c,
                                    bufferType type, int *theDim );
extern void freeComponentImage( typeComponentImage *c );

extern int imageToComponentImage( void *theBuf, bufferType theType,
                                  typeComponentImage *componentImage );

#ifdef __cplusplus
}
#endif

#endif
