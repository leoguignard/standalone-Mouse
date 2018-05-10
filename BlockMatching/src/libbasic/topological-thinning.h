/****************************************************
 * topological-thinning.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  2 mar 2016 22:02:15 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_thinning_h_
#define _topological_thinning_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>

extern void setVerboseInTopologicalThinning( int v );
extern void incrementVerboseInTopologicalThinning( );
extern void decrementVerboseInTopologicalThinning( );
extern void setDebugInTopologicalThinning( int v );
extern void incrementDebugInTopologicalThinning( );
extern void decrementDebugInTopologicalThinning( );
extern void setTimeInTopologicalThinning( int v );





typedef enum enumThickness {
  _04THICKNESS_ = 4,
  _08THICKNESS_ = 8,
  _06THICKNESS_ = 6,
  _18THICKNESS_ = 18,
  _26THICKNESS_ = 26
} enumThickness;
    


typedef enum enumTypeEndPoint {
  _SURFACE_,
  _PURE_SURFACE_,
  _CURVE_,
  _PURE_CURVE_,
  _NO_END_POINT_
} enumTypeEndPoint;



typedef enum enumTypeListOrder {
  _INCREASING_ORDER_,
  _DECREASING_ORDER_
} enumTypeListOrder;



typedef enum enumTypeChange {
  _FOREGROUND_TO_BACKGROUND_, /* thinning */
  _BACKGROUND_TO_FOREGROUND_  /* thickening */
} enumTypeChange;



typedef struct typeThinningParameters {
  /* list construction: bin size
   */
  int binLength;

  /* some control conditions
   * typically, for distance-based thinning we have
   *  - typeThickness = _26THICKNESS_
   *  - typeOrdering = _INCREASING_ORDER_
   *  - typeChanging = _FOREGROUND_TO_BACKGROUND_
   * for thickening, and distance computed from the object to be thickened
   *  - typeThickness = _06THICKNESS_
   *  - typeOrdering = _INCREASING_ORDER_
   *  - typeChanging = _BACKGROUND_TO_FOREGROUND_
   */
  enumTypeListOrder typeOrdering;
  enumTypeChange typeChanging;

  /* points to be changed
   */
  enumThickness typeThickness;
  enumTypeEndPoint typeEndPoint;

  /* end conditions are tested if
   * 1. valueBeforeEnding is reached (e.g. allows to avoid to have
   *    endpoints with small distances for distance-based thinning)
   * 2. the number of cycles (=loop over directions for a given
   *    value) is reached (more or less the same thing for plateaus
   *    of constant value)
   */
  int cyclesBeforeEnding;
  int valueBeforeEnding;

  /* another end condition is if the maximal number of iteration
   * is reached
   * one iteration is a cycle (loop over directions) for a given value
   */
  int maxIteration;

} typeThinningParameters;



extern void initTypeThinningParameters( typeThinningParameters *p );

extern int thinningThreshold( void *bufferIn, bufferType typeIn,
                              unsigned char *bufferOut, int *bufferDims,
                              float lowThreshold,
                              float highThreshold );

extern int chamferBasedThinning( unsigned char *theBuf,
                                 int *theDim,
                                 const int chamfer,
                                 typeThinningParameters *p );

extern int valueBasedThinning( unsigned char *theBuf,
                               void *theValues,
                               bufferType valuesType,
                               int *theDim,
                               typeThinningParameters *par );



#ifdef __cplusplus
}
#endif

#endif
