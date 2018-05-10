/****************************************************
 * topological-growing.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 29 fév 2016 18:28:00 CET
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _topological_growing_h_
#define _topological_growing_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <chamferdistance.h>
#include <topological-operations-common.h>


extern void setVerboseInTopologicalGrowing( int v );
extern void incrementVerboseInTopologicalGrowing();
extern void decrementVerboseInTopologicalGrowing();







typedef struct typeGrowingParameters {
  int maxIteration;
  int connectivity;
  typeChamferMask *theMask;
  enumTypeSort additionalSorting;
} typeGrowingParameters;

extern void initGrowingParameters( typeGrowingParameters *p );




extern int initGrowingImageFromDistance( unsigned char *resBuf,
				      unsigned short *theDistance,
				      int *theDim );

extern int valueBasedGrowing( unsigned char *resBuf,
			    unsigned short *theDistance,
			    unsigned short *thePropagation,
			    int *theDim,
          typeGrowingParameters *par );


#ifdef __cplusplus
}
#endif

#endif
