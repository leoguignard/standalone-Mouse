/*************************************************************************
 * cspline.h - Cubic splines
 *
 * $Id: cspline.h,v 1.2 2000/10/23 14:32:34 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Alexis Roche
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Thu Oct 12 12:08:12 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *	
 *	
 *	
 *
 */

#ifndef _cspline_h_
#define _cspline_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>
#include <cspline-coeff.h>


extern void setVerboseInCspline( int v );
extern void incrementVerboseInCspline( );
extern void decrementVerboseInCspline( );




extern int Reech3DCSpline4x4WithCoefficients( typeCSplineCoefficients *theCoeff,
		       r32 *resBuf,  /* result buffer */
		       int *resDim,  /* dimensions of this buffer */
		       double *mat,
		       int slice, int *derivative );

extern int Reech2DCSpline4x4WithCoefficients( typeCSplineCoefficients *theCoeff,
		       r32* resBuf,  /* result buffer */
		       int *resDim,  /* dimensions of this buffer */
		       double *mat,
		       int *derivative );

extern int ReechCSpline4x4( void* theBuf, bufferType theType, int *theDim,
			      void* resBuf, bufferType resType, int *resDim,
			      double *mat,
			      int *derivative );




#ifdef __cplusplus
}
#endif

#endif
