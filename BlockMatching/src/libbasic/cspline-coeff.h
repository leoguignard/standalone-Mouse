/*************************************************************************
 * cspline-coeff.h - Cubic splines coefficients
 *
 * Copyright (c) INRIA 2017
 *
 * AUTHOR:
 * Alexis Roche
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  4 mai 2017 11:11:32 CEST
 *
 * ADDITIONS, CHANGES
 *	
 *	
 *	
 *
 */

#ifndef _cspline_coeff_h_
#define _cspline_coeff_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <typedefs.h>


extern void setVerboseInCsplineCoeff( int v );
extern void incrementVerboseInCsplineCoeff( );
extern void decrementVerboseInCsplineCoeff( );

typedef struct {

  /* dimension du buffer
     theDim[0] -> dimension selon X
     theDim[1] -> dimension selon Y
     theDim[2] -> dimension selon Z
  */
  int theDim[3];
  /* buffer contenant les coefficients
   */
  float *theCoeff;

} typeCSplineCoefficients;


extern void InitTypeCSplineCoefficients( typeCSplineCoefficients *t );
extern void FreeTypeCSplineCoefficients( typeCSplineCoefficients **t );

extern typeCSplineCoefficients *ComputeCSplineCoefficients( void *theBuf,
						     bufferType theType,
						     int *theDim );




#ifdef __cplusplus
}
#endif

#endif
