/*************************************************************************
 * bal-transformation-inversion.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Sam 18 fév 2017 14:11:50 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_INVERSION_H
#define BAL_TRANSFORMATION_INVERSION_H

#ifdef __cplusplus
extern "C" {
#endif


#include <bal-transformation.h>


typedef enum enumVectorFieldInverseInitialization {
  ZERO,
  FORWARD_INTERPOLATION,
} enumVectorFieldInverseInitialization;


extern void BAL_SetVerboseInBalTransformationInversion( int v );
extern void BAL_IncrementVerboseInBalTransformationInversion(  );
extern void BAL_DecrementVerboseInBalTransformationInversion(  );
extern void BAL_SetDebugInBalTransformationInversion( int v );
extern void BAL_IncrementDebugInBalTransformationInversion(  );
extern void BAL_DecrementDebugInBalTransformationInversion(  );


/***************************************************
 *
 * transformation inversion
 *
 ***************************************************/

extern void BAL_SetDerivationSigmaForVectorFieldInversionInBalTransformationInversion( double s );
extern double BAL_GetDerivationSigmaForVectorFieldInversionInBalTransformationInversion();
extern void BAL_SetIterationsMaxForVectorFieldInversionInBalTransformationInversion( int i );
extern int BAL_GetIterationsMaxForVectorFieldInversionInBalTransformationInversion( );
extern void BAL_SetErrorMaxForVectorFieldInversionInBalTransformationInversion( double e );
extern double BAL_GetErrorMaxForVectorFieldInversionInBalTransformationInversion( );
extern void BAL_SetInitializationForVectorFieldInversionInBalTransformationInversion( enumVectorFieldInverseInitialization i );
extern enumVectorFieldInverseInitialization BAL_GetInitializationForVectorFieldInversionInBalTransformationInversion();
extern void BAL_SetForwardSigmaForVectorFieldInversionInBalTransformationInversion( double s );
extern double BAL_GetForwardSigmaForVectorFieldInversionInBalTransformationInversion();

extern void BAL_SetImageInverseErrorsForVectorFieldInversionInBalTransformationInversion( bal_image *i );

extern void BAL_SetConversionForVectorFieldInversionInBalTransformationInversion( int c );

extern int BAL_InverseTransformation( bal_transformation *theTrsf,
                                      bal_transformation *invTrsf );


#ifdef __cplusplus
}
#endif

#endif
