/*************************************************************************
 * bal-transformation-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_TOOLS_H
#define BAL_TRANSFORMATION_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string-tools.h>

#include <bal-transformation.h>
#include <bal-field.h>
#include <bal-estimator.h>



typedef enum enumTransformationInterpolation {
  NEAREST,
  LINEAR,
  CSPLINE
} enumTransformationInterpolation;



extern void BAL_SetVerboseInBalTransformationTools( int v );
extern void BAL_IncrementVerboseInBalTransformationTools(  );
extern void BAL_DecrementVerboseInBalTransformationTools(  );
extern void BAL_SetDebugInBalTransformationTools( int v );
extern void BAL_IncrementDebugInBalTransformationTools(  );
extern void BAL_DecrementDebugInBalTransformationTools(  );




/***************************************************
 *
 * transformation from one image to an other
 *
 ***************************************************/

extern int BAL_ComputeImageToImageTransformation( bal_image *subsampled_image,
                                                  bal_image *image_to_be_subsampled,
                                                  bal_transformation *subsampling_trsf );

extern int BAL_ComputeInitialTransformation( bal_image *refIm,
                                             bal_image *floIm,
                                             bal_transformation *theTrsf,
                                             enumInitialTransfo initial_transformation );



/***************************************************
 *
 * transformation from a pairing field
 *
 ***************************************************/

extern int BAL_ComputeIncrementalTransformation( bal_transformation *T,  
                                                 FIELD *field,
                                                 bal_estimator *estimator );

extern int BAL_ComputeTransformationResiduals( bal_transformation *T,  
                                               FIELD *field );



/***************************************************
 *
 * transformation use
 *
 ***************************************************/

/* point transformation
 */
extern int BAL_TransformFloatPoint( bal_floatPoint *thePt, bal_floatPoint *resPt, bal_transformation *theTr );

extern int BAL_TransformDoublePoint( bal_doublePoint *thePt, bal_doublePoint *resPt, bal_transformation *theTr );

/* image resampling
 */
extern int BAL_ResampleImage( bal_image *image, bal_image *resim, bal_transformation *theTr,
                              enumTransformationInterpolation interpolation );

extern int BAL_LinearResamplingCoefficients( bal_image *image, bal_image *resim,
                                             bal_transformation *theTr,
                                             enumTransformationInterpolation interpolation,
                                             int index );







/***************************************************
 *
 * transformation construction
 *
 ***************************************************/

extern void BAL_SetSinusoidAmplitudeForVectorFieldTransformation( double *a );
extern void BAL_SetSinusoidPeriodForVectorFieldTransformation( double *p );

extern int BAL_Sinusoid3DVectorField( bal_transformation *theTrsf );
extern int BAL_Sinusoid2DVectorField( bal_transformation *theTrsf );

/* random transformation 
 */

extern void BAL_SetMinAngleForRandomTransformation( double d );
extern void BAL_SetMaxAngleForRandomTransformation( double d );
extern void BAL_SetMinScaleForRandomTransformation( double d );
extern void BAL_SetMaxScaleForRandomTransformation( double d );
extern void BAL_SetMinShearForRandomTransformation( double d );
extern void BAL_SetMaxShearForRandomTransformation( double d );
extern void BAL_SetMinTranslationForRandomTransformation( double d );
extern void BAL_SetMaxTranslationForRandomTransformation( double d );

extern int BAL_Random2DTranslationMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DTranslationMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DTranslationScalingMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DTranslationScalingMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DRigidMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DRigidMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DSimilitudeMatrix( bal_transformation *theTrsf ) ;
extern int BAL_Random3DSimilitudeMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DAffineMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DAffineMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DVectorField( bal_transformation *theTrsf );
extern int BAL_Random3DVectorField( bal_transformation *theTrsf );
extern int BAL_SetTransformationToRandom( bal_transformation *t );

#ifdef __cplusplus
}
#endif

#endif
