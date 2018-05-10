/*************************************************************************
 * bal-transformation-copy.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 22 jan 2017 14:45:30 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_COPY_H
#define BAL_TRANSFORMATION_COPY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-transformation.h>



extern void BAL_SetVerboseInBalTransformationCopy( int v );\
extern void BAL_IncrementVerboseInBalTransformationCopy( );
extern void BAL_DecrementVerboseInBalTransformationCopy( );
extern void BAL_SetDebugInBalTransformationCopy( int v );
extern void BAL_IncrementDebugInBalTransformationCopy( );
extern void BAL_DecrementDebugInBalTransformationCopy( );



/***************************************************
 *
 * transformation copy
 *
 ***************************************************/

extern int BAL_CopyTransformation( bal_transformation *theTrsf,
                                   bal_transformation *resTrsf );
extern int BAL_CopyTransformationList( bal_transformationList *theTrsfs,
                                bal_transformationList *resTrsfs );



/***************************************************
 *
 * transformation conversion
 *
 ***************************************************/
/* unit conversion

   'theTrsf' is the transformation that goes from 'resim' to 'image'
   ie allows to resample 'image' into the geometry of 'resim'.
   When 'theTrsf' is the transformation issued from matching,
   'resim' is the reference image and 'image' the floating image.

*/

extern int BAL_ChangeTransformationToRealUnit( bal_image *image,
                                               bal_image *resim,
                                               bal_transformation *theTrsf,
                                               bal_transformation *resTrsf );

extern int BAL_ChangeTransformationToVoxelUnit( bal_image *image,
                                                bal_image *resim,
                                                bal_transformation *theTrsf,
                                                bal_transformation *resTrsf );

#ifdef __cplusplus
}
#endif

#endif
