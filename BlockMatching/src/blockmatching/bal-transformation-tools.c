/*************************************************************************
 * bal-transformation-tools.c -
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



#include <stdlib.h>
#include <stdio.h>
#include <math.h>



#include <chunks.h>
#include <cspline.h>
#include <reech4x4.h>
#include <reech4x4-coeff.h>
#include <reech-def.h>
#include <vtmalloc.h>

#include <bal-transformation-copy.h>
#include <bal-lineartrsf.h>
#include <bal-lineartrsf-tools.h>
#include <bal-vectorfield.h>

#include <bal-transformation-tools.h>

static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationTools( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationTools(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationTools(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationTools( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalTransformationTools(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationTools(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}









/*--------------------------------------------------
 *
 * TRANSFORMATION COMPUTATION
 *
 --------------------------------------------------*/




/* compute the subsampling matrix T (in 'real' coordinates)
   that allows to subsample the image 'image_to_be_subsampled'
   into the geometry of 'subsampled_image'

   It goes then from 'subsampled_image' to 'image_to_be_subsampled'
   and we have 'subsampled_image = image_to_be_subsampled o T'

   T = T {image_to_be_subsampled <- subsampled_image}
     = T {flo <- ref}

   C_flo = T o C_ref = C_ref + TRS
   => TRS = C_flo - C_ref

   It puts the center of the field of view of 'subsampled_image' on the 
   center of the field of view of 'image_to_be_subsampled'.
   Assuming that the origin is at the center of the upper left voxel, 
   the coordinates of the center of the FOV are 
   ((dx-1)/2, (dy-1)/2, (dz-1)/2 )
  
*/
/* made compliant to QForm matrix
 */
int BAL_ComputeImageToImageTransformation( bal_image *subsampled_image,
                                           bal_image *image_to_be_subsampled,
                                           bal_transformation *subsampling_trsf )
{
  char *proc ="BAL_ComputeImageToImageTransformation";
  _MATRIX mat;
  size_t i, v;
  bal_image *refIm = subsampled_image;
  bal_image *floIm = image_to_be_subsampled;
  double refCtr[3], refTrsfedCtr[3];
  double floCtr[3], floTrsfedCtr[3];
  float *vx, *vy, *vz;

  if ( _alloc_mat( &mat, 4, 4) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error when allocating matrix\n", proc );
    return( -1 );
  }
  for ( i=0; i<16; i++ ) mat.m[i] = 0.0;
  mat.m[0] = mat.m[5] = mat.m[10] = mat.m[15] = 1.0;

  refCtr[0] = (double)(refIm->ncols - 1) / 2.0;
  refCtr[1] = (double)(refIm->nrows - 1) / 2.0;
  refCtr[2] = (double)(refIm->nplanes - 1) / 2.0;

  refTrsfedCtr[0] = refIm->to_real.m[ 0] * refCtr[0]
                  + refIm->to_real.m[ 1] * refCtr[1]
                  + refIm->to_real.m[ 2] * refCtr[2]
                  + refIm->to_real.m[ 3];
  refTrsfedCtr[1] = refIm->to_real.m[ 4] * refCtr[0]
                  + refIm->to_real.m[ 5] * refCtr[1]
                  + refIm->to_real.m[ 6] * refCtr[2]
                  + refIm->to_real.m[ 7];
  refTrsfedCtr[2] = refIm->to_real.m[ 8] * refCtr[0]
                  + refIm->to_real.m[ 9] * refCtr[1]
                  + refIm->to_real.m[10] * refCtr[2]
                  + refIm->to_real.m[11];

  floCtr[0] = (double)(floIm->ncols - 1) / 2.0;
  floCtr[1] = (double)(floIm->nrows - 1) / 2.0;
  floCtr[2] = (double)(floIm->nplanes - 1) / 2.0;

  floTrsfedCtr[0] = floIm->to_real.m[ 0] * floCtr[0]
                  + floIm->to_real.m[ 1] * floCtr[1]
                  + floIm->to_real.m[ 2] * floCtr[2]
                  + floIm->to_real.m[ 3];
  floTrsfedCtr[1] = floIm->to_real.m[ 4] * floCtr[0]
                  + floIm->to_real.m[ 5] * floCtr[1]
                  + floIm->to_real.m[ 6] * floCtr[2]
                  + floIm->to_real.m[ 7];
  floTrsfedCtr[2] = floIm->to_real.m[ 8] * floCtr[0]
                  + floIm->to_real.m[ 9] * floCtr[1]
                  + floIm->to_real.m[10] * floCtr[2]
                  + floIm->to_real.m[11];


  mat.m[ 3] = floTrsfedCtr[0] - refTrsfedCtr[0];
  mat.m[ 7] = floTrsfedCtr[1] - refTrsfedCtr[1];
  mat.m[11] = floTrsfedCtr[2] - refTrsfedCtr[2];

  switch ( subsampling_trsf->type ) {

  default :

    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    if ( subsampling_trsf->mat.m ==  (double*)NULL ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: matrix is not allocated\n", proc );
      return( -1 );
    }

    for ( i=0; i<16; i++ ) 
      subsampling_trsf->mat.m[i] = mat.m[i];
    break;

  case VECTORFIELD_2D :

    if ( subsampling_trsf->vx.ncols != subsampled_image->ncols
         || subsampling_trsf->vx.nrows != subsampled_image->nrows
         || subsampling_trsf->vx.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: X component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    if ( subsampling_trsf->vy.ncols != subsampled_image->ncols
         || subsampling_trsf->vy.nrows != subsampled_image->nrows
         || subsampling_trsf->vy.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: Y component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    
    v = subsampling_trsf->vx.ncols * subsampling_trsf->vx.nrows * subsampling_trsf->vx.nplanes;
    vx = (float*)(subsampling_trsf->vx.data);
    vy = (float*)(subsampling_trsf->vy.data);

    for ( i=0; i<v; i++ ) {
      vx[i] =  mat.m[ 3];
      vy[i] =  mat.m[ 7];
    }
    break;

  case VECTORFIELD_3D :

    if ( subsampling_trsf->vx.ncols != subsampled_image->ncols
         || subsampling_trsf->vx.nrows != subsampled_image->nrows
         || subsampling_trsf->vx.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: X component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    if ( subsampling_trsf->vy.ncols != subsampled_image->ncols
         || subsampling_trsf->vy.nrows != subsampled_image->nrows
         || subsampling_trsf->vy.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: Y component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    
    if ( subsampling_trsf->vz.ncols != subsampled_image->ncols
         || subsampling_trsf->vz.nrows != subsampled_image->nrows
         || subsampling_trsf->vz.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: Z component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }

    v = subsampling_trsf->vx.ncols * subsampling_trsf->vx.nrows * subsampling_trsf->vx.nplanes;
    vx = (float*)(subsampling_trsf->vx.data);
    vy = (float*)(subsampling_trsf->vy.data);
    vz = (float*)(subsampling_trsf->vz.data);
    for ( i=0; i<v; i++ ) {
      vx[i] =  mat.m[ 3];
      vy[i] =  mat.m[ 7];
      vz[i] =  mat.m[11];
    }
    break;

  }
  _free_mat( &mat );

  /* the transformation has been calculated in real units
   */
  subsampling_trsf->transformation_unit = REAL_UNIT;

  return ( 1 );
}





int BAL_ComputeInitialTransformation( bal_image *refIm,
                                           bal_image *floIm,
                                           bal_transformation *theTrsf,
                                           enumInitialTransfo initial_transformation )
{
  char *proc = "BAL_ComputeImageToImageTransformation";

  switch( initial_transformation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such initial transformation type not handled yet\n", proc );
    return( -1 );
  case _BAL_IDENTITY_TRANSFORMATION_ :
    BAL_SetTransformationToIdentity( theTrsf );
    break;
  case _BAL_FOVCENTER_TRANSFORMATION_ :
    if ( BAL_ComputeImageToImageTransformation( refIm, floIm, theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute FOV center based transformation\n", proc );
      return( -1 );
    }
    break;
  }

  return( 1 );
}






/*--------------------------------------------------
 *
 * transformation from a pairing field
 *
 --------------------------------------------------*/



int BAL_ComputeIncrementalTransformation( bal_transformation *T,  FIELD *field, 
                                          bal_estimator *estimator )

{
  char *proc = "BAL_ComputeIncrementalTransformation";
  
  switch( T->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    /* le champ est en voxel,
       avec le point de depart dans l'image de reference
       on le passe en reel (important pour la transformation rigide)
    */
    BAL_ChangeFieldToRealUnit( field );
    
    if ( BAL_ComputeLinearTransformation( T, field, estimator ) <= 0 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: incremental linear transformation computation failed\n", proc );
      return( -1 );
    }
    
    /* the transformation has been calculated in real units
     */
    T->transformation_unit = REAL_UNIT;

    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    /* le champ est en voxel,
       avec le point de depart dans l'image de reference
    */
    BAL_ChangeFieldToVoxelUnit( field );

    if ( BAL_ComputeVectorFieldTransformation( T, field, estimator ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: incremental vector field transformation computation failed\n", proc );
      return( -1 );
    }

    /* the transformation has been calculated in voxel units
     */
    T->transformation_unit = VOXEL_UNIT;

    /* change into real units
     */
    if ( BAL_ChangeTransformationToRealUnit( &(T->vx), &(T->vx), T, T ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: error when translating transformation vector field from voxel to real\n", proc );
      return( -1 );
    }
    
  }

  return( 1 );

}





int BAL_ComputeTransformationResiduals( bal_transformation *T,  FIELD *field )
{ 
  char *proc = "BAL_ComputeTransformationResiduals";

  switch( T->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    /* check the units
    */
    switch ( T->transformation_unit ) {
    default :
      if ( _verbose_ ) {
        fprintf( stderr, "%s: weird, transformation type is not defined\n", proc );
        fprintf( stderr, "\t set it to REAL_UNIT, but this has to be fixed\n" );
      }
      T->transformation_unit = REAL_UNIT;
    case REAL_UNIT :
      BAL_ChangeFieldToRealUnit( field );
      break;
    case VOXEL_UNIT :
      BAL_ChangeFieldToVoxelUnit( field );
      break;
    }

    if ( BAL_LinearTrsf_Residuals( T, field ) <= 0 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: linear transformation residual computation failed\n", proc );
      return( -1 );
    }

    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    if ( T->transformation_unit != VOXEL_UNIT && field->unit != VOXEL_UNIT ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: units have to be checked\n", proc );
      return( -1 );
    }

    if ( BAL_VectorField_Residuals( T, field  ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: vector field residual transformation computation failed\n", proc );
      return( -1 );
    }
    
  }

  return( 1 );

}














/*--------------------------------------------------
 *
 * transformation use
 *
 --------------------------------------------------*/

int BAL_TransformFloatPoint( bal_floatPoint *thePt, bal_floatPoint *resPt, bal_transformation *theTr )
{
  char *proc = "BAL_TransformFloatPoint";
  double *mat;
  bal_floatPoint tmpPt;


  switch ( theTr->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    mat = theTr->mat.m;
    tmpPt.x = mat[ 0] * thePt->x + mat[ 1] * thePt->y + mat[ 2] * thePt->z + mat[ 3];
    tmpPt.y = mat[ 4] * thePt->x + mat[ 5] * thePt->y + mat[ 6] * thePt->z + mat[ 7];
    tmpPt.z = mat[ 8] * thePt->x + mat[ 9] * thePt->y + mat[10] * thePt->z + mat[11];

    resPt->x = tmpPt.x;
    resPt->y = tmpPt.y;
    resPt->z = tmpPt.z;

    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
      tmpPt = *thePt;
      resPt->x = tmpPt.x + BAL_GetXYZvalue( &(theTr->vx),
                                            tmpPt.x/theTr->vx.vx,
                                            tmpPt.y/theTr->vx.vy,
                                            tmpPt.z/theTr->vx.vz );
      resPt->y = tmpPt.y + BAL_GetXYZvalue( &(theTr->vy),
                                            tmpPt.x/theTr->vy.vx,
                                            tmpPt.y/theTr->vy.vy,
                                            tmpPt.z/theTr->vy.vz );
      if ( theTr->type == VECTORFIELD_2D )
          resPt->z = tmpPt.z;
      else {
          resPt->z = tmpPt.z + BAL_GetXYZvalue( &(theTr->vz),
                                                tmpPt.x/theTr->vz.vx,
                                                tmpPt.y/theTr->vz.vy,
                                                tmpPt.z/theTr->vz.vz );
      }
      break;


  }

  return( 1 );
}



int BAL_TransformDoublePoint( bal_doublePoint *thePt, bal_doublePoint *resPt, bal_transformation *theTr )
{
  char *proc = "BAL_TransformDoublePoint";
  double *mat;
  bal_doublePoint tmpPt;


  switch ( theTr->type ) {

  default : 
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    mat = theTr->mat.m;
    tmpPt.x = mat[ 0] * thePt->x + mat[ 1] * thePt->y + mat[ 2] * thePt->z + mat[ 3];
    tmpPt.y = mat[ 4] * thePt->x + mat[ 5] * thePt->y + mat[ 6] * thePt->z + mat[ 7];
    tmpPt.z = mat[ 8] * thePt->x + mat[ 9] * thePt->y + mat[10] * thePt->z + mat[11];

    resPt->x = tmpPt.x;
    resPt->y = tmpPt.y;
    resPt->z = tmpPt.z;

    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
      tmpPt = *thePt;
      resPt->x = tmpPt.x + BAL_GetXYZvalue( &(theTr->vx),
                                            tmpPt.x/theTr->vx.vx,
                                            tmpPt.y/theTr->vx.vy,
                                            tmpPt.z/theTr->vx.vz );
      resPt->y = tmpPt.y + BAL_GetXYZvalue( &(theTr->vy),
                                            tmpPt.x/theTr->vy.vx,
                                            tmpPt.y/theTr->vy.vy,
                                            tmpPt.z/theTr->vy.vz );
      if ( theTr->type == VECTORFIELD_2D )
          resPt->z = tmpPt.z;
      else {
          resPt->z = tmpPt.z + BAL_GetXYZvalue( &(theTr->vz),
                                                tmpPt.x/theTr->vz.vx,
                                                tmpPt.y/theTr->vz.vy,
                                                tmpPt.z/theTr->vz.vz );
      }
      break;
  }

  return( 1 );
}









/*--------------------------------------------------
 *
 * IMAGE RESAMPLING 
 *
 --------------------------------------------------*/





/*******************************************************************
 *
 * static resampling procedures
 *
 *******************************************************************/





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Reech3DCSpline4x4( bal_image *theIm, bal_image *resIm,
                                bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DCSpline4x4";
  int theDim[3];
  int resDim[3];
  int derivative[3] = {0, 0, 0};
  double *m = (double*)theTr->mat.m;

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  if ( ReechCSpline4x4( theIm->data, theIm->type, theDim,
                        resIm->data, resIm->type, resDim,
                        m, derivative ) != 1 ) {
      if ( _verbose_ )
        fprintf(stderr, "%s: some error occurs\n", proc );
      return( -1 );
  }
  return( 1 );
}





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Reech3DTriLin4x4( bal_image *theIm, bal_image *resIm,
                                bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin4x4";
  int theDim[3];
  int resDim[3];
  double *m = (double*)theTr->mat.m;

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case SSHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case USHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case FLOAT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DTriLin2DVectorField( bal_image *theIm, bal_image *resIm,
                                         bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin2DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols
       || resIm->nrows != theTr->vx.nrows
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)NULL;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech2DTriLinVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,
                                     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech2DTriLinVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech2DTriLinVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech2DTriLinVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DTriLin3DVectorField( bal_image *theIm, bal_image *resIm,
                                         bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin3DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols
       || resIm->nrows != theTr->vx.nrows
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)(theTr->vz.data);

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech3DTriLinVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,
                                     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech3DTriLinVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech3DTriLinVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech3DTriLinVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Reech3DNearest4x4( bal_image *theIm, bal_image *resIm,
                                bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest4x4";
  int theDim[3];
  int resDim[3];
  double *m = (double*)theTr->mat.m;

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case SSHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case USHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case FLOAT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DNearest2DVectorField( bal_image *theIm, bal_image *resIm,
                                         bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest2DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols
       || resIm->nrows != theTr->vx.nrows
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)NULL;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech2DNearestVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,
                                     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech2DNearestVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech2DNearestVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech2DNearestVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DNearest3DVectorField( bal_image *theIm, bal_image *resIm,
                                         bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest3DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols
       || resIm->nrows != theTr->vx.nrows
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)(theTr->vz.data);

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech3DNearestVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,
                                     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech3DNearestVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech3DNearestVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech3DNearestVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,
                                      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





/*******************************************************************
 *
 * resampling coefficient procedures
 *
 *******************************************************************/





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Coeff3DTriLin4x4( bal_image *theIm, bal_image *resIm,
                                 bal_transformation *theTr,
                                 int index )
{
  char *proc = "BAL_Coeff3DTriLin4x4";
  int theDim[3];
  int resDim[3];
  double *m = (double*)theTr->mat.m;

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  switch( resIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case FLOAT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Coeff2DTriLin4x4( theDim, (float*)(resIm->data), resDim, m, index );
    else
      Coeff3DTriLin4x4( theDim, (float*)(resIm->data), resDim, m, index );
    break;
  }
  return( 1 );
}



/*******************************************************************
 *
 * end of static resampling procedures
 *
 *******************************************************************/





/* Resample 'image' into the geometry of 'resim'
   theTr is then the transformation that goes from 'resim' to 'image'
*/

int BAL_ResampleImage( bal_image *image, bal_image *resim, 
                       bal_transformation *theTr,
                       enumTransformationInterpolation interpolation )
{
  char *proc = "BAL_ResampleImage";
  bal_transformation voxTr, tmpTr, *ptrTr;
  enumUnitTransfo unit;


  /* image to image transformation,
   * if no input transformation is given
   */
  BAL_InitTransformation( &tmpTr );
  if ( theTr == (bal_transformation *) NULL ) {
      if ( BAL_AllocTransformation( &tmpTr, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when allocating image to image transformation matrix\n", proc );
        return( -1 );
      }
      if ( BAL_ComputeImageToImageTransformation( resim, image, &tmpTr ) != 1 ) {
        if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when computing image to image transformation matrix\n", proc );
        return( -1 );
      }
      ptrTr = &tmpTr;
  }
  else {
      ptrTr = theTr;
  }


  /* test on transformation
   */
  if ( BAL_IsTransformationVectorField( ptrTr ) == 1 ) {
      if ( resim->nplanes != ptrTr->vx.nplanes
           || resim->nrows != ptrTr->vx.nrows
           || resim->ncols != ptrTr->vx.ncols ) {
        if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
        if ( _verbose_ ) {
          fprintf( stderr, "%s: result image and vector field transformation have different dimensions\n", proc );
          fprintf( stderr, "\t result image   = [%lu %lu %lu]\n", resim->ncols, resim->nrows, resim->nplanes );
          fprintf( stderr, "\t transformation = [%lu %lu %lu]\n", ptrTr->vx.ncols , ptrTr->vx.nrows, ptrTr->vx.nplanes );
        }
        return( -1 );
      }
      if ( 0.99 >= resim->vx/ptrTr->vx.vx || resim->vx/ptrTr->vx.vx >= 1.01
           || 0.99 >= resim->vy/ptrTr->vx.vy || resim->vy/ptrTr->vx.vy >= 1.01
           || 0.99 >= resim->vz/ptrTr->vx.vz || resim->vz/ptrTr->vx.vz >= 1.01 ) {
          if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
          if ( _verbose_ )
            fprintf( stderr, "%s: result image and vector field transformation have different voxel sizes\n", proc );
          return( -1 );
        }
  }

  /* test on images
   */
  if ( image->data == (void*)NULL || image->array == (void***)NULL ) {
    if ( _verbose_ )
        fprintf( stderr, "%s: input image is not allocated\n", proc );
    return( -1 );
  }
  if ( resim->data == (void*)NULL || resim->array == (void***)NULL ) {
    if ( _verbose_ )
        fprintf( stderr, "%s: result image is not allocated\n", proc );
    return( -1 );
  }


  /* passage a une transformation voxel
   * pour reechantillonnage
   */

  unit = ptrTr->transformation_unit;
  BAL_InitTransformation( &voxTr );

  switch ( ptrTr->type ) {

  default : 
    if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    
    if ( BAL_AllocTransformation( &voxTr, ptrTr->type, (bal_image *)NULL ) != 1 ) {
        if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when allocating voxel transformation matrix\n", proc );
        return( -1 );
    }
    if ( BAL_ChangeTransformationToVoxelUnit( image, resim, ptrTr, &voxTr ) != 1 ) {
        BAL_FreeTransformation( &voxTr );
        if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when translating transformation matrix from real to voxel\n", proc );
        return( -1 );
    }

    switch ( interpolation ) {
    default :
        BAL_FreeTransformation( &voxTr );
        if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
        return( -1 );
    case NEAREST :
        if ( BAL_Reech3DNearest4x4( image, resim, &voxTr ) != 1 ) {
            BAL_FreeTransformation( &voxTr );
            if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
          return( -1 );
        }
        break;
    case LINEAR :
        if ( BAL_Reech3DTriLin4x4( image, resim, &voxTr ) != 1 ) {
            BAL_FreeTransformation( &voxTr );
            if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
            BAL_FreeTransformation( &voxTr );
            return( -1 );
        }
        break;
    case CSPLINE :
        if ( BAL_Reech3DCSpline4x4( image, resim, &voxTr ) != 1 ) {
            BAL_FreeTransformation( &voxTr );
            if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
            return( -1 );
        }
        break;
    }

    BAL_FreeTransformation( &voxTr );
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    if ( resim->ncols != ptrTr->vx.ncols
       || resim->nrows != ptrTr->vx.nrows
       || resim->nplanes != ptrTr->vx.nplanes ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: output image should have the same dimensions than the transformation\n", proc );
      return( -1 );
    }

    if ( fabs( resim->vx - ptrTr->vx.vx) > ptrTr->vx.vx / 1000.0
         || fabs( resim->vy - ptrTr->vx.vy) > ptrTr->vx.vy / 1000.0
         || fabs( resim->vz - ptrTr->vx.vz) > ptrTr->vx.vz / 1000.0 ) {
           if ( _verbose_ )
        fprintf( stderr, "%s: output image should have the same voxel sizes than the transformation\n", proc );
      return( -1 );
    }

    if ( unit != VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToVoxelUnit( image, resim, ptrTr, ptrTr ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when translating transformation vector field from real to voxel\n", proc );
        return( -1 );
      }
    }

    switch ( ptrTr->type ) {
    default :
      break;
    case VECTORFIELD_2D :
      switch ( interpolation ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
        return( -1 );
      case NEAREST :
        if ( BAL_Reech3DNearest2DVectorField( image, resim, ptrTr ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
          return( -1 );
        }
        break;
      case LINEAR :
        if ( BAL_Reech3DTriLin2DVectorField( image, resim, ptrTr ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
          return( -1 );
        }
        break;
      }
      break;
    case VECTORFIELD_3D :
      switch ( interpolation ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
        return( -1 );
      case NEAREST :
        if ( BAL_Reech3DNearest3DVectorField( image, resim, ptrTr ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
          return( -1 );
        }
        break;
      case LINEAR :
        if ( BAL_Reech3DTriLin3DVectorField( image, resim, ptrTr ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
          return( -1 );
        }
        break;
      }
      break;
    }

    if ( unit != VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToRealUnit( image, resim, ptrTr, ptrTr ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when translating transformation vector field from voxel to real\n", proc );
        return( -1 );
      }
    }
    
    break;

  }

  if ( theTr == (bal_transformation *) NULL ) BAL_FreeTransformation( &tmpTr );

  return( 1 );
}




int BAL_LinearResamplingCoefficients( bal_image *image, bal_image *resim,
                                      bal_transformation *theTr,
                                      enumTransformationInterpolation interpolation,
                                      int index )
{
  char *proc = "BAL_LinearResamplingCoefficients";
  bal_transformation voxTr;
  bal_transformation *ptrTr = theTr;

  BAL_InitTransformation( &voxTr );

  switch ( ptrTr->type ) {

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    if ( BAL_AllocTransformation( &voxTr, ptrTr->type, (bal_image *)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when allocating voxel transformation matrix\n", proc );
        return( -1 );
    }
    if ( BAL_ChangeTransformationToVoxelUnit( image, resim, ptrTr, &voxTr ) != 1 ) {
        BAL_FreeTransformation( &voxTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when translating transformation matrix from real to voxel\n", proc );
        return( -1 );
    }

    switch ( interpolation ) {
    default :
        BAL_FreeTransformation( &voxTr );
        if ( _verbose_ )
          fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
        return( -1 );
    case LINEAR :
        if ( BAL_Coeff3DTriLin4x4( image, resim, &voxTr, index ) != 1 ) {
            BAL_FreeTransformation( &voxTr );
            if ( _verbose_ )
              fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
            BAL_FreeTransformation( &voxTr );
            return( -1 );
        }
        break;
    }

    BAL_FreeTransformation( &voxTr );
    break;

  }

  return( 1 );
}










/*--------------------------------------------------
 *
 * Non-Linear Transformations Construction
 *
 --------------------------------------------------*/



static double sinusoid_amplitude[3] = {2.0, 2.0, 2.0};
static double sinusoid_period[3] = {50.0, 50.0, 50.0};

void BAL_SetSinusoidAmplitudeForVectorFieldTransformation( double *a )
{
    int i;
    for ( i=0; i<3; i++ )
        sinusoid_amplitude[i] = a[i];
}

void BAL_SetSinusoidPeriodForVectorFieldTransformation( double *p )
{
    int i;
    for ( i=0; i<3; i++ )
        sinusoid_period[i] = p[i];
}

static int BAL_SinusoidVectorField( bal_transformation *theTrsf, int drawx, int drawy, int drawz )
{
  char *proc = "BAL_SinusoidVectorField";
  size_t i, x, y, z;
  float *xbuf, *ybuf, *zbuf;
  double wx, wy, wz;

  if ( BAL_IsTransformationVectorField( theTrsf ) != 1 ) {
      if ( _verbose_ ) {
          fprintf( stderr, "%s: input transformations are not vector field\n", proc );
      }
      return( -1 );
  }

  if ( drawz && drawy && drawx ) {
    xbuf = (float*)(theTrsf->vx.data);
    ybuf = (float*)(theTrsf->vy.data);
    zbuf = (float*)(theTrsf->vz.data);
    wx = 2 * 3.14159265 / sinusoid_period[0];
    wy = 2 * 3.14159265 / sinusoid_period[1];
    wz = 2 * 3.14159265 / sinusoid_period[2];
    for ( i=0, z=0; z<theTrsf->vx.nplanes; z++ )
    for ( y=0; y<theTrsf->vx.nrows; y++ )
    for ( x=0; x<theTrsf->vx.ncols; x++, i++ )
      xbuf[i] = ybuf[i] = zbuf[i] = sinusoid_amplitude[2] * sin( wz * (z-theTrsf->vz.nplanes/2) )
        * sinusoid_amplitude[1] * sin( wy * (y-theTrsf->vy.nrows/2) )
        * sinusoid_amplitude[0] * sin( wx * (x-theTrsf->vx.ncols/2) );
  }

  if ( !drawz && drawy && drawx ) {
    xbuf = (float*)(theTrsf->vx.data);
    ybuf = (float*)(theTrsf->vy.data);
    wx = 2 * 3.14159265 / sinusoid_period[0];
    wy = 2 * 3.14159265 / sinusoid_period[1];
    for ( i=0, z=0; z<theTrsf->vx.nplanes; z++ )
    for ( y=0; y<theTrsf->vx.nrows; y++ )
    for ( x=0; x<theTrsf->vx.ncols; x++, i++ )
      xbuf[i] = ybuf[i] = sinusoid_amplitude[1] * sin( wx * (y-theTrsf->vy.nrows/2) )
        * sinusoid_amplitude[0] * sin( wy * (x-theTrsf->vx.ncols/2) );
    if ( theTrsf->vz.data != (void*)NULL
         && theTrsf->vz.nplanes == theTrsf->vx.nplanes
         && theTrsf->vz.nrows == theTrsf->vx.nrows
         && theTrsf->vz.ncols == theTrsf->vx.ncols ) {
        zbuf = (float*)(theTrsf->vz.data);
        for ( i=0, z=0; z<theTrsf->vx.nplanes; z++ )
        for ( y=0; y<theTrsf->vx.nrows; y++ )
        for ( x=0; x<theTrsf->vx.ncols; x++, i++ )
          zbuf[i] = 0.0;
    }
  }


  return( 1 );

}


int BAL_Sinusoid3DVectorField( bal_transformation *theTrsf )
{
  if ( 0 ) {
    sinusoid_amplitude[0] = 1.3;
    sinusoid_amplitude[1] = 1.3;
    sinusoid_amplitude[2] = 1.3;
  }
  return( BAL_SinusoidVectorField( theTrsf, 1, 1, 1 ) );
}

int BAL_Sinusoid2DVectorField( bal_transformation *theTrsf )
{
  if ( 0 ) {
    sinusoid_amplitude[0] = 1.75;
    sinusoid_amplitude[1] = 1.75;
    sinusoid_amplitude[2] = 1.75;
  }
  return( BAL_SinusoidVectorField( theTrsf, 1, 1, 0 ) );
}





/*--------------------------------------------------
 *
 * Random Transformations
 *
 --------------------------------------------------*/

static double angle_interval[2] = {0.0, 1.57079}; /* 0, pi/2 */
static double scale_interval[2] = {0.7, 1.4};
static double shear_interval[2] = {0.0, 0.3};
static double translation_interval[2] = {-10, 10};

void BAL_SetMinAngleForRandomTransformation( double d )
{
    angle_interval[0] = d;
}

void BAL_SetMaxAngleForRandomTransformation( double d )
{
    angle_interval[1] = d;
}

void BAL_SetMinScaleForRandomTransformation( double d )
{
    scale_interval[0] = d;
}

void BAL_SetMaxScaleForRandomTransformation( double d )
{
    scale_interval[1] = d;
}

void BAL_SetMinShearForRandomTransformation( double d )
{
    shear_interval[0] = d;
}

void BAL_SetMaxShearForRandomTransformation( double d )
{
    shear_interval[1] = d;
}

void BAL_SetMinTranslationForRandomTransformation( double d )
{
    translation_interval[0] = d;
}

void BAL_SetMaxTranslationForRandomTransformation( double d )
{
    translation_interval[1] = d;
}





int BAL_Random2DTranslationMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _IdentityMatrix( m );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DTranslationMatrix( bal_transformation *theTrsf )
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _IdentityMatrix( m );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DTranslationScalingMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DScaleMatrix( m, scale_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DTranslationScalingMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DScaleMatrix( m, scale_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DRigidMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DRotationMatrix( m, angle_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DRigidMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DRotationMatrix( m, angle_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DSimilitudeMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DSimilitudeMatrix( m, angle_interval, scale_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DSimilitudeMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DSimilitudeMatrix( m, angle_interval, scale_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DAffineMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DAffineMatrix( m, angle_interval, scale_interval, shear_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DAffineMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DAffineMatrix( m, angle_interval, scale_interval, shear_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DVectorField( bal_transformation *theTrsf __attribute__ ((unused)) )
{
    char *proc = "BAL_Random2DVectorField";
    if ( _verbose_ ) {
        fprintf( stderr, "%s: not implemented yet\n", proc );
    }
    return( -1 );
}



int BAL_Random3DVectorField( bal_transformation *theTrsf __attribute__ ((unused)) )
{
    char *proc = "BAL_Random3DVectorField";
    if ( _verbose_ ) {
        fprintf( stderr, "%s: not implemented yet\n", proc );
    }
    return( -1 );
}



int BAL_SetTransformationToRandom( bal_transformation *theTrsf )
{
  char *proc = "BAL_SetTransformationToRandom";
  
  if ( _verbose_ >= 3 ) {
      fprintf( stderr, "%s: compute random transformation of type ", proc );
      BAL_PrintTransformationType( stderr, theTrsf->type );
  }

  switch( theTrsf->type ) {
    
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );

  case VECTORFIELD_2D :
    if (  BAL_Random2DVectorField( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D vector field\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;
    
  case VECTORFIELD_3D :
    if (  BAL_Random3DVectorField( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D vector field\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;
    
  case TRANSLATION_2D :
    if ( BAL_Random2DTranslationMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D translation matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_3D :
    if ( BAL_Random3DTranslationMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D translation matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_SCALING_2D :
    if ( BAL_Random2DTranslationScalingMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D translation and scaling matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_SCALING_3D :
    if ( BAL_Random3DTranslationScalingMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D translation and scaling matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case RIGID_2D :
    if ( BAL_Random2DRigidMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D rigid matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case RIGID_3D :
    if ( BAL_Random3DRigidMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D rigid matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case SIMILITUDE_2D :
    if ( BAL_Random2DSimilitudeMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D similitude matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case SIMILITUDE_3D :
    if ( BAL_Random3DSimilitudeMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D similitude matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case AFFINE_2D :
    if ( BAL_Random2DAffineMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 2D affine matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case AFFINE_3D :
    if ( BAL_Random3DAffineMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to generate 3D affine matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  }

  return( 1 );
}
