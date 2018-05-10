/*************************************************************************
 * bal-transformation-copy.c -
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

#include <stdlib.h>
#include <stdio.h>

#include <chunks.h>

#include <bal-transformation-copy.h>





static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationCopy( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationCopy( )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationCopy( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationCopy( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalTransformationCopy( )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationCopy( )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}










/************************************************************
 *
 * copy transformations
 *
 ************************************************************/



static int _Resample2DVectorfield( bal_image *theVx,
                                   bal_image *theVy,
                                   bal_image *resVx,
                                   bal_image *resVy )
{
  char *proc = "_Resample2DVectorfield";
  unsigned int i, j, k, n;
  double *res_to_r = resVx->to_real.m;
  double *the_to_z = theVx->to_voxel.m;
  double x, y, z, u, v, w;

  if ( resVx->type != resVy->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  switch( resVx->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );

  case FLOAT :
    {
      float *bufX = (float*)resVx->data;
      float *bufY = (float*)resVy->data;
      for ( n=0, k=0; k<resVx->nplanes; k++ )
      for ( j=0; j<resVx->nrows; j++ )
      for ( i=0; i<resVx->ncols; i++, n++ ) {

        switch( resVx->geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          x = res_to_r[ 0] * i;
          y =                    res_to_r[ 5] * j;
          z =                                       res_to_r[10] * k;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          x = res_to_r[ 0] * i                                       + res_to_r[ 3];
          y =                    res_to_r[ 5] * j                    + res_to_r[ 7];
          z =                                       res_to_r[10] * k + res_to_r[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          x = res_to_r[ 0] * i + res_to_r[ 1] * j + res_to_r[ 2] * k + res_to_r[ 3];
          y = res_to_r[ 4] * i + res_to_r[ 5] * j + res_to_r[ 6] * k + res_to_r[ 7];
          z = res_to_r[ 8] * i + res_to_r[ 9] * j + res_to_r[10] * k + res_to_r[11];
          break;
        }

        switch( theVx->geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          u = the_to_z[ 0] * x;
          v =                    the_to_z[ 5] * y;
          w =                                       the_to_z[10] * z;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          u = the_to_z[ 0] * x                                       + the_to_z[ 3];
          v =                    the_to_z[ 5] * y                    + the_to_z[ 7];
          w =                                       the_to_z[10] * z + the_to_z[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          u = the_to_z[ 0] * x + the_to_z[ 1] * y + the_to_z[ 2] * z + the_to_z[ 3];
          v = the_to_z[ 4] * x + the_to_z[ 5] * y + the_to_z[ 6] * z + the_to_z[ 7];
          w = the_to_z[ 8] * x + the_to_z[ 9] * y + the_to_z[10] * z + the_to_z[11];
          break;
        }

        bufX[n] = BAL_GetXYKvalue( theVx, u, v, (int)(w+0.5) );
        bufY[n] = BAL_GetXYKvalue( theVy, u, v, (int)(w+0.5) );
      }
    }
    break;
  }
  return( 1 );
}





static int _Resample3DVectorfield( bal_image *theVx,
                                   bal_image *theVy,
                                   bal_image *theVz,
                                   bal_image *resVx,
                                   bal_image *resVy,
                                   bal_image *resVz )
{
  char *proc = "_Resample3DVectorfield";
  unsigned int i, j, k, n;
  double *res_to_r = resVx->to_real.m;
  double *the_to_z = theVx->to_voxel.m;
  double x, y, z, u, v, w;

  if ( resVx->type != resVy->type || resVx->type != resVz->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  switch( resVx->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );

  case FLOAT :
    {
      float *bufX = (float*)resVx->data;
      float *bufY = (float*)resVy->data;
      float *bufZ = (float*)resVz->data;
      for ( n=0, k=0; k<resVx->nplanes; k++ )
      for ( j=0; j<resVx->nrows; j++ )
      for ( i=0; i<resVx->ncols; i++, n++ ) {

        switch( resVx->geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          x = res_to_r[ 0] * i;
          y =                    res_to_r[ 5] * j;
          z =                                       res_to_r[10] * k;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          x = res_to_r[ 0] * i                                       + res_to_r[ 3];
          y =                    res_to_r[ 5] * j                    + res_to_r[ 7];
          z =                                       res_to_r[10] * k + res_to_r[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          x = res_to_r[ 0] * i + res_to_r[ 1] * j + res_to_r[ 2] * k + res_to_r[ 3];
          y = res_to_r[ 4] * i + res_to_r[ 5] * j + res_to_r[ 6] * k + res_to_r[ 7];
          z = res_to_r[ 8] * i + res_to_r[ 9] * j + res_to_r[10] * k + res_to_r[11];
          break;
        }

        switch( theVx->geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          u = the_to_z[ 0] * x;
          v =                    the_to_z[ 5] * y;
          w =                                       the_to_z[10] * z;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          u = the_to_z[ 0] * x                                       + the_to_z[ 3];
          v =                    the_to_z[ 5] * y                    + the_to_z[ 7];
          w =                                       the_to_z[10] * z + the_to_z[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          u = the_to_z[ 0] * x + the_to_z[ 1] * y + the_to_z[ 2] * z + the_to_z[ 3];
          v = the_to_z[ 4] * x + the_to_z[ 5] * y + the_to_z[ 6] * z + the_to_z[ 7];
          w = the_to_z[ 8] * x + the_to_z[ 9] * y + the_to_z[10] * z + the_to_z[11];
          break;
        }

        bufX[n] = BAL_GetXYZvalue( theVx, u, v, w );
        bufY[n] = BAL_GetXYZvalue( theVy, u, v, w );
        bufZ[n] = BAL_GetXYZvalue( theVz, u, v, w );
      }
    }
    break;
  }
  return( 1 );
}










/* made compliant to QForm matrix
 */
int BAL_CopyTransformation( bal_transformation *theTrsf,
                            bal_transformation *resTrsf )
{
  char *proc="BAL_CopyTransformation";
  size_t i, j, k, v;
  double *m = (double*)NULL;
  float ***vecx, ***vecy, ***vecz;
  bal_transformation tmpTrsf;

  /* nothing to do
   */
  if ( theTrsf == resTrsf ) return( 1 );


  switch ( theTrsf->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet for input transformation\n", proc );
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

    switch( resTrsf->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such type not handled yet for output transformation (matrix case)\n", proc );
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

      v = theTrsf->mat.c * theTrsf->mat.l;
      for ( i=0; i<v; i++ ) resTrsf->mat.m[i] = theTrsf->mat.m[i];
      resTrsf->type = theTrsf->type;
      resTrsf->transformation_unit = theTrsf->transformation_unit;

      break;

    case VECTORFIELD_2D :

      BAL_InitTransformation( &tmpTrsf );
      if ( BAL_AllocTransformation( &tmpTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
        return( -1 );
      }
      v = theTrsf->mat.c * theTrsf->mat.l;
      for ( i=0; i<v; i++ ) tmpTrsf.mat.m[i] = theTrsf->mat.m[i];
      tmpTrsf.type = theTrsf->type;
      tmpTrsf.transformation_unit = theTrsf->transformation_unit;

      switch ( theTrsf->transformation_unit ) {
      default :
        BAL_FreeTransformation( &tmpTrsf );
        if ( _verbose_ )
          fprintf( stderr, "%s: such transformation unit no handled yet (matrix -> 2D vector field)\n", proc );
        return( -1 );

      case VOXEL_UNIT :
        BAL_FreeTransformation( &tmpTrsf );
        if ( _verbose_ )
          fprintf( stderr, "%s: such case no handled yet (matrix (voxel unit) -> 2D vector field)\n", proc );
        return( -1 );

      case REAL_UNIT :
        /* v_R( M_{ref,R} ) = ( T_{flo<-ref,R} - Id ) o H_{ref,R<-Z} M_{ref,Z}
         *  0  1  2  3
         *  4  5  6  7
         *  8  9 10 11
         * 12 13 14 15
         */
        tmpTrsf.mat.m[ 0] -= 1.0;
        tmpTrsf.mat.m[ 5] -= 1.0;
        tmpTrsf.mat.m[10] -= 1.0;
        /* tmpTrsf o resTrsf->vx.to_real = tmpTrsf
         */
        _mult_mat( &(tmpTrsf.mat), &(resTrsf->vx.to_real), &(tmpTrsf.mat) );

        m = tmpTrsf.mat.m;
        vecx = (float***)resTrsf->vx.array;
        vecy = (float***)resTrsf->vy.array;

        for ( k=0; k<resTrsf->vx.nplanes; k++ )
        for ( j=0; j<resTrsf->vx.nrows; j++ )
        for ( i=0; i<resTrsf->vx.ncols; i++  ) {
          vecx[k][j][i] = m[ 0] * i + m[ 1] * j + m[ 2] * k + m[ 3];
          vecy[k][j][i] = m[ 4] * i + m[ 5] * j + m[ 6] * k + m[ 7];
        }
        resTrsf->transformation_unit = REAL_UNIT;
        break;
      }

      BAL_FreeTransformation( &tmpTrsf );

      break;

    case VECTORFIELD_3D :

      BAL_InitTransformation( &tmpTrsf );
      if ( BAL_AllocTransformation( &tmpTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
        return( -1 );
      }
      v = theTrsf->mat.c * theTrsf->mat.l;
      for ( i=0; i<v; i++ ) tmpTrsf.mat.m[i] = theTrsf->mat.m[i];
      tmpTrsf.type = theTrsf->type;
      tmpTrsf.transformation_unit = theTrsf->transformation_unit;

      switch ( theTrsf->transformation_unit ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such transformation unit no handled yet (matrix -> 3D vector case)\n", proc );
        return( -1 );

      case VOXEL_UNIT :
        BAL_FreeTransformation( &tmpTrsf );
        if ( _verbose_ )
          fprintf( stderr, "%s: such case no handled yet (matrix (voxel unit) -> 2D vector field)\n", proc );
        return( -1 );

      case REAL_UNIT :
        /* v_R( M_{ref,R} ) = ( T_{flo<-ref,R} - Id ) o H_{ref,R<-Z} M_{ref,Z}
         *  0  1  2  3
         *  4  5  6  7
         *  8  9 10 11
         * 12 13 14 15
         */
        tmpTrsf.mat.m[ 0] -= 1.0;
        tmpTrsf.mat.m[ 5] -= 1.0;
        tmpTrsf.mat.m[10] -= 1.0;
        /* tmpTrsf o resTrsf->vx.to_real = tmpTrsf
         */
        _mult_mat( &(tmpTrsf.mat), &(resTrsf->vx.to_real), &(tmpTrsf.mat) );

        m = tmpTrsf.mat.m;
        vecx = (float***)resTrsf->vx.array;
        vecy = (float***)resTrsf->vy.array;
        vecz = (float***)resTrsf->vz.array;

        for ( k=0; k<resTrsf->vx.nplanes; k++ )
        for ( j=0; j<resTrsf->vx.nrows; j++ )
        for ( i=0; i<resTrsf->vx.ncols; i++ ) {
          vecx[k][j][i] = m[ 0] * i + m[ 1] * j + m[ 2] * k + m[ 3];
          vecy[k][j][i] = m[ 4] * i + m[ 5] * j + m[ 6] * k + m[ 7];
          vecz[k][j][i] = m[ 8] * i + m[ 9] * j + m[10] * k + m[11];
        }
        resTrsf->transformation_unit = REAL_UNIT;
        break;
      }

      BAL_FreeTransformation( &tmpTrsf );

      break;

    }

    break;
    /* end of matrix case for theTrsf
     */

  case VECTORFIELD_2D :

    switch( resTrsf->type ) {
    default :
      if ( _verbose_ ) {
        fprintf( stderr, "%s: such type not handled yet for output transformation (2D vector field)\n", proc );
      }
      return( -1 );

    case VECTORFIELD_3D :

      if( BAL_FillImage( &(resTrsf->vz), 0.0 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to set Z component of vector field to 0\n", proc );
        return( -1 );
      }

    case VECTORFIELD_2D :

      if ( theTrsf->vx.ncols == resTrsf->vx.ncols
           && theTrsf->vx.nrows == resTrsf->vx.nrows
           && theTrsf->vx.nplanes == resTrsf->vx.nplanes ) {

        if ( BAL_CopyImage( &(theTrsf->vx), &(resTrsf->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy X component of vector field\n", proc );
          return( -1 );
        }
        if ( BAL_CopyImage( &(theTrsf->vy), &(resTrsf->vy) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy Y component of vector field\n", proc );
          return( -1 );
        }
        resTrsf->transformation_unit = theTrsf->transformation_unit;

      }
      else {
        switch( theTrsf->transformation_unit ) {
        default :
        case VOXEL_UNIT :
          if ( _verbose_ )
            fprintf( stderr, "%s: such unit not handled yet\n", proc );
          return( -1 );
        case REAL_UNIT :
            if ( _Resample2DVectorfield( &(theTrsf->vx), &(theTrsf->vy),
                                         &(resTrsf->vx), &(resTrsf->vy) ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to resample 2D vector field\n", proc );
            return( -1 );
          }
          resTrsf->transformation_unit = theTrsf->transformation_unit;
        }
      }

    }

    break;
    /* end of VECTORFIELD_2D case for theTrsf
     */

  case VECTORFIELD_3D :

    switch( resTrsf->type ) {
    default :
      if ( _verbose_ ) {
        fprintf( stderr, "%s: such type not handled yet for output transformation (3D vector field)\n", proc );
      }
      return( -1 );

    case VECTORFIELD_3D :

      if ( theTrsf->vx.ncols == resTrsf->vx.ncols
           && theTrsf->vx.nrows == resTrsf->vx.nrows
           && theTrsf->vx.nplanes == resTrsf->vx.nplanes ) {

        if ( BAL_CopyImage( &(theTrsf->vx), &(resTrsf->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy X component of vector field\n", proc );
          return( -1 );
        }
        if ( BAL_CopyImage( &(theTrsf->vy), &(resTrsf->vy) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy Y component of vector field\n", proc );
          return( -1 );
        }
        if ( BAL_CopyImage( &(theTrsf->vz), &(resTrsf->vz) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy Z component of vector field\n", proc );
          return( -1 );
        }
        resTrsf->transformation_unit = theTrsf->transformation_unit;
      }
      else {
        switch( theTrsf->transformation_unit ) {
        default :
        case VOXEL_UNIT :
          if ( _verbose_ )
            fprintf( stderr, "%s: such unit not handled yet\n", proc );
          return( -1 );
        case REAL_UNIT :
          if ( _Resample3DVectorfield( &(theTrsf->vx), &(theTrsf->vy),
                                       &(theTrsf->vz), &(resTrsf->vx),
                                       &(resTrsf->vy), &(resTrsf->vz) ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to resample 3D vector field\n", proc );
            return( -1 );
          }
          resTrsf->transformation_unit = theTrsf->transformation_unit;
        }
      }

    }

    break;
    /* end of VECTORFIELD_3D case for theTrsf
     */
  }

  return( 1 );
}










int BAL_CopyTransformationList( bal_transformationList *theTrsfs,
                                bal_transformationList *resTrsfs )
{
  char *proc = "BAL_CopyTransformationList";
  int i;

  /* nothing to do
   */
  if ( theTrsfs == resTrsfs ) return( 1 );

  if ( theTrsfs->n_trsfs > resTrsfs->n_allocated_trsfs ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: not enough allocated transformations in the destination list\n", proc );
    return( -1 );
  }

  for ( i=0; i<theTrsfs->n_trsfs; i++ ) {
    if ( BAL_DoesTransformationExist( &(theTrsfs->data[i]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: input transformation #%d has not been allocated or filled\n", proc, i );
      return( -1 );
    }
    if ( BAL_DoesTransformationExist( &(resTrsfs->data[i]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: result transformation #%d has not been allocated or filled\n", proc, i );
      return( -1 );
    }
    if ( BAL_CopyTransformation( &(theTrsfs->data[i]), &(resTrsfs->data[i]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when copying transformation #%d\n", proc, i );
      return( -1 );
    }
  }

  resTrsfs->n_trsfs = theTrsfs->n_trsfs;

  return( 1 );
}








/************************************************************
 *
 *
 *
 ************************************************************/


/*--------------------------------------------------
 *
 * transformation conversion
 *
 --------------------------------------------------*/



typedef struct _UnitTransformationChangeParam {
  bal_image *floIm;
  bal_image *refIm;
  bal_transformation *theTrsf;
  bal_transformation *resTrsf;
  bal_transformation *tmpTrsf;
  bal_imageGeometry tmpGeometry;
} _UnitTransformationChangeParam;





static bal_imageGeometry _combinedGeometry( bal_imageGeometry g1, bal_imageGeometry g2 )
{
  switch( g1 ) {
  default :
  case _BAL_UNKNOWN_GEOMETRY_ :
  case _BAL_HOMOTHETY_GEOMETRY_ :
    switch( g2 ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      return( _BAL_HOMOTHETY_GEOMETRY_ );
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      return( _BAL_TRANSLATION_GEOMETRY_ );
      break;
    case _BAL_QFORM_GEOMETRY_ :
      return( _BAL_QFORM_GEOMETRY_ );
      break;
    }
    break;
  case _BAL_TRANSLATION_GEOMETRY_ :
    switch( g2 ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
    case _BAL_TRANSLATION_GEOMETRY_ :
      return( _BAL_TRANSLATION_GEOMETRY_ );
      break;
    case _BAL_QFORM_GEOMETRY_ :
      return( _BAL_QFORM_GEOMETRY_ );
      break;
    }
    break;
  case _BAL_QFORM_GEOMETRY_ :
    return( _BAL_QFORM_GEOMETRY_ );
    break;
  }

  /* default
   * should not reach this point
   */
  return( _BAL_QFORM_GEOMETRY_ );
}





/*
   theTr is the transformation that goes from 'resim' to 'image'

   - linear case:
     we have M_floIm = T * M_ref, with M_floIm and M_ref in voxels
     we go into real coordinates with
     H_flo * M_flo = H_flo * T * M_ref
                       = H_flo * T * H^(-1)_ref * H_ref * M_ref
     thus H_flo * T * H^(-1)_ref is the transformation in real coordinates

   - vector field case:
     we have M_flo = M_ref + V(M_ref)
     we go into real coordinates with
     H_flo * M_flo = H_flo * M_ref + H_flo * V(M_ref)
                       = H_ref * M_ref + (H_flo - H_ref) * M_ref + H_flo * V(M_ref)
     thus, the displacement vector in real coordinates is
     (H_flo - H_ref) * M_ref + H_flo * V(M_ref)
*/



static void *_Change2DVectorFieldToRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *floIm = p->floIm;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float x, y;
  double *h = floIm->to_real.m;
  double *m = p->tmpTrsf->mat.m;

  size_t ncols = (theTrsf->vx).ncols; /* dimx */
  size_t nrows = (theTrsf->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;


  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  /*  0  1  2  3
   *  4  5  6  7
   *  8  9 10 11
   * 12 13 14 15
   */
  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    switch( floIm->geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
    case _BAL_TRANSLATION_GEOMETRY_ :
      x = h[ 0] * theVecX[k][j][i]                           ;
      y =                            h[ 5] * theVecY[k][j][i];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      x = h[ 0] * theVecX[k][j][i] + h[ 1] * theVecY[k][j][i];
      y = h[ 4] * theVecX[k][j][i] + h[ 5] * theVecY[k][j][i];
      break;
    }
    switch( p->tmpGeometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i                      + x;
      resVecY[k][j][i] =             m [ 5] * j         + y;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i              + m[ 3] + x;
      resVecY[k][j][i] =             m [ 5] * j + m[ 7] + y;
      break;
    case _BAL_QFORM_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i + m [ 1] * j + m[ 3] + x;
      resVecY[k][j][i] = m[ 4] * i + m [ 5] * j + m[ 7] + y;
      break;
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_Change3DVectorFieldToRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *floIm = p->floIm;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***theVecZ = (float***)(theTrsf->vz.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float ***resVecZ = (float***)(resTrsf->vz.array);
  float x, y, z;
  double *h = floIm->to_real.m;
  double *m = p->tmpTrsf->mat.m;

  size_t ncols = (theTrsf->vx).ncols; /* dimx */
  size_t nrows = (theTrsf->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;


  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  /*  0  1  2  3
   *  4  5  6  7
   *  8  9 10 11
   * 12 13 14 15
   */
  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    switch( floIm->geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
    case _BAL_TRANSLATION_GEOMETRY_ :
      x = h[ 0] * theVecX[k][j][i]                                                      ;
      y =                            h[ 5] * theVecY[k][j][i]                           ;
      z =                                                       h[10] * theVecZ[k][j][i];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      x = h[ 0] * theVecX[k][j][i] + h[ 1] * theVecY[k][j][i] + h[ 2] * theVecZ[k][j][i];
      y = h[ 4] * theVecX[k][j][i] + h[ 5] * theVecY[k][j][i] + h[ 6] * theVecZ[k][j][i];
      z = h[ 8] * theVecX[k][j][i] + h[ 9] * theVecY[k][j][i] + h[10] * theVecZ[k][j][i];
      break;
    }
    switch( p->tmpGeometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i                                   + x;
      resVecY[k][j][i] =             m [ 5] * j                      + y;
      resVecZ[k][j][i] =                          m [10] * k         + z;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i                           + m[ 3] + x;
      resVecY[k][j][i] =             m [ 5] * j              + m[ 7] + y;
      resVecZ[k][j][i] =                          m [10] * k + m[11] + z;
      break;
    case _BAL_QFORM_GEOMETRY_ :
      resVecX[k][j][i] = m[ 0] * i + m [ 1] * j + m [ 2] * k + m[ 3] + x;
      resVecY[k][j][i] = m[ 4] * i + m [ 5] * j + m [ 6] * k + m[ 7] + y;
      resVecZ[k][j][i] = m[ 8] * i + m [ 9] * j + m [10] * k + m[11] + z;
      break;
    }
  }
  chunk->ret = 1;
  return( (void*)NULL );
}






/* theTrsf allows to resample 'floIm' into 'refIm' thus goes
 * from 'refIm' to 'floIm'
*/
int BAL_ChangeTransformationToRealUnit( bal_image *floIm,
                                        bal_image *refIm,
                                        bal_transformation *theTrsf,
                                        bal_transformation *resTrsf )
{
  char *proc = "BAL_ChangeTransformationToRealUnit";
  bal_transformation tmpTrsf;

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _UnitTransformationChangeParam p;

  if ( floIm == (bal_image*)NULL || refIm == (bal_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL template images\n", proc );
    return( -1 );
  }

  if ( theTrsf->type != resTrsf->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations have different types\n", proc );
    return( -1 );
  }

  /* nothing to do, but copy the transformation
   */
  if ( theTrsf->transformation_unit == REAL_UNIT ) {
    if ( BAL_CopyTransformation( theTrsf, resTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy transformation\n", proc );
      return( -1 );
    }
    return( 1 );
  }

  /* something to do
   */

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
    return( -1 );
  }

  /* preparing parallelism
   */
  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    first = 0;
    last = theTrsf->vx.nplanes * theTrsf->vx.nrows * theTrsf->vx.ncols - 1;
    initChunks( &chunks );
    if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
    }

    /* floIm->to_real - refIm->to_real = tmpTrsf
     */
    _sub_mat( &(floIm->to_real), &(refIm->to_real), &(tmpTrsf.mat) );

    p.floIm = floIm;
    p.refIm = refIm;
    p.theTrsf = theTrsf;
    p.resTrsf = resTrsf;
    p.tmpTrsf = &tmpTrsf;
    p.tmpGeometry = _combinedGeometry( floIm->geometry, refIm->geometry );

    for ( i=0; i<chunks.n_allocated_chunks; i++ )
      chunks.data[i].parameters = (void*)(&p);
    break;
  }



  switch ( theTrsf->type ) {
  default :
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

    /* theTrsf o refIm->to_voxel = tmpTrsf
     */
    _mult_mat( &(theTrsf->mat), &(refIm->to_voxel), &(tmpTrsf.mat) );

    /* floIm->to_real o tmpTrsf = resTrsf
     */
    _mult_mat( &(floIm->to_real), &(tmpTrsf.mat), &(resTrsf->mat) );

    break;

  case VECTORFIELD_2D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols
         || resTrsf->vx.nrows != theTrsf->vx.nrows
         || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
                 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
                 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      freeChunks( &chunks );
      return( -1 );
    }

    if ( processChunks( &_Change2DVectorFieldToRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to change 2D vector field into real unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case VECTORFIELD_3D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols
         || resTrsf->vx.nrows != theTrsf->vx.nrows
         || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
                 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
                 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      freeChunks( &chunks );
      return( -1 );
    }

    if ( processChunks( &_Change3DVectorFieldToRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to change 3D vector field into real unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  BAL_FreeTransformation( &tmpTrsf );

  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    freeChunks( &chunks );
    break;
  }

  resTrsf->transformation_unit = REAL_UNIT;

  return( 1 );
}






/*
   theTr is the transformation that goes from 'resim' to 'image'

   - linear case:
     we have H_flo * M_flo = T * H_ref * M_ref,
        with M_flo and M_ref in voxel coordinates
     we go into voxel coordinates with
     M_flo = H^(-1)_flo * T * H_ref * M_ref,
     thus H^(-1)_flo * T * H_ref is the transformation in voxel coordinates

   - vector field case:
     we have H_flo * M_flo = H_ref * M_ref + V(M_ref)
     we go into voxel coordinates with
     M_flo = H^(-1)_flo * H_ref * M_ref + H^(-1)_flo * V(M_ref)
             = M_ref + ( H^(-1)_flo * H_ref - Id ) * M_ref + H^(-1)_flo * V(M_ref)
     thus, the displacement vector in real coordinates is
     ( H^(-1)_flo * H_ref - Id ) * M_ref + H^(-1)_flo * V(M_ref)
*/



static void *_Change2DVectorFieldToVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *floIm = p->floIm;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float x, y;
  double *h = floIm->to_voxel.m;
  double *m = p->tmpTrsf->mat.m;

  size_t ncols = (theTrsf->vx).ncols; /* dimx */
  size_t nrows = (theTrsf->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;


  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  /*  0  1  2  3
   *  4  5  6  7
   *  8  9 10 11
   * 12 13 14 15
   */
  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    x = h[ 0] * theVecX[k][j][i] + h[ 1] * theVecY[k][j][i];
    y = h[ 4] * theVecX[k][j][i] + h[ 5] * theVecY[k][j][i];
    resVecX[k][j][i] = m[ 0] * i + m [ 1] * j + m[ 3] + x;
    resVecY[k][j][i] = m[ 4] * i + m [ 5] * j + m[ 7] + y;
  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_Change3DVectorFieldToVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *floIm = p->floIm;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***theVecZ = (float***)(theTrsf->vz.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float ***resVecZ = (float***)(resTrsf->vz.array);
  float x, y, z;
  double *h = floIm->to_voxel.m;
  double *m = p->tmpTrsf->mat.m;

  size_t ncols = (theTrsf->vx).ncols; /* dimx */
  size_t nrows = (theTrsf->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;


  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  /*  0  1  2  3
   *  4  5  6  7
   *  8  9 10 11
   * 12 13 14 15
   */
  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    x = h[ 0] * theVecX[k][j][i] + h[ 1] * theVecY[k][j][i] + h[ 2] * theVecZ[k][j][i];
    y = h[ 4] * theVecX[k][j][i] + h[ 5] * theVecY[k][j][i] + h[ 6] * theVecZ[k][j][i];
    z = h[ 8] * theVecX[k][j][i] + h[ 9] * theVecY[k][j][i] + h[10] * theVecZ[k][j][i];
    resVecX[k][j][i] = m[ 0] * i + m [ 1] * j + m [ 2] * k + m[ 3] + x;
    resVecY[k][j][i] = m[ 4] * i + m [ 5] * j + m [ 6] * k + m[ 7] + y;
    resVecZ[k][j][i] = m[ 8] * i + m [ 9] * j + m [10] * k + m[11] + z;
  }
  chunk->ret = 1;
  return( (void*)NULL );
}







/* theTrsf allows to resample 'floIm' into 'refIm' thus goes
 * from 'refIm' to 'floIm'
*/
int BAL_ChangeTransformationToVoxelUnit( bal_image *floIm,
                                         bal_image *refIm,
                                         bal_transformation *theTrsf,
                                         bal_transformation *resTrsf )
{
  char *proc = "BAL_ChangeTransformationToVoxelUnit";
  bal_transformation tmpTrsf;

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _UnitTransformationChangeParam p;

  if ( floIm == (bal_image*)NULL || refIm == (bal_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL template images\n", proc );
    return( -1 );
  }

  if ( theTrsf->type != resTrsf->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations have different types\n", proc );
    return( -1 );
  }

  /* nothing to do, but copy the transformation
   */
  if ( theTrsf->transformation_unit == VOXEL_UNIT ) {
    if ( BAL_CopyTransformation( theTrsf, resTrsf ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy transformation\n", proc );
      return( -1 );
    }
    return( 1 );
  }

  /* something to do
   */

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
    return( -1 );
  }

  /* preparing parallelism for vectorfield
   */
  initChunks( &chunks );
  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    first = 0;
    last = theTrsf->vx.nplanes * theTrsf->vx.nrows * theTrsf->vx.ncols - 1;
    if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
    }

    /* floIm->to_voxel o refIm->to_real = tmpTrsf
     */
    _mult_mat( &(floIm->to_voxel), &(refIm->to_real), &(tmpTrsf.mat) );
    /* tmpTrsf = ( floIm->to_voxel o refIm->to_real - Id )
     */
    tmpTrsf.mat.m[ 0] -= 1.0;
    tmpTrsf.mat.m[ 5] -= 1.0;
    tmpTrsf.mat.m[10] -= 1.0;

    p.floIm = floIm;
    p.refIm = refIm;
    p.theTrsf = theTrsf;
    p.resTrsf = resTrsf;
    p.tmpTrsf = &tmpTrsf;
    p.tmpGeometry = _combinedGeometry( floIm->geometry, refIm->geometry );

    for ( i=0; i<chunks.n_allocated_chunks; i++ )
      chunks.data[i].parameters = (void*)(&p);
    break;
  }

  switch ( theTrsf->type ) {
  default :
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

    /* theTrsf o refIm->to_real = tmpTrsf
     */
    _mult_mat( &(theTrsf->mat), &(refIm->to_real), &(tmpTrsf.mat) );

    /* floIm->to_voxel o tmpTrsf = resTrsf
     */
    _mult_mat( &(floIm->to_voxel), &(tmpTrsf.mat), &(resTrsf->mat) );

    break;

  case VECTORFIELD_2D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols
         || resTrsf->vx.nrows != theTrsf->vx.nrows
         || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
                 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
                 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      freeChunks( &chunks );
      return( -1 );
    }

     if ( processChunks( &_Change2DVectorFieldToVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to change 2D vector field into voxel unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case VECTORFIELD_3D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols
         || resTrsf->vx.nrows != theTrsf->vx.nrows
         || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
                 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
                 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      freeChunks( &chunks );
      return( -1 );
    }

    if ( processChunks( &_Change3DVectorFieldToVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to change3D vector field into voxel unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  BAL_FreeTransformation( &tmpTrsf );

  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    freeChunks( &chunks );
    break;
  }

  resTrsf->transformation_unit = VOXEL_UNIT;

  return( 1 );

}

