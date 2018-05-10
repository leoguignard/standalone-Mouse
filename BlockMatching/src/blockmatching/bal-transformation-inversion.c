/*************************************************************************
 * bal-transformation-inversion.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Sam 18 fév 2017 14:13:33 CET
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
#include <vtmalloc.h>

#include <bal-transformation-copy.h>
#include <bal-transformation-inversion.h>

static int _verbose_ = 1;
static int _debug_ = 0;


void BAL_SetVerboseInBalTransformationInversion( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationInversion(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationInversion(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationInversion( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalTransformationInversion(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationInversion(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}



/*************************************************************
 *
 * static functions
 *
 ************************************************************/


static int BAL_Inverse2DVectorField( bal_transformation *theTrsf,
                                     bal_transformation *invTrsf );
static int BAL_Inverse3DVectorField( bal_transformation *theTrsf,
                                     bal_transformation *invTrsf );


/*--------------------------------------------------
 *
 * TRANSFORMATION INVERSION
 *
 --------------------------------------------------*/


/* this is mainly for tests
 * this allows to compare inversion of 'voxels' vector field
 * to the inversion of 'real' vector field
 */
static int _use_conversion_for_real_vector_field_ = 0;

void BAL_SetConversionForVectorFieldInversionInBalTransformationInversion( int c )
{
  _use_conversion_for_real_vector_field_ = c;
}



int BAL_InverseTransformation( bal_transformation *theTrsf,
                               bal_transformation *invTrsf )
{
  char *proc = "BAL_InverseTransformation";
  double *theMat = theTrsf->mat.m;
  double *invMat = invTrsf->mat.m;
  enumUnitTransfo transformation_unit;

  if ( theTrsf == invTrsf ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: output transformation has to be different from the input one\n", proc );
    }
    return( -1 );
  }



  if ( BAL_DoesTransformationExist( invTrsf ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: inverse transformation does not exist or is not allocated\n", proc );
    }
    return( -1 );
  }


  if ( _debug_ ) {
    fprintf( stderr, "%s: will invert transformation\n", proc );
    fprintf( stderr, "\t - type = " );
    BAL_PrintTransformationType( stderr, theTrsf->type );
    fprintf( stderr, "\t - unit = " );
    BAL_PrintTransformationUnit( stderr, theTrsf->transformation_unit );
  }

  /* ensure that result transformation has the same unit
   * than the input transformation
   */
  invTrsf->transformation_unit =  theTrsf->transformation_unit;


  switch ( theTrsf->type ) {

  default :

    if ( _verbose_ ) {
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
      BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
    }
    return( -1 );

  case TRANSLATION_2D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 3] = (- theMat[ 3]);
    invMat[ 7] = (- theMat[ 7]);
    break;

  case TRANSLATION_3D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 3] = (- theMat[ 3]);
    invMat[ 7] = (- theMat[ 7]);
    invMat[11] = (- theMat[11]);
    break;

  case TRANSLATION_SCALING_2D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 0] = 1.0/theMat[ 0];
    invMat[ 5] = 1.0/theMat[ 5];
    invMat[ 3] = (- theMat[ 3])/theMat[ 0];
    invMat[ 7] = (- theMat[ 7])/theMat[ 5];
    break;

  case TRANSLATION_SCALING_3D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 0] = 1.0/theMat[ 0];
    invMat[ 5] = 1.0/theMat[ 5];
    invMat[10] = 1.0/theMat[10];
    invMat[ 3] = (- theMat[ 3])/theMat[ 0];
    invMat[ 7] = (- theMat[ 7])/theMat[ 5];
    invMat[11] = (- theMat[11])/theMat[10];
    break;

  case RIGID_2D :
    /* l'inverse d'une matrice de rotation R est sa transposee R^t
       la translation inverse est -R^t * T
     */
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 1] = theMat[ 4];
    invMat[ 4] = theMat[ 1];
    invMat[ 3] = ( - (theMat[ 0]*theMat[ 3] + theMat[ 4]*theMat[ 7]) );
    invMat[ 7] = ( - (theMat[ 1]*theMat[ 3] + theMat[ 5]*theMat[ 7]) );
    break;

  case RIGID_3D :
    /* l'inverse d'une matrice de rotation R est sa transposee R^t
       la translation inverse est -R^t * T
     */
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 1] = theMat[ 4];
    invMat[ 2] = theMat[ 8];
    invMat[ 4] = theMat[ 1];
    invMat[ 6] = theMat[ 9];
    invMat[ 8] = theMat[ 2];
    invMat[ 9] = theMat[ 6];
    invMat[ 3] = ( - (theMat[ 0]*theMat[ 3] + theMat[ 4]*theMat[ 7] + theMat[ 8]*theMat[11]) );
    invMat[ 7] = ( - (theMat[ 1]*theMat[ 3] + theMat[ 5]*theMat[ 7] + theMat[ 9]*theMat[11]) );
    invMat[11] = ( - (theMat[ 2]*theMat[ 3] + theMat[ 6]*theMat[ 7] + theMat[10]*theMat[11]) );
    break;

  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    if ( InverseMat4x4( theMat, invMat ) != 4 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: input transformation is not invertible (linear case)\n", proc );
        BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }
    break;

  case VECTORFIELD_2D :
    if ( _use_conversion_for_real_vector_field_ && theTrsf->transformation_unit == REAL_UNIT ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: 2D vector field is converted into voxel units for inversion\n", proc );
      transformation_unit = theTrsf->transformation_unit;
      (void)BAL_ChangeTransformationToVoxelUnit( &(invTrsf->vx), &(theTrsf->vx), theTrsf, theTrsf );
      invTrsf->transformation_unit = VOXEL_UNIT;
    }

    if ( BAL_Inverse2DVectorField( theTrsf, invTrsf ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: input transformation is not invertible (2D vector field)\n", proc );
        BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }

    if ( _use_conversion_for_real_vector_field_ && transformation_unit == REAL_UNIT ) {
      (void)BAL_ChangeTransformationToRealUnit( &(invTrsf->vx), &(theTrsf->vx), theTrsf, theTrsf );
      (void)BAL_ChangeTransformationToRealUnit( &(theTrsf->vx), &(invTrsf->vx), invTrsf, invTrsf );
    }
    break;

  case VECTORFIELD_3D :
    if ( _use_conversion_for_real_vector_field_ && theTrsf->transformation_unit == REAL_UNIT ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: 3D vector field is converted into voxel units for inversion\n", proc  );
      transformation_unit = theTrsf->transformation_unit;
      (void)BAL_ChangeTransformationToVoxelUnit( &(invTrsf->vx), &(theTrsf->vx), theTrsf, theTrsf );
      invTrsf->transformation_unit = VOXEL_UNIT;
    }

    if ( BAL_Inverse3DVectorField( theTrsf, invTrsf ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: input transformation is not invertible (3D vector field)\n", proc );
        BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }

    if ( _use_conversion_for_real_vector_field_ && transformation_unit == REAL_UNIT ) {
      (void)BAL_ChangeTransformationToRealUnit( &(invTrsf->vx), &(theTrsf->vx), theTrsf, theTrsf );
      (void)BAL_ChangeTransformationToRealUnit( &(theTrsf->vx), &(invTrsf->vx), invTrsf, invTrsf );
    }
    break;

  }

  return( 1 );

}








/********************************************************************************
 *
 * vector field inversion
 *
 ********************************************************************************/



static double alpha = 0.5;
static double derivationSigma = 2.0;
static double ERRMAX = 0.05;
static int ITERMAX = 10;

static enumVectorFieldInverseInitialization initialization = FORWARD_INTERPOLATION;
static double forwardSigma = 2.0;

void BAL_SetDerivationSigmaForVectorFieldInversionInBalTransformationInversion( double s )
{
    derivationSigma = s;
}

double BAL_GetDerivationSigmaForVectorFieldInversionInBalTransformationInversion()
{
    return( derivationSigma );
}

void BAL_SetIterationsMaxForVectorFieldInversionInBalTransformationInversion( int i )
{
    ITERMAX = i;
}

int BAL_GetIterationsMaxForVectorFieldInversionInBalTransformationInversion( )
{
    return( ITERMAX );
}

void BAL_SetErrorMaxForVectorFieldInversionInBalTransformationInversion( double e )
{
    ERRMAX = e;
}

double BAL_GetErrorMaxForVectorFieldInversionInBalTransformationInversion( )
{
    return( ERRMAX );
}

static bal_image *imInverseErrors = (bal_image *)NULL;

void BAL_SetImageInverseErrorsForVectorFieldInversionInBalTransformationInversion( bal_image *i )
{
    imInverseErrors = i;
}

void BAL_SetInitializationForVectorFieldInversionInBalTransformationInversion( enumVectorFieldInverseInitialization i )
{
    initialization = i;
}

enumVectorFieldInverseInitialization BAL_GetInitializationForVectorFieldInversionInBalTransformationInversion()
{
    return( initialization );
}

void BAL_SetForwardSigmaForVectorFieldInversionInBalTransformationInversion( double s )
{
    forwardSigma = s;
}

double BAL_GetForwardSigmaForVectorFieldInversionInBalTransformationInversion()
{
    return( forwardSigma );
}




static void _get2DMatrixVector( double x, double y, int k,
                                float ***xdx, float ***xdy, float ***ydx, float ***ydy,
                                float ***vecx, float*** vecy,
                                int *dim,
                                double *mat, double *vec )
{
  int ix, iy;
  double dx, dy;
  double a00, a01, a10, a11;

  ix = (int)x;
  iy = (int)y;
  dx = x - ix;
  dy = y - iy;

  if ( ix < 0 ) {
    ix = 0;
    dx = 0;
  }
  if ( iy < 0 ) {
    iy = 0;
    dy = 0;
  }

  if ( ix >= dim[0] - 1 ) {
    ix = dim[0]-2;
    dx = 1;
  }
  if ( iy >= dim[1] - 1 ) {
    iy = dim[1]-2;
    dy = 1;
  }

  /* res = (1-dx) * ( (1-dy) * buf[k][iy][ix]
                      + (dy) * buf[k][iy+1][ix] )
           + (dx) * ( (1-dy) * buf[k][iy][ix+1]
                             + (dy) * buf[k][iy+1][ix+1] );
         =   a00 * buf[k][iy][ix]
           + a10 * buf[k][iy+1][ix]
           + a01 * buf[k][iy][ix+1]
                  + a11 * buf[k][iy+1][ix+1];
     a00 = (1-dx) * (1-dy);
     a10 = (1-dx) *    dy;
     a01 =    dx  * (1-dy);
     a11 =    dx  *    dy
  */
  a11 = dx * dy;
  a01 = dx - a11;
  a10 = dy - a11;
  a00 = 1 - dx - a10;

  mat[0] = a00 * xdx[k][iy][ix] + a10 * xdx[k][iy+1][ix] + a01 * xdx[k][iy][ix+1]  + a11 * xdx[k][iy+1][ix+1];
  mat[1] = a00 * xdy[k][iy][ix] + a10 * xdy[k][iy+1][ix] + a01 * xdy[k][iy][ix+1]  + a11 * xdy[k][iy+1][ix+1];
  mat[2] = a00 * ydx[k][iy][ix] + a10 * ydx[k][iy+1][ix] + a01 * ydx[k][iy][ix+1]  + a11 * ydx[k][iy+1][ix+1];
  mat[3] = a00 * ydy[k][iy][ix] + a10 * ydy[k][iy+1][ix] + a01 * ydy[k][iy][ix+1]  + a11 * ydy[k][iy+1][ix+1];

  vec[0] = a00 * vecx[k][iy][ix] + a10 * vecx[k][iy+1][ix] + a01 * vecx[k][iy][ix+1]  + a11 * vecx[k][iy+1][ix+1];
  vec[1] = a00 * vecy[k][iy][ix] + a10 * vecy[k][iy+1][ix] + a01 * vecy[k][iy][ix+1]  + a11 * vecy[k][iy+1][ix+1];
}





static void _get3DMatrixVector( double x, double y, double z,
                                float ***xdx, float ***xdy, float ***xdz,
                                float ***ydx, float ***ydy, float ***ydz,
                                float ***zdx, float ***zdy, float ***zdz,
                                float ***vecx, float*** vecy, float*** vecz,
                                int *dim,
                                double *mat, double *vec )
{
  int ix, iy, iz;
  double dx, dy, dz;
  double a000, a001, a010, a011, a100, a101, a110, a111;
  double dxdy, dxdz, dydz;

  ix = (int)x;
  iy = (int)y;
  iz = (int)z;
  dx = x - ix;
  dy = y - iy;
  dz = z - iz;

  if ( ix < 0 ) {
    ix = 0;
    dx = 0;
  }
  if ( iy < 0 ) {
    iy = 0;
    dy = 0;
  }
  if ( iz < 0 ) {
    iz = 0;
    dz = 0;
  }

  if ( ix >= dim[0] - 1 ) {
    ix = dim[0]-2;
    dx = 1;
  }
  if ( iy >= dim[1] - 1 ) {
    iy = dim[1]-2;
    dy = 1;
  }
  if ( iz >= dim[2] - 1 ) {
    iz = dim[2]-2;
    dz = 1;
  }

  /* res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                  + (dz)*buf[iz+1][iy][ix] )
                       + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                  + (dz)*buf[iz+1][iy+1][ix] ) )
            + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                  + (dz)*buf[iz+1][iy][ix+1] )
                       + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                  + (dz)*buf[iz+1][iy+1][ix+1] ) );

         =   a000 * buf[iz][iy][ix]
           + a100 * buf[iz+1][iy][ix]
           + a010 * buf[iz][iy+1][ix]
           + a110 * buf[iz+1][iy+1][ix]
           + a001 * buf[iz][iy][ix+1]
           + a101 * buf[iz+1][iy][ix+1]
           + a011 * buf[iz][iy+1][ix+1]
           + a111 * buf[iz+1][iy+1][ix+1]

     a000 = (1-dx) * (1-dy) * (1-dz);
     a100 = (1-dx) * (1-dy) *    dz;
     a010 = (1-dx) *    dy  * (1-dz);
     a110 = (1-dx) *    dy  *    dz;
     a001 =    dx  * (1-dy) * (1-dz);
     a101 =    dx  * (1-dy) *    dz;
     a011 =    dx  *    dy  * (1-dz);
     a111 =    dx  *    dy  *    dz;

  */

  dxdy = dx*dy;
  dxdz = dx*dz;
  dydz = dy*dz;

  a111 = dxdy * dz;
  a011 = dxdy - a111;
  a101 = dxdz - a111;
  a110 = dydz - a111;
  a001 = dx - dxdy - a101;
  a010 = dy - dydz - a011;
  a100 = dz - dxdz - a110;
  a000 = 1 - dy -dz + dydz - a001;


  mat[0] = a000 * xdx[iz][iy][ix] + a100 * xdx[iz+1][iy][ix] + a010 * xdx[iz][iy+1][ix]  + a110 * xdx[iz+1][iy+1][ix]
    + a001 * xdx[iz][iy][ix+1] + a101 * xdx[iz+1][iy][ix+1] + a011 * xdx[iz][iy+1][ix+1]  + a111 * xdx[iz+1][iy+1][ix+1];
  mat[1] = a000 * xdy[iz][iy][ix] + a100 * xdy[iz+1][iy][ix] + a010 * xdy[iz][iy+1][ix]  + a110 * xdy[iz+1][iy+1][ix]
    + a001 * xdy[iz][iy][ix+1] + a101 * xdy[iz+1][iy][ix+1] + a011 * xdy[iz][iy+1][ix+1]  + a111 * xdy[iz+1][iy+1][ix+1];
  mat[2] = a000 * xdz[iz][iy][ix] + a100 * xdz[iz+1][iy][ix] + a010 * xdz[iz][iy+1][ix]  + a110 * xdz[iz+1][iy+1][ix]
    + a001 * xdz[iz][iy][ix+1] + a101 * xdz[iz+1][iy][ix+1] + a011 * xdz[iz][iy+1][ix+1]  + a111 * xdz[iz+1][iy+1][ix+1];

  mat[3] = a000 * ydx[iz][iy][ix] + a100 * ydx[iz+1][iy][ix] + a010 * ydx[iz][iy+1][ix]  + a110 * ydx[iz+1][iy+1][ix]
    + a001 * ydx[iz][iy][ix+1] + a101 * ydx[iz+1][iy][ix+1] + a011 * ydx[iz][iy+1][ix+1]  + a111 * ydx[iz+1][iy+1][ix+1];
  mat[4] = a000 * ydy[iz][iy][ix] + a100 * ydy[iz+1][iy][ix] + a010 * ydy[iz][iy+1][ix]  + a110 * ydy[iz+1][iy+1][ix]
    + a001 * ydy[iz][iy][ix+1] + a101 * ydy[iz+1][iy][ix+1] + a011 * ydy[iz][iy+1][ix+1]  + a111 * ydy[iz+1][iy+1][ix+1];
  mat[5] = a000 * ydz[iz][iy][ix] + a100 * ydz[iz+1][iy][ix] + a010 * ydz[iz][iy+1][ix]  + a110 * ydz[iz+1][iy+1][ix]
    + a001 * ydz[iz][iy][ix+1] + a101 * ydz[iz+1][iy][ix+1] + a011 * ydz[iz][iy+1][ix+1]  + a111 * ydz[iz+1][iy+1][ix+1];

  mat[6] = a000 * zdx[iz][iy][ix] + a100 * zdx[iz+1][iy][ix] + a010 * zdx[iz][iy+1][ix]  + a110 * zdx[iz+1][iy+1][ix]
    + a001 * zdx[iz][iy][ix+1] + a101 * zdx[iz+1][iy][ix+1] + a011 * zdx[iz][iy+1][ix+1]  + a111 * zdx[iz+1][iy+1][ix+1];
  mat[7] = a000 * zdy[iz][iy][ix] + a100 * zdy[iz+1][iy][ix] + a010 * zdy[iz][iy+1][ix]  + a110 * zdy[iz+1][iy+1][ix]
    + a001 * zdy[iz][iy][ix+1] + a101 * zdy[iz+1][iy][ix+1] + a011 * zdy[iz][iy+1][ix+1]  + a111 * zdy[iz+1][iy+1][ix+1];
  mat[8] = a000 * zdz[iz][iy][ix] + a100 * zdz[iz+1][iy][ix] + a010 * zdz[iz][iy+1][ix]  + a110 * zdz[iz+1][iy+1][ix]
    + a001 * zdz[iz][iy][ix+1] + a101 * zdz[iz+1][iy][ix+1] + a011 * zdz[iz][iy+1][ix+1]  + a111 * zdz[iz+1][iy+1][ix+1];

  vec[0] = a000 * vecx[iz][iy][ix] + a100 * vecx[iz+1][iy][ix] + a010 * vecx[iz][iy+1][ix]  + a110 * vecx[iz+1][iy+1][ix]
    + a001 * vecx[iz][iy][ix+1] + a101 * vecx[iz+1][iy][ix+1] + a011 * vecx[iz][iy+1][ix+1]  + a111 * vecx[iz+1][iy+1][ix+1];
  vec[1] = a000 * vecy[iz][iy][ix] + a100 * vecy[iz+1][iy][ix] + a010 * vecy[iz][iy+1][ix]  + a110 * vecy[iz+1][iy+1][ix]
    + a001 * vecy[iz][iy][ix+1] + a101 * vecy[iz+1][iy][ix+1] + a011 * vecy[iz][iy+1][ix+1]  + a111 * vecy[iz+1][iy+1][ix+1];
  vec[2] = a000 * vecz[iz][iy][ix] + a100 * vecz[iz+1][iy][ix] + a010 * vecz[iz][iy+1][ix]  + a110 * vecz[iz+1][iy+1][ix]
    + a001 * vecz[iz][iy][ix+1] + a101 * vecz[iz+1][iy][ix+1] + a011 * vecz[iz][iy+1][ix+1]  + a111 * vecz[iz+1][iy+1][ix+1];

}





typedef struct _TransformationInversionParam {

    int theTrsfDim[3];
    int invTrsfDim[3];

    bal_imageGeometry theGeometry;
    bal_imageGeometry invGeometry;

    _MATRIX *the_to_voxel;
    _MATRIX *inv_to_real;

    float ***arrXdX;
    float ***arrXdY;
    float ***arrXdZ;
    float ***arrYdX;
    float ***arrYdY;
    float ***arrYdZ;
    float ***arrZdX;
    float ***arrZdY;
    float ***arrZdZ;

    float ***arrInvX;
    float ***arrInvY;
    float ***arrInvZ;

    float ***arrTrsfX;
    float ***arrTrsfY;
    float ***arrTrsfZ;

    u8 ***theErrors;

    int ndivergence;
    int nnonconvergence;
} _TransformationInversionParam;





static void _initTransformationInversionParam( _TransformationInversionParam *p )
{
    p->theTrsfDim[0] = 0;
    p->theTrsfDim[1] = 0;
    p->theTrsfDim[2] = 0;

    p->invTrsfDim[0] = 0;
    p->invTrsfDim[1] = 0;
    p->invTrsfDim[2] = 0;

    p->theGeometry = _BAL_UNKNOWN_GEOMETRY_;
    p->invGeometry = _BAL_UNKNOWN_GEOMETRY_;


    p->the_to_voxel = (_MATRIX*)NULL;
    p->inv_to_real = (_MATRIX*)NULL;

    p->arrXdX = (float***)NULL;
    p->arrXdY = (float***)NULL;
    p->arrXdZ = (float***)NULL;
    p->arrYdX = (float***)NULL;
    p->arrYdY = (float***)NULL;
    p->arrYdZ = (float***)NULL;
    p->arrZdX = (float***)NULL;
    p->arrZdY = (float***)NULL;
    p->arrZdZ = (float***)NULL;

    p->arrInvX = (float***)NULL;
    p->arrInvY = (float***)NULL;
    p->arrInvZ = (float***)NULL;

    p->arrTrsfX = (float***)NULL;
    p->arrTrsfY = (float***)NULL;
    p->arrTrsfZ = (float***)NULL;

    p->theErrors = (u8***)NULL;

    p->ndivergence = 0;
    p->nnonconvergence = 0;
}





static void *_Inverse2DVoxelMatrix( void *par )
{
    char *proc = "_Inverse2DVoxelMatrix";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;


    size_t dimx = p->theTrsfDim[0];
    size_t dimy = p->theTrsfDim[1];
    size_t dimz = p->theTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;

    double tmpMat[4], det;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    /* add Identity to V.Nabla^t
       and invert (Id + V.Nabla^t)
       = ( xdx + 1     xdy     )
         ( ydx         ydy + 1 )

       The inverse is computed with
       | a11 a12 |-1             |  a22 -a12 |
       | a21 a22 |    =  1/DET * | -a21  a11 |

       with DET  =  a11 a22- a12 a21

     */

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                arrXdX[k][j][i] += 1.0;
                arrYdY[k][j][i] += 1.0;

                det = arrXdX[k][j][i]*arrYdY[k][j][i] - arrXdY[k][j][i] * arrYdX[k][j][i];
                tmpMat[0] = arrYdY[k][j][i] / det;
                tmpMat[1] = (- arrXdY[k][j][i]) / det;
                tmpMat[2] = (- arrYdX[k][j][i]) / det;
                tmpMat[3] = arrXdX[k][j][i] / det;

                arrXdX[k][j][i] = tmpMat[0];
                arrXdY[k][j][i] = tmpMat[1];
                arrYdX[k][j][i] = tmpMat[2];
                arrYdY[k][j][i] = tmpMat[3];
            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse3DVoxelMatrix( void *par )
{
    char *proc = "_Inverse3DVoxelMatrix";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrXdZ = p->arrXdZ;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;
    float ***arrYdZ = p->arrYdZ;
    float ***arrZdX = p->arrZdX;
    float ***arrZdY = p->arrZdY;
    float ***arrZdZ = p->arrZdZ;


    size_t dimx = p->theTrsfDim[0];
    size_t dimy = p->theTrsfDim[1];
    size_t dimz = p->theTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;

    double tmpMat[9], det;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    /* add Identity to V.Nabla^t
       and invert (Id + V.Nabla^t)
       = ( xdx + 1     xdy         xdz     )
         ( ydx         ydy + 1     ydz     )
         ( zdx         zdy         zdz + 1 )

       The inverse is computed with
       | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
       | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
       | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

       with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

     */

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                arrXdX[k][j][i] += 1.0;
                arrYdY[k][j][i] += 1.0;
                arrZdZ[k][j][i] += 1.0;

                det = arrXdX[k][j][i] * ( arrYdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrYdZ[k][j][i] )
                  - arrYdX[k][j][i] * ( arrXdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrXdZ[k][j][i] )
                  + arrZdX[k][j][i] * ( arrXdY[k][j][i]*arrYdZ[k][j][i] - arrYdY[k][j][i]*arrXdZ[k][j][i] );
                tmpMat[0] = ( arrYdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrYdZ[k][j][i] ) / det;
                tmpMat[1] = ( - arrXdY[k][j][i]*arrZdZ[k][j][i] + arrZdY[k][j][i]*arrXdZ[k][j][i] ) / det;
                tmpMat[2] = ( arrXdY[k][j][i]*arrYdZ[k][j][i] - arrYdY[k][j][i]*arrXdZ[k][j][i] ) / det;
                tmpMat[3] = ( - arrYdX[k][j][i]*arrZdZ[k][j][i] + arrZdX[k][j][i]*arrYdZ[k][j][i] ) / det;
                tmpMat[4] = ( arrXdX[k][j][i]*arrZdZ[k][j][i] - arrZdX[k][j][i]*arrXdZ[k][j][i] ) / det;
                tmpMat[5] = ( - arrXdX[k][j][i]*arrYdZ[k][j][i] + arrYdX[k][j][i]*arrXdZ[k][j][i] ) / det;
                tmpMat[6] = ( arrYdX[k][j][i]*arrZdY[k][j][i] - arrZdX[k][j][i]*arrYdY[k][j][i] ) / det;
                tmpMat[7] = ( - arrXdX[k][j][i]*arrZdY[k][j][i] + arrZdX[k][j][i]*arrXdY[k][j][i] ) / det;
                tmpMat[8] = ( arrXdX[k][j][i]*arrYdY[k][j][i] - arrYdX[k][j][i]*arrXdY[k][j][i] ) / det;

                arrXdX[k][j][i] = tmpMat[0];
                arrXdY[k][j][i] = tmpMat[1];
                arrXdZ[k][j][i] = tmpMat[2];
                arrYdX[k][j][i] = tmpMat[3];
                arrYdY[k][j][i] = tmpMat[4];
                arrYdZ[k][j][i] = tmpMat[5];
                arrZdX[k][j][i] = tmpMat[6];
                arrZdY[k][j][i] = tmpMat[7];
                arrZdZ[k][j][i] = tmpMat[8];
            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse2DRealMatrix( void *par )
{
    char *proc = "_Inverse2DRealMatrix";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;


    size_t dimx = p->theTrsfDim[0];
    size_t dimy = p->theTrsfDim[1];
    size_t dimz = p->theTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;

    double tmp[4], mat[4], det;
    double *to_voxel = p->the_to_voxel->m;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    /* add Identity to V.Nabla^t
       and invert (Id + V.Nabla^t)
       = ( xdx + 1     xdy     )
         ( ydx         ydy + 1 )

       The inverse is computed with
       | a11 a12 |-1             |  a22 -a12 |
       | a21 a22 |    =  1/DET * | -a21  a11 |

       with DET  =  a11 a22- a12 a21

                 ( 0  1  2 )
                 ( 4  5  6 )
                 ( 8  9 10 )
       ( 0 1 )
       ( 2 3 )

     */

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                tmp[0] = arrXdX[k][j][i] + 1.0;
                tmp[1] = arrXdY[k][j][i];

                tmp[2] = arrYdX[k][j][i];
                tmp[3] = arrYdY[k][j][i] + 1.0;

                switch( p->theGeometry ) {
                default :
                case _BAL_UNKNOWN_GEOMETRY_ :
                case _BAL_HOMOTHETY_GEOMETRY_ :
                case _BAL_TRANSLATION_GEOMETRY_ :
                    mat[0] = tmp[0] * to_voxel[ 0];
                    mat[1] =                         tmp[1] * to_voxel[ 5];

                    mat[2] = tmp[2] * to_voxel[ 0];
                    mat[3] =                         tmp[3] * to_voxel[ 5];
                    break;
                case _BAL_QFORM_GEOMETRY_ :
                    mat[0] = tmp[0] * to_voxel[ 0] + tmp[1] * to_voxel[ 4];
                    mat[1] = tmp[0] * to_voxel[ 1] + tmp[1] * to_voxel[ 5];

                    mat[2] = tmp[2] * to_voxel[ 0] + tmp[3] * to_voxel[ 4];
                    mat[3] = tmp[2] * to_voxel[ 1] + tmp[3] * to_voxel[ 5];
                    break;
                }

                det = mat[0] * mat[3] - mat[1] * mat[2];

                tmp[0] = (  mat[3] ) / det;
                tmp[1] = (- mat[1] ) / det;
                tmp[2] = (- mat[2] ) / det;
                tmp[3] = (  mat[0] ) / det;

                arrXdX[k][j][i] = tmp[0];
                arrXdY[k][j][i] = tmp[1];
                arrYdX[k][j][i] = tmp[2];
                arrYdY[k][j][i] = tmp[3];;
            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse3DRealMatrix( void *par )
{
    char *proc = "_Inverse3DRealMatrix";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrXdZ = p->arrXdZ;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;
    float ***arrYdZ = p->arrYdZ;
    float ***arrZdX = p->arrZdX;
    float ***arrZdY = p->arrZdY;
    float ***arrZdZ = p->arrZdZ;


    size_t dimx = p->theTrsfDim[0];
    size_t dimy = p->theTrsfDim[1];
    size_t dimz = p->theTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;

    double tmp[9], mat[9], det;
    double *to_voxel = p->the_to_voxel->m;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    /* add Identity to V.Nabla^t
       and invert (Id + V.Nabla^t)
       = ( xdx + 1     xdy         xdz     )
         ( ydx         ydy + 1     ydz     )
         ( zdx         zdy         zdz + 1 )

       The inverse is computed with
       | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
       | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
       | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

       with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

                 ( 0  1  2 )
                 ( 4  5  6 )
                 ( 8  9 10 )
       ( 0 1 2 )
       ( 3 4 5 )
       ( 6 7 8 )

     */

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                tmp[0] = arrXdX[k][j][i] + 1.0;
                tmp[1] = arrXdY[k][j][i];
                tmp[2] = arrXdZ[k][j][i];

                tmp[3] = arrYdX[k][j][i];
                tmp[4] = arrYdY[k][j][i] + 1.0;
                tmp[5] = arrYdZ[k][j][i];

                tmp[6] = arrZdX[k][j][i];
                tmp[7] = arrZdY[k][j][i];
                tmp[8] = arrZdZ[k][j][i] + 1.0;

                switch( p->theGeometry ) {
                default :
                case _BAL_UNKNOWN_GEOMETRY_ :
                case _BAL_HOMOTHETY_GEOMETRY_ :
                case _BAL_TRANSLATION_GEOMETRY_ :
                    mat[0] = tmp[0] * to_voxel[ 0];
                    mat[1] =                         tmp[1] * to_voxel[ 5];
                    mat[2] =                                                 tmp[2] * to_voxel[10];

                    mat[3] = tmp[3] * to_voxel[ 0];
                    mat[4] =                         tmp[4] * to_voxel[ 5];
                    mat[5] =                                                 tmp[5] * to_voxel[10];

                    mat[6] = tmp[6] * to_voxel[ 0];
                    mat[7] =                         tmp[7] * to_voxel[ 5];
                    mat[8] =                                                 tmp[8] * to_voxel[10];
                    break;
                case _BAL_QFORM_GEOMETRY_ :
                    mat[0] = tmp[0] * to_voxel[ 0] + tmp[1] * to_voxel[ 4] + tmp[2] * to_voxel[ 8];
                    mat[1] = tmp[0] * to_voxel[ 1] + tmp[1] * to_voxel[ 5] + tmp[2] * to_voxel[ 9];
                    mat[2] = tmp[0] * to_voxel[ 2] + tmp[1] * to_voxel[ 6] + tmp[2] * to_voxel[10];

                    mat[3] = tmp[3] * to_voxel[ 0] + tmp[4] * to_voxel[ 4] + tmp[5] * to_voxel[ 8];
                    mat[4] = tmp[3] * to_voxel[ 1] + tmp[4] * to_voxel[ 5] + tmp[5] * to_voxel[ 9];
                    mat[5] = tmp[3] * to_voxel[ 2] + tmp[4] * to_voxel[ 6] + tmp[5] * to_voxel[10];

                    mat[6] = tmp[6] * to_voxel[ 0] + tmp[7] * to_voxel[ 4] + tmp[8] * to_voxel[ 8];
                    mat[7] = tmp[6] * to_voxel[ 1] + tmp[7] * to_voxel[ 5] + tmp[8] * to_voxel[ 9];
                    mat[8] = tmp[6] * to_voxel[ 2] + tmp[7] * to_voxel[ 6] + tmp[8] * to_voxel[10];
                    break;
                }

                det = mat[0] * ( mat[4]*mat[8] - mat[5]*mat[7] )
                    - mat[3] * ( mat[1]*mat[8] - mat[2]*mat[7] )
                    + mat[6] * ( mat[1]*mat[5] - mat[2]*mat[4] );
                tmp[0] = ( mat[4]*mat[8] - mat[5]*mat[7] ) / det;
                tmp[1] = ( mat[2]*mat[7] - mat[1]*mat[8] ) / det;
                tmp[2] = ( mat[1]*mat[5] - mat[2]*mat[4] ) / det;

                tmp[3] = ( mat[5]*mat[6] - mat[3]*mat[8] ) / det;
                tmp[4] = ( mat[0]*mat[8] - mat[2]*mat[6] ) / det;
                tmp[5] = ( mat[2]*mat[3] - mat[0]*mat[5] ) / det;

                tmp[6] = ( mat[3]*mat[7] - mat[4]*mat[6] ) / det;
                tmp[7] = ( mat[1]*mat[6] - mat[0]*mat[7] ) / det;
                tmp[8] = ( mat[0]*mat[4] - mat[1]*mat[3] ) / det;


                arrXdX[k][j][i] = tmp[0];
                arrXdY[k][j][i] = tmp[1];
                arrXdZ[k][j][i] = tmp[2];
                arrYdX[k][j][i] = tmp[3];
                arrYdY[k][j][i] = tmp[4];
                arrYdZ[k][j][i] = tmp[5];
                arrZdX[k][j][i] = tmp[6];
                arrZdY[k][j][i] = tmp[7];
                arrZdZ[k][j][i] = tmp[8];
            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse2DVoxelVectorField( void *par )
{
    char *proc = "_Inverse2DVoxelVectorField";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;

    float ***arrInvX = p->arrInvX;
    float ***arrInvY = p->arrInvY;

    float ***arrTrsfX = p->arrTrsfX;
    float ***arrTrsfY = p->arrTrsfY;

    u8 ***theErrors = p->theErrors;

    size_t dimx = p->invTrsfDim[0];
    size_t dimy = p->invTrsfDim[1];
    size_t dimz = p->invTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;
    int iterations;

    double tmpMat[4], tmpVec[2];
    double x, y;
    double dx, dy;
    double pe = 0, e = 0;

    float pInvX, pInvY;

    p->ndivergence = 0;
    p->nnonconvergence = 0;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                e = pe = 0.0;

                for ( iterations = 0; iterations< ITERMAX; iterations++ ) {

                  /* M + I(M) in voxel
                   */
                  x = i + arrInvX[k][j][i];
                  y = j + arrInvY[k][j][i];

                  /* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
                     and the vector v at M + I(M)
                  */

                  _get2DMatrixVector( x, y, k,
                                      arrXdX, arrXdY,
                                      arrYdX, arrYdY,
                                      arrTrsfX, arrTrsfY,
                                      p->theTrsfDim, tmpMat, tmpVec );
                  /* I(M) + v(M+I(M))
                   */
                  x = arrInvX[k][j][i] + tmpVec[0];
                  y = arrInvY[k][j][i] + tmpVec[1];

                  /* first tests on |I(M) + v(M+I(M))|
                   * - is it small enough?
                   * - does it increase ? -> retrieve previous values
                   */
                  e = fabs( x ) + fabs( y );
                  if ( e <= ERRMAX )
                    break;
                  if ( iterations > 0 && e > pe ) {
                      p->ndivergence ++;
                      arrInvX[k][j][i] = pInvX;
                      arrInvY[k][j][i] = pInvY;
                      break;
                  }
                  pe = e;

                  /* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
                     is the optimal variation of I(M)
                  */
                  dx = tmpMat[0] * x + tmpMat[1] * y;
                  dy = tmpMat[2] * x + tmpMat[3] * y;

                  pInvX = arrInvX[k][j][i];
                  pInvY = arrInvY[k][j][i];

                  arrInvX[k][j][i] -= alpha * dx;
                  arrInvY[k][j][i] -= alpha * dy;

                }

                if ( iterations >= ITERMAX ) {
                    p->nnonconvergence ++;
                }

                if ( iterations >= ITERMAX || (iterations > 0 && e > pe) ) {
                    if ( theErrors != (u8***)NULL )
                        theErrors[k][j][i] = 255;
                }

            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse3DVoxelVectorField( void *par )
{
    char *proc = "_Inverse3DVoxelVectorField";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrXdZ = p->arrXdZ;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;
    float ***arrYdZ = p->arrYdZ;
    float ***arrZdX = p->arrZdX;
    float ***arrZdY = p->arrZdY;
    float ***arrZdZ = p->arrZdZ;

    float ***arrInvX = p->arrInvX;
    float ***arrInvY = p->arrInvY;
    float ***arrInvZ = p->arrInvZ;

    float ***arrTrsfX = p->arrTrsfX;
    float ***arrTrsfY = p->arrTrsfY;
    float ***arrTrsfZ = p->arrTrsfZ;

    u8 ***theErrors = p->theErrors;

    size_t dimx = p->invTrsfDim[0];
    size_t dimy = p->invTrsfDim[1];
    size_t dimz = p->invTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;
    int iterations;

    double tmpMat[9], tmpVec[3];
    double x, y, z;
    double dx, dy, dz;
    double pe = 0, e = 0;

    float pInvX, pInvY, pInvZ;

    p->ndivergence = 0;
    p->nnonconvergence = 0;

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                e = pe = 0.0;

                for ( iterations = 0; iterations< ITERMAX; iterations++ ) {

                  /* M + I(M) in voxel
                   */
                  x = i + arrInvX[k][j][i];
                  y = j + arrInvY[k][j][i];
                  z = k + arrInvZ[k][j][i];

                  /* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
                     and the vector v at M + I(M)
                  */

                  _get3DMatrixVector( x, y, z,
                                      arrXdX, arrXdY, arrXdZ,
                                      arrYdX, arrYdY, arrYdZ,
                                      arrZdX, arrZdY, arrZdZ,
                                      arrTrsfX, arrTrsfY, arrTrsfZ,
                                      p->theTrsfDim, tmpMat, tmpVec );
                  /* I(M) + v(M+I(M))
                   */
                  x = arrInvX[k][j][i] + tmpVec[0];
                  y = arrInvY[k][j][i] + tmpVec[1];
                  z = arrInvZ[k][j][i] + tmpVec[2];

                  /* first tests on |I(M) + v(M+I(M))|
                   * - is it small enough?
                   * - does it increase ? -> retrieve previous values
                   */
                  e = fabs( x ) + fabs( y ) + fabs( z );
                  if ( e <= ERRMAX )
                    break;
                  if ( iterations > 0 && e > pe ) {
                      p->ndivergence ++;
                      arrInvX[k][j][i] = pInvX;
                      arrInvY[k][j][i] = pInvY;
                      arrInvZ[k][j][i] = pInvZ;
                      break;
                  }
                  pe = e;

                  /* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
                     is the optimal variation of I(M)
                  */
                  dx = tmpMat[0] * x + tmpMat[1] * y + tmpMat[2] * z;
                  dy = tmpMat[3] * x + tmpMat[4] * y + tmpMat[5] * z;
                  dz = tmpMat[6] * x + tmpMat[7] * y + tmpMat[8] * z;

                  pInvX = arrInvX[k][j][i];
                  pInvY = arrInvY[k][j][i];
                  pInvZ = arrInvZ[k][j][i];

                  arrInvX[k][j][i] -= alpha * dx;
                  arrInvY[k][j][i] -= alpha * dy;
                  arrInvZ[k][j][i] -= alpha * dz;

                }

                if ( iterations >= ITERMAX ) {
                    p->nnonconvergence ++;
                }

                if ( iterations >= ITERMAX || (iterations > 0 && e > pe) ) {
                    if ( theErrors != (u8***)NULL )
                        theErrors[k][j][i] = 255;
                }

            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse2DRealVectorField( void *par )
{
    char *proc = "_Inverse2DRealVectorField";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;

    float ***arrInvX = p->arrInvX;
    float ***arrInvY = p->arrInvY;

    float ***arrTrsfX = p->arrTrsfX;
    float ***arrTrsfY = p->arrTrsfY;

    u8 ***theErrors = p->theErrors;

    size_t dimx = p->invTrsfDim[0];
    size_t dimy = p->invTrsfDim[1];
    size_t dimz = p->invTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;
    int iterations;

    double tmpMat[9], tmpVec[3];
    double *to_voxel = p->the_to_voxel->m;
    double *to_real = p->inv_to_real->m;
    double conv[16];

    bal_imageGeometry geometry;

    double x, y;
    double dx, dy;
    double pe = 0, e = 0;

    float pInvX, pInvY;

    p->ndivergence = 0;
    p->nnonconvergence = 0;

    switch( p->invGeometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
        switch( p->theGeometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          geometry = _BAL_HOMOTHETY_GEOMETRY_;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          geometry = _BAL_TRANSLATION_GEOMETRY_;
          break;
        case _BAL_QFORM_GEOMETRY_ :
          geometry = _BAL_QFORM_GEOMETRY_;
          break;
        }
        break;
    case _BAL_TRANSLATION_GEOMETRY_ :
        switch( p->theGeometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
        case _BAL_TRANSLATION_GEOMETRY_ :
          geometry = _BAL_TRANSLATION_GEOMETRY_;
          break;
        case _BAL_QFORM_GEOMETRY_ :
          geometry = _BAL_QFORM_GEOMETRY_;
          break;
        }
        break;
    case _BAL_QFORM_GEOMETRY_ :
        geometry = _BAL_QFORM_GEOMETRY_;
        break;
    }

    /* \mathbf{H}_{ref,\mathbb{Z} \leftarrow \mathbb{R}}
       \times \mathbf{H}_{flo,\mathbb{R} \leftarrow \mathbb{Z}}
       = to_voxel * to_real
                       (  0  1  2  3 )
                       (  4  5  6  7 )
                       (  8  9 10 11 )
                       ( 12 13 14 15 )
       (  0  1  2  3 )
       (  4  5  6  7 )
       (  8  9 10 11 )
       ( 12 13 14 15 )
     */
    for ( i=0; i<16; i++ ) conv[i] = 0;
    conv[0] = conv[5] = conv[10] = conv[15] = 1.0;
    switch ( geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0];
      conv[ 5] = to_voxel[ 5] * to_real[ 5];
      conv[10] = to_voxel[10] * to_real[10];
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0];
      conv[ 3] = to_voxel[ 0] * to_real[ 3] + to_voxel[ 3];
      conv[ 5] = to_voxel[ 5] * to_real[ 5];
      conv[ 7] = to_voxel[ 5] * to_real[ 7] + to_voxel[ 7];
      conv[10] = to_voxel[10] * to_real[10];
      conv[11] = to_voxel[10] * to_real[11] + to_voxel[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0]
               + to_voxel[ 1] * to_real[ 4]
               + to_voxel[ 2] * to_real[ 8];
      conv[ 1] = to_voxel[ 0] * to_real[ 1]
               + to_voxel[ 1] * to_real[ 5]
               + to_voxel[ 2] * to_real[ 9];
      conv[ 2] = to_voxel[ 0] * to_real[ 2]
               + to_voxel[ 1] * to_real[ 6]
               + to_voxel[ 2] * to_real[10];
      conv[ 3] = to_voxel[ 0] * to_real[ 3]
               + to_voxel[ 1] * to_real[ 7]
               + to_voxel[ 2] * to_real[11]
               + to_voxel[ 3];
      conv[ 4] = to_voxel[ 4] * to_real[ 0]
               + to_voxel[ 5] * to_real[ 4]
               + to_voxel[ 6] * to_real[ 8];
      conv[ 5] = to_voxel[ 4] * to_real[ 1]
               + to_voxel[ 5] * to_real[ 5]
               + to_voxel[ 6] * to_real[ 9];
      conv[ 6] = to_voxel[ 4] * to_real[ 2]
               + to_voxel[ 5] * to_real[ 6]
               + to_voxel[ 6] * to_real[10];
      conv[ 7] = to_voxel[ 4] * to_real[ 3]
               + to_voxel[ 5] * to_real[ 7]
               + to_voxel[ 6] * to_real[11]
               + to_voxel[ 7];
      conv[ 8] = to_voxel[ 8] * to_real[ 0]
               + to_voxel[ 9] * to_real[ 4]
               + to_voxel[10] * to_real[ 8];
      conv[ 9] = to_voxel[ 8] * to_real[ 1]
               + to_voxel[ 9] * to_real[ 5]
               + to_voxel[10] * to_real[ 9];
      conv[10] = to_voxel[ 8] * to_real[ 2]
               + to_voxel[ 9] * to_real[ 6]
               + to_voxel[10] * to_real[10];
      conv[11] = to_voxel[ 8] * to_real[ 3]
               + to_voxel[ 9] * to_real[ 7]
               + to_voxel[10] * to_real[11]
               + to_voxel[11];
      break;
    }

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                e = pe = 0.0;

                for ( iterations = 0; iterations< ITERMAX; iterations++ ) {

                  /* M + I(M) in voxel
                   */
                  switch ( geometry ) {
                  default :
                  case _BAL_UNKNOWN_GEOMETRY_ :
                  case _BAL_HOMOTHETY_GEOMETRY_ :
                    x = conv[ 0] * i;
                    y = conv[ 5] * j;
                    break;
                  case _BAL_TRANSLATION_GEOMETRY_ :
                    x = conv[ 0] * i + conv[ 3];
                    y = conv[ 5] * j + conv[ 7];
                    break;
                  case _BAL_QFORM_GEOMETRY_ :
                    x = conv[ 0] * i + conv[ 1] * j + conv[ 2] * k + conv[ 3];
                    y = conv[ 4] * i + conv[ 5] * j + conv[ 6] * k + conv[ 7];
                    break;
                  }
                  switch ( p->theGeometry ) {
                  default :
                  case _BAL_UNKNOWN_GEOMETRY_ :
                  case _BAL_HOMOTHETY_GEOMETRY_ :
                  case _BAL_TRANSLATION_GEOMETRY_ :
                    x += to_voxel[ 0] * arrInvX[k][j][i];
                    y += to_voxel[ 5] * arrInvY[k][j][i];
                    break;
                  case _BAL_QFORM_GEOMETRY_ :
                    x += to_voxel[ 0] * arrInvX[k][j][i]
                       + to_voxel[ 1] * arrInvY[k][j][i];
                    y += to_voxel[ 4] * arrInvX[k][j][i]
                       + to_voxel[ 5] * arrInvY[k][j][i];
                    break;
                  }

                  /* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
                     and the vector v at M + I(M)
                  */

                  _get2DMatrixVector( x, y, k,
                                      arrXdX, arrXdY,
                                      arrYdX, arrYdY,
                                      arrTrsfX, arrTrsfY,
                                      p->theTrsfDim, tmpMat, tmpVec );

                  /* I(M) + v(M+I(M))
                   */
                  x = arrInvX[k][j][i] + tmpVec[0];
                  y = arrInvY[k][j][i] + tmpVec[1];

                  /* first tests on |I(M) + v(M+I(M))|
                   * - is it small enough?
                   * - does it increase ? -> retrieve previous values
                   */
                  e = fabs( x ) + fabs( y );
                  if ( e <= ERRMAX )
                    break;
                  if ( iterations > 0 && e > pe ) {
                      p->ndivergence ++;
                      arrInvX[k][j][i] = pInvX;
                      arrInvY[k][j][i] = pInvY;
                      break;
                  }
                  pe = e;

                  /* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
                     is the optimal variation of I(M)
                  */
                  dx = tmpMat[0] * x + tmpMat[1] * y;
                  dy = tmpMat[3] * x + tmpMat[4] * y;

                  pInvX = arrInvX[k][j][i];
                  pInvY = arrInvY[k][j][i];

                  arrInvX[k][j][i] -= alpha * dx;
                  arrInvY[k][j][i] -= alpha * dy;

                }

                if ( iterations >= ITERMAX ) {
                    p->nnonconvergence ++;
                }

                if ( iterations >= ITERMAX || (iterations > 0 && e > pe) ) {
                    if ( theErrors != (u8***)NULL )
                        theErrors[k][j][i] = 255;
                }

            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}





static void *_Inverse3DRealVectorField( void *par )
{
    char *proc = "_Inverse3DRealVectorField";
    typeChunk *chunk = (typeChunk *)par;
    void *parameter = chunk->parameters;
    size_t first = chunk->first;
    size_t last = chunk->last;

    _TransformationInversionParam *p = (_TransformationInversionParam*)parameter;

    float ***arrXdX = p->arrXdX;
    float ***arrXdY = p->arrXdY;
    float ***arrXdZ = p->arrXdZ;
    float ***arrYdX = p->arrYdX;
    float ***arrYdY = p->arrYdY;
    float ***arrYdZ = p->arrYdZ;
    float ***arrZdX = p->arrZdX;
    float ***arrZdY = p->arrZdY;
    float ***arrZdZ = p->arrZdZ;

    float ***arrInvX = p->arrInvX;
    float ***arrInvY = p->arrInvY;
    float ***arrInvZ = p->arrInvZ;

    float ***arrTrsfX = p->arrTrsfX;
    float ***arrTrsfY = p->arrTrsfY;
    float ***arrTrsfZ = p->arrTrsfZ;

    u8 ***theErrors = p->theErrors;

    size_t dimx = p->invTrsfDim[0];
    size_t dimy = p->invTrsfDim[1];
    size_t dimz = p->invTrsfDim[2];
    size_t dimxy = dimx * dimy;

    size_t i, j, k;
    size_t ifirst, jfirst, kfirst;
    size_t ilast, jlast, klast;
    size_t iend, jend;
    int iterations;

    double tmpMat[9], tmpVec[3];
    double *to_voxel = p->the_to_voxel->m;
    double *to_real = p->inv_to_real->m;
    double conv[16];

    bal_imageGeometry geometry;

    double x, y, z;
    double dx, dy, dz;
    double pe = 0, e = 0;

    float pInvX, pInvY, pInvZ;

    p->ndivergence = 0;
    p->nnonconvergence = 0;

    switch( p->invGeometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
        switch( p->theGeometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          geometry = _BAL_HOMOTHETY_GEOMETRY_;
          break;
        case _BAL_TRANSLATION_GEOMETRY_ :
          geometry = _BAL_TRANSLATION_GEOMETRY_;
          break;
        case _BAL_QFORM_GEOMETRY_ :
          geometry = _BAL_QFORM_GEOMETRY_;
          break;
        }
        break;
    case _BAL_TRANSLATION_GEOMETRY_ :
        switch( p->theGeometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
        case _BAL_TRANSLATION_GEOMETRY_ :
          geometry = _BAL_TRANSLATION_GEOMETRY_;
          break;
        case _BAL_QFORM_GEOMETRY_ :
          geometry = _BAL_QFORM_GEOMETRY_;
          break;
        }
        break;
    case _BAL_QFORM_GEOMETRY_ :
        geometry = _BAL_QFORM_GEOMETRY_;
        break;
    }

    /* \mathbf{H}_{ref,\mathbb{Z} \leftarrow \mathbb{R}}
       \times \mathbf{H}_{flo,\mathbb{R} \leftarrow \mathbb{Z}}
       = to_voxel * to_real
                       (  0  1  2  3 )
                       (  4  5  6  7 )
                       (  8  9 10 11 )
                       ( 12 13 14 15 )
       (  0  1  2  3 )
       (  4  5  6  7 )
       (  8  9 10 11 )
       ( 12 13 14 15 )
     */
    for ( i=0; i<16; i++ ) conv[i] = 0;
    conv[0] = conv[5] = conv[10] = conv[15] = 1.0;
    switch ( geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0];
      conv[ 5] = to_voxel[ 5] * to_real[ 5];
      conv[10] = to_voxel[10] * to_real[10];
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0];
      conv[ 3] = to_voxel[ 0] * to_real[ 3] + to_voxel[ 3];
      conv[ 5] = to_voxel[ 5] * to_real[ 5];
      conv[ 7] = to_voxel[ 5] * to_real[ 7] + to_voxel[ 7];
      conv[10] = to_voxel[10] * to_real[10];
      conv[11] = to_voxel[10] * to_real[11] + to_voxel[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      conv[ 0] = to_voxel[ 0] * to_real[ 0]
               + to_voxel[ 1] * to_real[ 4]
               + to_voxel[ 2] * to_real[ 8];
      conv[ 1] = to_voxel[ 0] * to_real[ 1]
               + to_voxel[ 1] * to_real[ 5]
               + to_voxel[ 2] * to_real[ 9];
      conv[ 2] = to_voxel[ 0] * to_real[ 2]
               + to_voxel[ 1] * to_real[ 6]
               + to_voxel[ 2] * to_real[10];
      conv[ 3] = to_voxel[ 0] * to_real[ 3]
               + to_voxel[ 1] * to_real[ 7]
               + to_voxel[ 2] * to_real[11]
               + to_voxel[ 3];
      conv[ 4] = to_voxel[ 4] * to_real[ 0]
               + to_voxel[ 5] * to_real[ 4]
               + to_voxel[ 6] * to_real[ 8];
      conv[ 5] = to_voxel[ 4] * to_real[ 1]
               + to_voxel[ 5] * to_real[ 5]
               + to_voxel[ 6] * to_real[ 9];
      conv[ 6] = to_voxel[ 4] * to_real[ 2]
               + to_voxel[ 5] * to_real[ 6]
               + to_voxel[ 6] * to_real[10];
      conv[ 7] = to_voxel[ 4] * to_real[ 3]
               + to_voxel[ 5] * to_real[ 7]
               + to_voxel[ 6] * to_real[11]
               + to_voxel[ 7];
      conv[ 8] = to_voxel[ 8] * to_real[ 0]
               + to_voxel[ 9] * to_real[ 4]
               + to_voxel[10] * to_real[ 8];
      conv[ 9] = to_voxel[ 8] * to_real[ 1]
               + to_voxel[ 9] * to_real[ 5]
               + to_voxel[10] * to_real[ 9];
      conv[10] = to_voxel[ 8] * to_real[ 2]
               + to_voxel[ 9] * to_real[ 6]
               + to_voxel[10] * to_real[10];
      conv[11] = to_voxel[ 8] * to_real[ 3]
               + to_voxel[ 9] * to_real[ 7]
               + to_voxel[10] * to_real[11]
               + to_voxel[11];
      break;
    }

    k = kfirst = first / (dimxy);
    j = jfirst = (first - kfirst*(dimxy)) / dimx;
    i = ifirst = (first - kfirst*(dimxy) - jfirst*dimx);

    klast = last / (dimxy);
    jlast = (last - klast*(dimxy)) / dimx;
    ilast = (last - klast*(dimxy) - jlast*dimx);

    for ( ; k<=klast; k++, j=0 ) {
        if ( _verbose_ >= 4 ) {
          fprintf( stderr, " ... processing slice #%3lu/%3lu in %s\n", k, dimz, proc );
        }
        jend = (k==klast) ? jlast+1 : dimy;
        for ( ; j<jend; j++, i=0 ) {
            iend = (j==jlast && k==klast) ? ilast+1 : dimx;
            for ( ; i<iend; i++ ) {

                e = pe = 0.0;

                for ( iterations = 0; iterations< ITERMAX; iterations++ ) {

                  /* M + I(M) in voxel
                   */
                  switch ( geometry ) {
                  default :
                  case _BAL_UNKNOWN_GEOMETRY_ :
                  case _BAL_HOMOTHETY_GEOMETRY_ :
                    x = conv[ 0] * i;
                    y = conv[ 5] * j;
                    z = conv[10] * k;
                    break;
                  case _BAL_TRANSLATION_GEOMETRY_ :
                    x = conv[ 0] * i + conv[ 3];
                    y = conv[ 5] * j + conv[ 7];
                    z = conv[10] * k + conv[11];
                    break;
                  case _BAL_QFORM_GEOMETRY_ :
                    x = conv[ 0] * i + conv[ 1] * j + conv[ 2] * k + conv[ 3];
                    y = conv[ 4] * i + conv[ 5] * j + conv[ 6] * k + conv[ 7];
                    z = conv[ 8] * i + conv[ 9] * j + conv[10] * k + conv[11];
                    break;
                  }
                  switch ( p->theGeometry ) {
                  default :
                  case _BAL_UNKNOWN_GEOMETRY_ :
                  case _BAL_HOMOTHETY_GEOMETRY_ :
                  case _BAL_TRANSLATION_GEOMETRY_ :
                    x += to_voxel[ 0] * arrInvX[k][j][i];
                    y += to_voxel[ 5] * arrInvY[k][j][i];
                    z += to_voxel[10] * arrInvZ[k][j][i];
                    break;
                  case _BAL_QFORM_GEOMETRY_ :
                    x += to_voxel[ 0] * arrInvX[k][j][i]
                       + to_voxel[ 1] * arrInvY[k][j][i]
                       + to_voxel[ 2] * arrInvZ[k][j][i];
                    y += to_voxel[ 4] * arrInvX[k][j][i]
                       + to_voxel[ 5] * arrInvY[k][j][i]
                       + to_voxel[ 6] * arrInvZ[k][j][i];
                    z += to_voxel[ 8] * arrInvX[k][j][i]
                       + to_voxel[ 9] * arrInvY[k][j][i]
                       + to_voxel[10] * arrInvZ[k][j][i];
                    break;
                  }

                  /* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
                     and the vector v at M + I(M)
                  */

                  _get3DMatrixVector( x, y, z,
                                      arrXdX, arrXdY, arrXdZ,
                                      arrYdX, arrYdY, arrYdZ,
                                      arrZdX, arrZdY, arrZdZ,
                                      arrTrsfX, arrTrsfY, arrTrsfZ,
                                      p->theTrsfDim, tmpMat, tmpVec );
                  /* I(M) + v(M+I(M))
                   */
                  x = arrInvX[k][j][i] + tmpVec[0];
                  y = arrInvY[k][j][i] + tmpVec[1];
                  z = arrInvZ[k][j][i] + tmpVec[2];

                  /* first tests on |I(M) + v(M+I(M))|
                   * - is it small enough?
                   * - does it increase ? -> retrieve previous values
                   */
                  e = fabs( x ) + fabs( y ) + fabs( z );
                  if ( e <= ERRMAX )
                    break;
                  if ( iterations > 0 && e > pe ) {
                      p->ndivergence ++;
                      arrInvX[k][j][i] = pInvX;
                      arrInvY[k][j][i] = pInvY;
                      arrInvZ[k][j][i] = pInvZ;
                      break;
                  }
                  pe = e;

                  /* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
                     is the optimal variation of I(M)
                  */
                  dx = tmpMat[0] * x + tmpMat[1] * y + tmpMat[2] * z;
                  dy = tmpMat[3] * x + tmpMat[4] * y + tmpMat[5] * z;
                  dz = tmpMat[6] * x + tmpMat[7] * y + tmpMat[8] * z;

                  pInvX = arrInvX[k][j][i];
                  pInvY = arrInvY[k][j][i];
                  pInvZ = arrInvZ[k][j][i];

                  arrInvX[k][j][i] -= alpha * dx;
                  arrInvY[k][j][i] -= alpha * dy;
                  arrInvZ[k][j][i] -= alpha * dz;

                }

                if ( iterations >= ITERMAX ) {
                    p->nnonconvergence ++;
                }

                if ( iterations >= ITERMAX || (iterations > 0 && e > pe) ) {
                    if ( theErrors != (u8***)NULL )
                        theErrors[k][j][i] = 255;
                }

            }
        }
    }

    chunk->ret = 1;
    return( (void*)NULL );
}









/********************************************************************************
 *
 * vector field inversion initialization
 *
 ********************************************************************************/



#define _UPDATE_2DINITIALIZATION( _X_, _Y_, _Z_, _W_ ) { \
    arrInvX[_Z_][_Y_][_X_] -= (_W_) * arrTrsfX[k][j][i]; \
    arrInvY[_Z_][_Y_][_X_] -= (_W_) * arrTrsfY[k][j][i]; \
    arrWght[_Z_][_Y_][_X_] += (_W_);                     \
}

#define _UPDATE_3DINITIALIZATION( _X_, _Y_, _Z_, _W_ ) { \
    arrInvX[_Z_][_Y_][_X_] -= (_W_) * arrTrsfX[k][j][i]; \
    arrInvY[_Z_][_Y_][_X_] -= (_W_) * arrTrsfY[k][j][i]; \
    arrInvZ[_Z_][_Y_][_X_] -= (_W_) * arrTrsfZ[k][j][i]; \
    arrWght[_Z_][_Y_][_X_] += (_W_);                     \
}





static void _ForwardInterpolateInverse2DVoxelVectorField( bal_transformation *theTrsf,
                                                    bal_transformation *invTrsf,
                                                    bal_image *imWeight )
{
    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrTrsfX = (float***)theTrsf->vx.array;
    float ***arrTrsfY = (float***)theTrsf->vy.array;
    float ***arrWght = (float***)imWeight->array;
    size_t i, j, k;
    float x, y;
    int ix, iy;
    float dx, dy;

    for ( k=0; k<invTrsf->vx.nplanes; k++ )
    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ )
        arrInvY[k][j][i] = arrInvX[k][j][i] = arrWght[k][j][i] = 0.0;

    for ( k=0; k<theTrsf->vx.nplanes; k++ )
    for ( j=0; j<theTrsf->vx.nrows; j++ )
    for ( i=0; i<theTrsf->vx.ncols; i++ ) {
        x = i + arrTrsfX[k][j][i];
        y = j + arrTrsfY[k][j][i];

        if ( y < -1.0 ) continue;

        if ( y < 0.0 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, 0, k, 1 );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, 0, k, 1-dx );
                _UPDATE_2DINITIALIZATION( ix+1, 0, k, dx );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, 0, k, 1 );
                continue;
            }
            continue;
        }

        iy = (int)y;
        dy = y - iy;

        if ( iy < (int)invTrsf->vx.nrows-1 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, iy, k, 1-dy );
                _UPDATE_2DINITIALIZATION( 0, iy+1, k, dy );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, (1-dx) * (1-dy) );
                _UPDATE_2DINITIALIZATION( ix+1, iy, k, dx * (1-dy) );
                _UPDATE_2DINITIALIZATION( ix, iy+1, k, (1-dx) * dy );
                _UPDATE_2DINITIALIZATION( ix+1, iy+1, k, dx * dy );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1-dy );
                _UPDATE_2DINITIALIZATION( ix, iy+1, k, dy );
                continue;
            }
            continue;
        }

        if ( iy == (int)invTrsf->vx.nrows-1 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, iy, k, 1 );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1-dx );
                _UPDATE_2DINITIALIZATION( ix+1, iy, k, dx );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1 );
                continue;
            }
            continue;
        }
        continue;
    }
}





static void _ForwardInterpolateInverse3DVoxelVectorField( bal_transformation *theTrsf,
                                                    bal_transformation *invTrsf,
                                                    bal_image *imWeight )
{
    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrInvZ = (float***)invTrsf->vz.array;
    float ***arrTrsfX = (float***)theTrsf->vx.array;
    float ***arrTrsfY = (float***)theTrsf->vy.array;
    float ***arrTrsfZ = (float***)theTrsf->vz.array;
    float ***arrWght = (float***)imWeight->array;
    size_t i, j, k;
    float x, y, z;
    int ix, iy, iz;
    float dx, dy, dz;

    for ( k=0; k<invTrsf->vx.nplanes; k++ )
    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ )
        arrInvZ[k][j][i] = arrInvY[k][j][i] = arrInvX[k][j][i] = arrWght[k][j][i] = 0.0;

    for ( k=0; k<theTrsf->vx.nplanes; k++ )
    for ( j=0; j<theTrsf->vx.nrows; j++ )
    for ( i=0; i<theTrsf->vx.ncols; i++ ) {
        x = i + arrTrsfX[k][j][i];
        y = j + arrTrsfY[k][j][i];
        z = k + arrTrsfZ[k][j][i];

        if ( z < -1.0 ) continue;

        if ( z < 0.0 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, 0, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, 0, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, 0, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, 0, 1 );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, 0, 1-dy );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, 0, dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, (1-dx) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, 0, dx * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, 0, (1-dx) * dy );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, 0, dx * dy );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1-dy );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, 0, dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, 0, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, 0, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1 );
                    continue;
                }

                continue;
            }
            continue;
        }

        iz = (int)z;
        dz = z - iz;

        if ( iz < (int)invTrsf->vx.nplanes-1 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( 0, 0, iz+1, dz );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, (1-dz) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz, (1-dz) * dx );
                    _UPDATE_3DINITIALIZATION( ix, 0, iz+1, dz * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz+1, dz * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( ix, 0, iz+1, dz );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, (1-dz) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz, (1-dz) * dy );
                    _UPDATE_3DINITIALIZATION( 0, iy, iz+1, dz * (1-dy) );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz+1, dz * dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dy) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, (1-dz) * (1-dy) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dz) * dy * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz, (1-dz) * dy * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dy) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz+1, dz * (1-dy) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz+1, dz * dy * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz+1, dz * dy * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dz) * dy );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz+1, dz * dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( 0, iy, iz+1, dz );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, (1-dz) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz+1, dz * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz );
                    continue;
                }
                continue;
            }
            continue;
        }

        if ( iz == (int)invTrsf->vx.nplanes-1 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, iz, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1 );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1-dy );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz, dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dx) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, dx * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dx) * dy );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz, dx * dy );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dy );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1 );
                    continue;
                }
                continue;
            }
            continue;
        }
        continue;
    }
}




static void _ForwardInterpolateInverse2DRealVectorField( bal_transformation *theTrsf,
                                                    bal_transformation *invTrsf,
                                                    bal_image *imWeight )
{
    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrTrsfX = (float***)theTrsf->vx.array;
    float ***arrTrsfY = (float***)theTrsf->vy.array;
    float ***arrWght = (float***)imWeight->array;
    size_t i, j, k;
    float u, v;
    float x, y;
    int ix, iy;
    float dx, dy;

    double *to_real  = theTrsf->vx.to_real.m;
    double *to_voxel = invTrsf->vx.to_voxel.m;

    for ( k=0; k<invTrsf->vx.nplanes; k++ )
    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ )
        arrInvY[k][j][i] = arrInvX[k][j][i] = arrWght[k][j][i] = 0.0;

    for ( k=0; k<theTrsf->vx.nplanes; k++ )
    for ( j=0; j<theTrsf->vx.nrows; j++ )
    for ( i=0; i<theTrsf->vx.ncols; i++ ) {

        /* (  0  1  2  3 )
         * (  4  5  6  7 )
         * (  8  9 10 11 )
         * ( 12 13 14 15 )
         */
        switch ( theTrsf->vx.geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          u = to_real[ 0] * i;
          v = to_real[ 5] * j;
        case _BAL_TRANSLATION_GEOMETRY_ :
          u = to_real[ 0] * i + to_real[ 3];
          v = to_real[ 5] * j + to_real[ 7];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          u = to_real[ 0] * i + to_real[ 1] * j + to_real[ 2] * k + to_real[ 3];
          v = to_real[ 4] * i + to_real[ 5] * j + to_real[ 6] * k + to_real[ 7];
          break;
        }

        u += arrTrsfX[k][j][i];
        v += arrTrsfY[k][j][i];

        switch ( invTrsf->vx.geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          x = to_voxel[ 0] * u;
          y = to_voxel[ 5] * v;
        case _BAL_TRANSLATION_GEOMETRY_ :
          x = to_voxel[ 0] * u + to_voxel[ 3];
          y = to_voxel[ 5] * v + to_voxel[ 7];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          x = to_voxel[ 0] * u + to_voxel[ 1] * v + to_voxel[ 3];
          y = to_voxel[ 4] * u + to_voxel[ 5] * v + to_voxel[ 7];
          break;
        }

        if ( y < -1.0 ) continue;

        if ( y < 0.0 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, 0, k, 1 );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, 0, k, 1-dx );
                _UPDATE_2DINITIALIZATION( ix+1, 0, k, dx );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, 0, k, 1 );
                continue;
            }
            continue;
        }

        iy = (int)y;
        dy = y - iy;

        if ( iy < (int)invTrsf->vx.nrows-1 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, iy, k, 1-dy );
                _UPDATE_2DINITIALIZATION( 0, iy+1, k, dy );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, (1-dx) * (1-dy) );
                _UPDATE_2DINITIALIZATION( ix+1, iy, k, dx * (1-dy) );
                _UPDATE_2DINITIALIZATION( ix, iy+1, k, (1-dx) * dy );
                _UPDATE_2DINITIALIZATION( ix+1, iy+1, k, dx * dy );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1-dy );
                _UPDATE_2DINITIALIZATION( ix, iy+1, k, dy );
                continue;
            }
            continue;
        }

        if ( iy == (int)invTrsf->vx.nrows-1 ) {
            if ( x < -1.0 ) continue;
            if ( x < 0.0 ) {
                _UPDATE_2DINITIALIZATION( 0, iy, k, 1 );
                continue;
            }
            ix = (int)x;
            dx = x - ix;
            if ( ix < (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1-dx );
                _UPDATE_2DINITIALIZATION( ix+1, iy, k, dx );
                continue;
            }
            if ( ix == (int)invTrsf->vx.ncols-1 ) {
                _UPDATE_2DINITIALIZATION( ix, iy, k, 1 );
                continue;
            }
            continue;
        }
        continue;
    }
}





static void _ForwardInterpolateInverse3DRealVectorField( bal_transformation *theTrsf,
                                                    bal_transformation *invTrsf,
                                                    bal_image *imWeight )
{
    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrInvZ = (float***)invTrsf->vz.array;
    float ***arrTrsfX = (float***)theTrsf->vx.array;
    float ***arrTrsfY = (float***)theTrsf->vy.array;
    float ***arrTrsfZ = (float***)theTrsf->vz.array;
    float ***arrWght = (float***)imWeight->array;
    size_t i, j, k;
    float u, v, w;
    float x, y, z;
    int ix, iy, iz;
    float dx, dy, dz;

    double *to_real  = theTrsf->vx.to_real.m;
    double *to_voxel = invTrsf->vx.to_voxel.m;

    for ( k=0; k<invTrsf->vx.nplanes; k++ )
    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ )
        arrInvZ[k][j][i] = arrInvY[k][j][i] = arrInvX[k][j][i] = arrWght[k][j][i] = 0.0;

    for ( k=0; k<theTrsf->vx.nplanes; k++ )
    for ( j=0; j<theTrsf->vx.nrows; j++ )
    for ( i=0; i<theTrsf->vx.ncols; i++ ) {

        /* (  0  1  2  3 )
         * (  4  5  6  7 )
         * (  8  9 10 11 )
         * ( 12 13 14 15 )
         */
        switch ( theTrsf->vx.geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          u = to_real[ 0] * i;
          v = to_real[ 5] * j;
          w = to_real[10] * k;
        case _BAL_TRANSLATION_GEOMETRY_ :
          u = to_real[ 0] * i + to_real[ 3];
          v = to_real[ 5] * j + to_real[ 7];
          w = to_real[10] * k + to_real[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          u = to_real[ 0] * i + to_real[ 1] * j + to_real[ 2] * k + to_real[ 3];
          v = to_real[ 4] * i + to_real[ 5] * j + to_real[ 6] * k + to_real[ 7];
          w = to_real[ 8] * i + to_real[ 9] * j + to_real[10] * k + to_real[11];
          break;
        }

        u += arrTrsfX[k][j][i];
        v += arrTrsfY[k][j][i];
        w += arrTrsfZ[k][j][i];

        switch ( invTrsf->vx.geometry ) {
        default :
        case _BAL_UNKNOWN_GEOMETRY_ :
        case _BAL_HOMOTHETY_GEOMETRY_ :
          x = to_voxel[ 0] * u;
          y = to_voxel[ 5] * v;
          z = to_voxel[10] * w;
        case _BAL_TRANSLATION_GEOMETRY_ :
          x = to_voxel[ 0] * u + to_voxel[ 3];
          y = to_voxel[ 5] * v + to_voxel[ 7];
          z = to_voxel[10] * w + to_voxel[11];
          break;
        case _BAL_QFORM_GEOMETRY_ :
          x = to_voxel[ 0] * u + to_voxel[ 1] * v + to_voxel[ 2] * w + to_voxel[ 3];
          y = to_voxel[ 4] * u + to_voxel[ 5] * v + to_voxel[ 6] * w + to_voxel[ 7];
          z = to_voxel[ 8] * u + to_voxel[ 9] * v + to_voxel[10] * w + to_voxel[11];
          break;
        }

        if ( z < -1.0 ) continue;

        if ( z < 0.0 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, 0, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, 0, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, 0, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, 0, 1 );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, 0, 1-dy );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, 0, dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, (1-dx) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, 0, dx * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, 0, (1-dx) * dy );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, 0, dx * dy );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1-dy );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, 0, dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, 0, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, 0, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, 0, 1 );
                    continue;
                }

                continue;
            }
            continue;
        }

        iz = (int)z;
        dz = z - iz;

        if ( iz < (int)invTrsf->vx.nplanes-1 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( 0, 0, iz+1, dz );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, (1-dz) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz, (1-dz) * dx );
                    _UPDATE_3DINITIALIZATION( ix, 0, iz+1, dz * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz+1, dz * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( ix, 0, iz+1, dz );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, (1-dz) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz, (1-dz) * dy );
                    _UPDATE_3DINITIALIZATION( 0, iy, iz+1, dz * (1-dy) );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz+1, dz * dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dy) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, (1-dz) * (1-dy) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dz) * dy * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz, (1-dz) * dy * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dy) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz+1, dz * (1-dy) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz+1, dz * dy * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz+1, dz * dy * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dz) * dy );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz+1, dz * dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( 0, iy, iz+1, dz );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dz) * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, (1-dz) * dx );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz * (1-dx) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz+1, dz * dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dz );
                    _UPDATE_3DINITIALIZATION( ix, iy, iz+1, dz );
                    continue;
                }
                continue;
            }
            continue;
        }

        if ( iz == (int)invTrsf->vx.nplanes-1 ) {
            if ( y < -1.0 ) continue;
            if ( y < 0.0 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, 0, iz, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, 0, iz, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, 0, iz, 1 );
                    continue;
                }
                continue;
            }
            iy = (int)y;
            dy = y - iy;
            if ( iy < (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1-dy );
                    _UPDATE_3DINITIALIZATION( 0, iy+1, iz, dy );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, (1-dx) * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, dx * (1-dy) );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, (1-dx) * dy );
                    _UPDATE_3DINITIALIZATION( ix+1, iy+1, iz, dx * dy );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dy );
                    _UPDATE_3DINITIALIZATION( ix, iy+1, iz, dy );
                    continue;
                }
                continue;
            }
            if ( iy == (int)invTrsf->vx.nrows-1 ) {
                if ( x < -1.0 ) continue;
                if ( x < 0.0 ) {
                    _UPDATE_3DINITIALIZATION( 0, iy, iz, 1 );
                    continue;
                }
                ix = (int)x;
                dx = x - ix;
                if ( ix < (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1-dx );
                    _UPDATE_3DINITIALIZATION( ix+1, iy, iz, dx );
                    continue;
                }
                if ( ix == (int)invTrsf->vx.ncols-1 ) {
                    _UPDATE_3DINITIALIZATION( ix, iy, iz, 1 );
                    continue;
                }
                continue;
            }
            continue;
        }
        continue;
    }
}





/********************************************************************************
 *
 * vector field inversion initialization
 *
 ********************************************************************************/



static int BAL_Inverse2DVectorFieldInitialization( bal_transformation *theTrsf,
                                               bal_transformation *invTrsf )
{
    char *proc = "BAL_Inverse2DVectorFieldInitialization";

    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrWght;
    bal_image imWeight;
    size_t i, j, k;
    bal_doublePoint theSigma;

    switch( initialization ) {
    default :
        if ( _verbose_ )
            fprintf( stderr, "%s: such case not handled yet\n", proc );
        return( -1 );
    case ZERO :
        for ( k=0; k<invTrsf->vx.nplanes; k++ )
        for ( j=0; j<invTrsf->vx.nrows; j++ )
        for ( i=0; i<invTrsf->vx.ncols; i++ )
            arrInvY[k][j][i] = arrInvX[k][j][i] = 0.0;
        break;
    case FORWARD_INTERPOLATION :
        if ( BAL_AllocScalarImageFromImage( &imWeight, (char*)NULL, &(invTrsf->vx), FLOAT ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to allocate weight image\n", proc );
            return( -1 );
        }

        switch( theTrsf->transformation_unit ) {
        default :
        case UNDEF_UNIT :
          BAL_FreeImage( &imWeight );
          if ( _verbose_ )
              fprintf( stderr, "%s: unknown transformation unit\n", proc );
          return( -1 );
         case VOXEL_UNIT :
          _ForwardInterpolateInverse2DVoxelVectorField( theTrsf, invTrsf, &imWeight );
          break;
        case REAL_UNIT :
          _ForwardInterpolateInverse2DRealVectorField( theTrsf, invTrsf, &imWeight );
          break;
        }

        theSigma.x = theSigma.y = theSigma.z = forwardSigma;

        if ( BAL_SmoothImage( &(invTrsf->vx), &theSigma ) != 1
             || BAL_SmoothImage( &(invTrsf->vy), &theSigma ) != 1
             || BAL_SmoothImage( &imWeight , &theSigma ) != 1 ) {
            BAL_FreeImage( &imWeight );
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to smooth images\n", proc );
            return( -1 );
        }

        arrWght = (float***)imWeight.array;
        for ( k=0; k<invTrsf->vx.nplanes; k++ )
        for ( j=0; j<invTrsf->vx.nrows; j++ )
        for ( i=0; i<invTrsf->vx.ncols; i++ ) {
            if ( arrWght[k][j][i] > 0.01 ) {
                arrInvX[k][j][i] /= arrWght[k][j][i];
                arrInvY[k][j][i] /= arrWght[k][j][i];
            }
        }
        BAL_FreeImage( &imWeight );

    }
    return( 1 );
}





static int BAL_Inverse3DVectorFieldInitialization( bal_transformation *theTrsf,
                                               bal_transformation *invTrsf )
{
    char *proc = "BAL_Inverse3DVectorFieldInitialization";

    float ***arrInvX = (float***)invTrsf->vx.array;
    float ***arrInvY = (float***)invTrsf->vy.array;
    float ***arrInvZ = (float***)invTrsf->vz.array;
    float ***arrWght;
    bal_image imWeight;
    size_t i, j, k;
    bal_doublePoint theSigma;

    switch( initialization ) {
    default :
        if ( _verbose_ )
            fprintf( stderr, "%s: such case not handled yet\n", proc );
        return( -1 );
    case ZERO :
        for ( k=0; k<invTrsf->vx.nplanes; k++ )
        for ( j=0; j<invTrsf->vx.nrows; j++ )
        for ( i=0; i<invTrsf->vx.ncols; i++ )
            arrInvZ[k][j][i] = arrInvY[k][j][i] = arrInvX[k][j][i] = 0.0;
        break;
    case FORWARD_INTERPOLATION :
        if ( BAL_AllocScalarImageFromImage( &imWeight, (char*)NULL, &(invTrsf->vx), FLOAT ) != 1 ) {
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to allocate weight image\n", proc );
            return( -1 );
        }

        switch( theTrsf->transformation_unit ) {
        default :
        case UNDEF_UNIT :
          BAL_FreeImage( &imWeight );
          if ( _verbose_ )
              fprintf( stderr, "%s: unknown transformation unit\n", proc );
          return( -1 );
         case VOXEL_UNIT :
          _ForwardInterpolateInverse3DVoxelVectorField( theTrsf, invTrsf, &imWeight );
          break;
        case REAL_UNIT :
          _ForwardInterpolateInverse3DRealVectorField( theTrsf, invTrsf, &imWeight );
          break;
        }

        theSigma.x = theSigma.y = theSigma.z = forwardSigma;

        if ( BAL_SmoothImage( &(invTrsf->vx), &theSigma ) != 1
             || BAL_SmoothImage( &(invTrsf->vy), &theSigma ) != 1
             || BAL_SmoothImage( &(invTrsf->vz), &theSigma ) != 1
             || BAL_SmoothImage( &imWeight , &theSigma ) != 1 ) {
            BAL_FreeImage( &imWeight );
            if ( _verbose_ )
                fprintf( stderr, "%s: unable to smooth images\n", proc );
            return( -1 );
        }

        arrWght = (float***)imWeight.array;
        for ( k=0; k<invTrsf->vx.nplanes; k++ )
        for ( j=0; j<invTrsf->vx.nrows; j++ )
        for ( i=0; i<invTrsf->vx.ncols; i++ ) {
            if ( arrWght[k][j][i] > 0.01 ) {
                arrInvX[k][j][i] /= arrWght[k][j][i];
                arrInvY[k][j][i] /= arrWght[k][j][i];
                arrInvZ[k][j][i] /= arrWght[k][j][i];
            }
        }
        BAL_FreeImage( &imWeight );

    }
    return( 1 );
}





/* On utilise l'algorithme de Newton, comme P. Cachier

   principe: soit I le champ inverse et V le champ direct,
   on a M' = M + V(M)
   M = M' + I(M')
   soit M' = M' + I(M') + V(M' +I(M'))
   donc  I(M) + V(M +I(M)) = 0

   On cherche a minimiser f =  I(M) + V(M +I(M))

   on ajoute une variation delta (d) a I(M) on a
   f = I(M) + d + V(M + I(M) + d)
   V(M'+d) = V(M') + (V.N^t)(M') d avec N l'operateur nabla
   donc f = I(M) + d + V(M + I(M)) +  (V.N^t)(M + I(M)) d
   qui s'annule pour
   d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )

   Question: faut-il faire une adaptation selon les tailles de pixel/voxel ?

*/





static int BAL_Inverse2DVectorField( bal_transformation *theTrsf,
                                     bal_transformation *invTrsf )
{
  char *proc = "BAL_Inverse2DVectorField";

  bal_image theXdX, theXdY, theYdX, theYdY;
  bal_doublePoint theSigma;
  u8 ***theErrors = (u8 ***)NULL;
  size_t i, j, k;

  int ndivergence = 0;
  int nnonconvergence = 0;

  typeChunks chunks;
  size_t first = 0;
  size_t last;
  _TransformationInversionParam *p = (_TransformationInversionParam*)NULL;
  int n;



  if ( _debug_ ) fprintf( stderr, "%s: entry\n", proc );



  /* check parameters
   */

  theSigma.x = theSigma.y = theSigma.z = derivationSigma;

  if ( invTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: input and inverse transformation should have the same z dimension\n", proc );
      return( -1 );
  }

  if ( imInverseErrors != (bal_image*)NULL ) {
      if ( imInverseErrors->ncols != invTrsf->vx.ncols || imInverseErrors->nrows != invTrsf->vx.nrows
           || imInverseErrors->nplanes != invTrsf->vx.nplanes || imInverseErrors->vdim != invTrsf->vx.vdim ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: error image and inverse transformation should have the same dimensions\n", proc );
      }
      else if ( imInverseErrors->type != UCHAR )
      {
          if ( _verbose_ )
              fprintf( stderr, "%s: error image should be of 'unsigned char' type\n", proc );
      }
      else {
          theErrors = (u8***)(imInverseErrors->array);
          for ( k=0; k<imInverseErrors->nplanes; k++ )
          for ( j=0; j<imInverseErrors->nrows; j++ )
          for ( i=0; i<imInverseErrors->ncols; i++ ) {
            theErrors[k][j][i] = 0;
          }
      }
  }



  /* allocation of derivative images
   */

  if ( BAL_AllocScalarImageFromImage( &theXdX, "XderivativeOfX.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of X component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theXdY, "YderivativeOfX.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of X component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theYdX, "XderivativeOfY.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Y component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theYdY, "YderivativeOfY.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Y component\n", proc );
    return( -1 );
  }





  /* calculation of derivatives
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: derivatives computation\n", proc );
  }


  if ( BAL_2DDerivativesOfImage( &(theTrsf->vx),
                                 &theXdX, &theXdY,
                                 &theSigma ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of X component\n", proc );
    return( -1 );
  }

  if ( BAL_2DDerivativesOfImage( &(theTrsf->vy),
                                 &theYdX, &theYdY,
                                 &theSigma ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Y component\n", proc );
    return( -1 );
  }



  /* add Identity to V.Nabla^t
     and invert (Id + V.Nabla^t)
     = ( xdx + 1     xdy     )
       ( ydx         ydy + 1 )

     The inverse is computed with
     | a11 a12 |-1             |  a22 -a12 |
     | a21 a22 |    =  1/DET * | -a21  a11 |

     with DET  =  a11 a22- a12 a21

   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: matrices computation\n", proc );
  }

  /* preparing parallelism for matrix inversion ( (Id + V.Nabla^t) )
   * arrXdX, arrYdY are in the geometry of theTrsf
   */

  first = 0;
  last = (theTrsf->vx).nplanes * (theTrsf->vx).nrows * (theTrsf->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
  }

  p = (_TransformationInversionParam*)vtmalloc( chunks.n_allocated_chunks * sizeof(_TransformationInversionParam),
                                                "p", proc );
  if ( p == (_TransformationInversionParam*)NULL ) {
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
      return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
      _initTransformationInversionParam( &(p[n]) );

      p[n].theTrsfDim[0] = theTrsf->vx.ncols;
      p[n].theTrsfDim[1] = theTrsf->vx.nrows;
      p[n].theTrsfDim[2] = theTrsf->vx.nplanes;

      p[n].theGeometry = theTrsf->vx.geometry;
      p[n].invGeometry = invTrsf->vx.geometry;

      p[n].the_to_voxel = &(theTrsf->vx.to_voxel);
      p[n].inv_to_real  = &(invTrsf->vx.to_real);

      p[n].arrXdX = (float***)(theXdX.array);
      p[n].arrXdY = (float***)(theXdY.array);
      p[n].arrYdX = (float***)(theYdX.array);
      p[n].arrYdY = (float***)(theYdY.array);

      chunks.data[n].parameters = (void*)(&(p[n]));
  }

  switch ( theTrsf->transformation_unit ) {
  default :
  case UNDEF_UNIT :
    if ( _verbose_ )
        fprintf( stderr, "%s: unknown transformation unit (2D case, matrices inversion)\n", proc );
    vtfree( p );
    freeChunks( &chunks );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    return( -1 );
  case VOXEL_UNIT :
    if ( processChunks( &_Inverse2DVoxelMatrix, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert matrices (2D case, voxel)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  case REAL_UNIT :
    if ( processChunks( &_Inverse2DRealMatrix, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert matrices (2D case, real)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  }

  vtfree( p );
  freeChunks( &chunks );



  /* initialization of inverse vector field
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: initialization of inverse vector field\n", proc );
  }

  if ( BAL_Inverse2DVectorFieldInitialization( theTrsf, invTrsf ) != 1 ) {
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
     if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize inverse vector field\n", proc );
      return( -1 );
  }



  /* preparing parallelism for vector field inversion
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: vectors inversion\n", proc );
  }

  first = 0;
  last = (invTrsf->vx).nplanes * (invTrsf->vx).nrows * (invTrsf->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
     if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
  }

  p = (_TransformationInversionParam*)vtmalloc( chunks.n_allocated_chunks * sizeof(_TransformationInversionParam),
                                                "p", proc );
  if ( p == (_TransformationInversionParam*)NULL ) {
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
      return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
      _initTransformationInversionParam( &(p[n]) );

      p[n].theTrsfDim[0] = theTrsf->vx.ncols;
      p[n].theTrsfDim[1] = theTrsf->vx.nrows;
      p[n].theTrsfDim[2] = theTrsf->vx.nplanes;

      p[n].invTrsfDim[0] = invTrsf->vx.ncols;
      p[n].invTrsfDim[1] = invTrsf->vx.nrows;
      p[n].invTrsfDim[2] = invTrsf->vx.nplanes;

      p[n].theGeometry = theTrsf->vx.geometry;
      p[n].invGeometry = invTrsf->vx.geometry;

      p[n].the_to_voxel = &(theTrsf->vx.to_voxel);
      p[n].inv_to_real  = &(invTrsf->vx.to_real);

      p[n].arrXdX = (float***)(theXdX.array);
      p[n].arrXdY = (float***)(theXdY.array);
      p[n].arrYdX = (float***)(theYdX.array);
      p[n].arrYdY = (float***)(theYdY.array);

      p[n].arrInvX = (float***)(invTrsf->vx.array);
      p[n].arrInvY = (float***)(invTrsf->vy.array);

      p[n].arrTrsfX = (float***)(theTrsf->vx.array);
      p[n].arrTrsfY = (float***)(theTrsf->vy.array);

      if ( theErrors != (u8***)NULL )
          p[n].theErrors = theErrors;

      p[n].ndivergence = 0;
      p[n].nnonconvergence = 0;

      chunks.data[n].parameters = (void*)(&(p[n]));
  }

  /* inversion
   */

  switch ( theTrsf->transformation_unit ) {
  default :
  case UNDEF_UNIT :
    if ( _verbose_ )
        fprintf( stderr, "%s: unknown transformation unit (2D case, matrices inversion)\n", proc );
    vtfree( p );
    freeChunks( &chunks );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    return( -1 );
  case VOXEL_UNIT :
    if ( processChunks( &_Inverse2DVoxelVectorField, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert vector field (2D case, voxel)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
     return( -1 );
    }
    break;
  case REAL_UNIT :
    if ( processChunks( &_Inverse2DRealVectorField, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert vector field (2D case, real)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
     return( -1 );
    }
    break;
  }

  for ( ndivergence=0, nnonconvergence=0, n=0; n<chunks.n_allocated_chunks; n++ ) {
      ndivergence += p[n].ndivergence;
      nnonconvergence += p[n].nnonconvergence;
  }

  vtfree( p );
  freeChunks( &chunks );



  if ( _verbose_ && (ndivergence > 0 || nnonconvergence > 0) ) {
    fprintf( stderr, "%s: divergence: %d/%lu, non-convergence: %d/%lu\n", proc,
             ndivergence,
             invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols,
             nnonconvergence,
             invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols );
  }

  BAL_FreeImage( &theYdY );
  BAL_FreeImage( &theYdX );
  BAL_FreeImage( &theXdY );
  BAL_FreeImage( &theXdX );

  return( 1 );

}





static int BAL_Inverse3DVectorField( bal_transformation *theTrsf,
                                     bal_transformation *invTrsf )
{
  char *proc = "BAL_Inverse3DVectorField";

  bal_image theXdX, theXdY, theXdZ, theYdX, theYdY, theYdZ, theZdX, theZdY, theZdZ;
  bal_doublePoint theSigma;
  u8 ***theErrors = (u8 ***)NULL;
  size_t i, j, k;

  int ndivergence = 0;
  int nnonconvergence = 0;

  typeChunks chunks;
  size_t first = 0;
  size_t last;
  _TransformationInversionParam *p = (_TransformationInversionParam*)NULL;
  int n;



  if ( _debug_ >= 2 ) fprintf( stderr, "%s: entry\n", proc );

  if ( _debug_ ) {
      fprintf( stderr, "\n" );
      BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      fprintf( stderr, "\n" );
  }



  /* check parameters
   */

  theSigma.x = theSigma.y = theSigma.z = derivationSigma;

  if ( imInverseErrors != (bal_image*)NULL ) {
      if ( imInverseErrors->ncols != invTrsf->vx.ncols || imInverseErrors->nrows != invTrsf->vx.nrows
           || imInverseErrors->nplanes != invTrsf->vx.nplanes || imInverseErrors->vdim != invTrsf->vx.vdim ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: error image and inverse transformation should have the same dimensions\n", proc );
      }
      else if ( imInverseErrors->type != UCHAR )
      {
          if ( _verbose_ )
              fprintf( stderr, "%s: error image should be of 'unsigned char' type\n", proc );
      }
      else {
          theErrors = (u8***)(imInverseErrors->array);
          for ( k=0; k<imInverseErrors->nplanes; k++ )
          for ( j=0; j<imInverseErrors->nrows; j++ )
          for ( i=0; i<imInverseErrors->ncols; i++ ) {
            theErrors[k][j][i] = 0;
          }
      }
  }



  /* allocation of derivative images
   */

  if ( BAL_AllocScalarImageFromImage( &theXdX, "XderivativeOfX.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of X component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theXdY, "YderivativeOfX.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of X component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theXdZ, "ZderivativeOfX.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of X component\n", proc );
    return( -1 );
  }

  if ( _debug_ >= 2 ) fprintf( stderr, "%s: allocation of X dZ done\n", proc );

  if ( BAL_AllocScalarImageFromImage( &theYdX, "XderivativeOfY.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Y component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theYdY, "YderivativeOfY.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Y component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theYdZ, "ZderivativeOfY.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of Y component\n", proc );
    return( -1 );
  }

  if ( _debug_ >= 2 ) fprintf( stderr, "%s: allocation of Y dZ done\n", proc );

  if ( BAL_AllocScalarImageFromImage( &theZdX, "XderivativeOfZ.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Z component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theZdY, "YderivativeOfZ.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Z component\n", proc );
    return( -1 );
  }

  if ( BAL_AllocScalarImageFromImage( &theZdZ, "ZderivativeOfZ.inr", &(theTrsf->vx), FLOAT ) != 1 ) {
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of Z component\n", proc );
    return( -1 );
  }

  if ( _debug_ >= 2 ) fprintf( stderr, "%s: allocation of Z dZ done\n", proc );




  /* calculation of derivatives
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: derivatives computation\n", proc );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vx),
                                 &theXdX, &theXdY, &theXdZ,
                                 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of X component\n", proc );
    return( -1 );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vy),
                                 &theYdX, &theYdY, &theYdZ,
                                 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Y component\n", proc );
    return( -1 );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vz),
                                 &theZdX, &theZdY, &theZdZ,
                                 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Z component\n", proc );
    return( -1 );
  }



  /* add Identity to V.Nabla^t
     and invert (Id + V.Nabla^t)
     = ( xdx + 1     xdy         xdz     )
       ( ydx         ydy + 1     ydz     )
       ( zdx         zdy         zdz + 1 )

     The inverse is computed with
     | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
     | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
     | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

     with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: matrices computation\n", proc );
  }

  /* preparing parallelism for matrix inversion ( (Id + V.Nabla^t) )
   * arrXdX, arrYdY, arrZdZ are in the geometry of theTrsf
   */

  first = 0;
  last = (theTrsf->vx).nplanes * (theTrsf->vx).nrows * (theTrsf->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
  }

  p = (_TransformationInversionParam*)vtmalloc( chunks.n_allocated_chunks * sizeof(_TransformationInversionParam),
                                                "p", proc );
  if ( p == (_TransformationInversionParam*)NULL ) {
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
      return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
      _initTransformationInversionParam( &(p[n]) );

      p[n].theTrsfDim[0] = theTrsf->vx.ncols;
      p[n].theTrsfDim[1] = theTrsf->vx.nrows;
      p[n].theTrsfDim[2] = theTrsf->vx.nplanes;

      p[n].theGeometry = theTrsf->vx.geometry;
      p[n].invGeometry = invTrsf->vx.geometry;

      p[n].the_to_voxel = &(theTrsf->vx.to_voxel);
      p[n].inv_to_real  = &(invTrsf->vx.to_real);

      p[n].arrXdX = (float***)(theXdX.array);
      p[n].arrXdY = (float***)(theXdY.array);
      p[n].arrXdZ = (float***)(theXdZ.array);
      p[n].arrYdX = (float***)(theYdX.array);
      p[n].arrYdY = (float***)(theYdY.array);
      p[n].arrYdZ = (float***)(theYdZ.array);
      p[n].arrZdX = (float***)(theZdX.array);
      p[n].arrZdY = (float***)(theZdY.array);
      p[n].arrZdZ = (float***)(theZdZ.array);

      chunks.data[n].parameters = (void*)(&(p[n]));
  }

  switch ( theTrsf->transformation_unit ) {
  default :
  case UNDEF_UNIT :
    if ( _verbose_ )
        fprintf( stderr, "%s: unknown transformation unit (3D case, matrices inversion)\n", proc );
    vtfree( p );
    freeChunks( &chunks );
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    return( -1 );
  case VOXEL_UNIT :
    if ( processChunks( &_Inverse3DVoxelMatrix, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to invert matrices (3D case, voxel)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  case REAL_UNIT :
    if ( processChunks( &_Inverse3DRealMatrix, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to invert matrices (3D case, real)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  }

  vtfree( p );
  freeChunks( &chunks );



  /* initialization of inverse vector field
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: initialization of inverse vector field\n", proc );
  }

  if ( BAL_Inverse3DVectorFieldInitialization( theTrsf, invTrsf ) != 1 ) {
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
     if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize inverse vector field\n", proc );
      return( -1 );
  }

  if ( _debug_ ) {
      fprintf( stderr, "\n" );
      BAL_PrintTransformation( stderr, invTrsf, "initial inverse transformation" );
      fprintf( stderr, "\n" );
  }

  /* preparing parallelism for vector field inversion
   */

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, " ... %s: vectors inversion\n", proc );
  }

  first = 0;
  last = (invTrsf->vx).nplanes * (invTrsf->vx).nrows * (invTrsf->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
  }

  p = (_TransformationInversionParam*)vtmalloc( chunks.n_allocated_chunks * sizeof(_TransformationInversionParam),
                                                "p", proc );
  if ( p == (_TransformationInversionParam*)NULL ) {
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
      return( -1 );
  }

  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {
      _initTransformationInversionParam( &(p[n]) );

      p[n].theTrsfDim[0] = theTrsf->vx.ncols;
      p[n].theTrsfDim[1] = theTrsf->vx.nrows;
      p[n].theTrsfDim[2] = theTrsf->vx.nplanes;

      p[n].invTrsfDim[0] = invTrsf->vx.ncols;
      p[n].invTrsfDim[1] = invTrsf->vx.nrows;
      p[n].invTrsfDim[2] = invTrsf->vx.nplanes;

      p[n].theGeometry = theTrsf->vx.geometry;
      p[n].invGeometry = invTrsf->vx.geometry;

      p[n].the_to_voxel = &(theTrsf->vx.to_voxel);
      p[n].inv_to_real  = &(invTrsf->vx.to_real);

      p[n].arrXdX = (float***)(theXdX.array);
      p[n].arrXdY = (float***)(theXdY.array);
      p[n].arrXdZ = (float***)(theXdZ.array);
      p[n].arrYdX = (float***)(theYdX.array);
      p[n].arrYdY = (float***)(theYdY.array);
      p[n].arrYdZ = (float***)(theYdZ.array);
      p[n].arrZdX = (float***)(theZdX.array);
      p[n].arrZdY = (float***)(theZdY.array);
      p[n].arrZdZ = (float***)(theZdZ.array);

      p[n].arrInvX = (float***)(invTrsf->vx.array);
      p[n].arrInvY = (float***)(invTrsf->vy.array);
      p[n].arrInvZ = (float***)(invTrsf->vz.array);

      p[n].arrTrsfX = (float***)(theTrsf->vx.array);
      p[n].arrTrsfY = (float***)(theTrsf->vy.array);
      p[n].arrTrsfZ = (float***)(theTrsf->vz.array);

      if ( theErrors != (u8***)NULL )
          p[n].theErrors = theErrors;

      p[n].ndivergence = 0;
      p[n].nnonconvergence = 0;

      chunks.data[n].parameters = (void*)(&(p[n]));
  }


  switch ( theTrsf->transformation_unit ) {
  default :
  case UNDEF_UNIT :
    if ( _verbose_ )
        fprintf( stderr, "%s: unknown transformation unit (3D case, matrices inversion)\n", proc );
    vtfree( p );
    freeChunks( &chunks );
    BAL_FreeImage( &theZdZ );
    BAL_FreeImage( &theZdY );
    BAL_FreeImage( &theZdX );
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    return( -1 );
  case VOXEL_UNIT :
    if ( processChunks( &_Inverse3DVoxelVectorField, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert vector field (3D case, voxel)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  case REAL_UNIT :
    if ( processChunks( &_Inverse3DRealVectorField, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert vector field (3D case, real)\n", proc );
      vtfree( p );
      freeChunks( &chunks );
      BAL_FreeImage( &theZdZ );
      BAL_FreeImage( &theZdY );
      BAL_FreeImage( &theZdX );
      BAL_FreeImage( &theYdZ );
      BAL_FreeImage( &theYdY );
      BAL_FreeImage( &theYdX );
      BAL_FreeImage( &theXdZ );
      BAL_FreeImage( &theXdY );
      BAL_FreeImage( &theXdX );
      return( -1 );
    }
    break;
  }

  for ( ndivergence=0, nnonconvergence=0, n=0; n<chunks.n_allocated_chunks; n++ ) {
      ndivergence += p[n].ndivergence;
      nnonconvergence += p[n].nnonconvergence;
  }

  vtfree( p );
  freeChunks( &chunks );

  if ( _debug_ ) {
      fprintf( stderr, "\n" );
      BAL_PrintTransformation( stderr, invTrsf, "result inverse transformation" );
      fprintf( stderr, "\n" );
  }


  if ( _verbose_ && (ndivergence > 0 || nnonconvergence > 0) ) {
    fprintf( stderr, "%s: divergence: %d/%lu, non-convergence: %d/%lu\n", proc,
             ndivergence,
             invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols,
             nnonconvergence,
             invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols );
  }

  BAL_FreeImage( &theZdZ );
  BAL_FreeImage( &theZdY );
  BAL_FreeImage( &theZdX );
  BAL_FreeImage( &theYdZ );
  BAL_FreeImage( &theYdY );
  BAL_FreeImage( &theYdX );
  BAL_FreeImage( &theXdZ );
  BAL_FreeImage( &theXdY );
  BAL_FreeImage( &theXdX );

  return( 1 );

}








