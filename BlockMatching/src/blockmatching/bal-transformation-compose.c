/*************************************************************************
 * bal-transformation-compose.c -
 *
 * $Id$
 *
 * Composeright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 22 jan 2017 18:53:01 CET
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
#include <bal-transformation-compose.h>





static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationCompose( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationCompose( )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationCompose( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalTransformationCompose( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalTransformationCompose( )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalTransformationCompose( )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}











/*--------------------------------------------------
 *
 * TRANSFORMATION COMPOSITION  (static operations)
 * t1 o t2 = res
 * K<-J o J<-I = K<-I
 *
 --------------------------------------------------*/



typedef struct _TransformationCompositionParam {
  bal_transformation *tkj; /* t1 */
  bal_transformation *tji; /* t2 */
  bal_transformation *tki; /* res */
} _TransformationCompositionParam;





static void *_2DVectorField2DVectorFieldCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double xj, yj;
  double xk, yk;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */
    xj = i + vjix[k][j][i];
    yj = j + vjiy[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = xj + BAL_GetXYKvalue( &(tkj->vx), xj, yj, k );
    yk = yj + BAL_GetXYKvalue( &(tkj->vy), xj, yj, k );

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_2DVectorField2DVectorFieldCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mi_to_r = tji->vx.to_real.m;
  double *mj_to_z = tkj->vx.to_voxel.m;
  double xir, yir, zir;
  double xjr, yjr, zjr;
  double xjz, yjz, zjz;
  double xkr, ykr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tji->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      zir =                                   mi_to_r[10] * k;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      zir =                                   mi_to_r[10] * k + mi_to_r[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      zir = mi_to_r[ 8] * i + mi_to_r[ 9] * j + mi_to_r[10] * k + mi_to_r[11];
      break;
    }

    /* (xjr, yjr, zjr) = M_J is the transformed point (xir, yir, zir) by tji
       in real coordinates
    */
    xjr = xir + vjix[k][j][i];
    yjr = yir + vjiy[k][j][i];
    zjr = zir;

    /* (xjz, yjz, zjz) = M_J is the same point in voxel coordinates
    */
    switch( tkj->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr;
      yjz =                     mj_to_z[ 5] * yjr;
      zjz =                                         mj_to_z[10] * zjr;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr                                         + mj_to_z[ 3];
      yjz =                     mj_to_z[ 5] * yjr                     + mj_to_z[ 7];
      zjz =                                         mj_to_z[10] * zjr + mj_to_z[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr + mj_to_z[ 1] * yjr + mj_to_z[ 2] * zjr + mj_to_z[ 3];
      yjz = mj_to_z[ 4] * xjr + mj_to_z[ 5] * yjr + mj_to_z[ 6] * zjr + mj_to_z[ 7];
      zjz = mj_to_z[ 8] * xjr + mj_to_z[ 9] * yjr + mj_to_z[10] * zjr + mj_to_z[11];
      break;
    }

    /* (xkr, ykr, zkr) is the transformed point (xjr, yjr, zjr) by tkj
       in real coordinates
    */
    xkr = xjr + BAL_GetXYKvalue( &(tkj->vx), xjz, yjz, (int)(zjz + 0.5) );
    ykr = yjr + BAL_GetXYKvalue( &(tkj->vy), xjz, yjz, (int)(zjz + 0.5) );

    /* the displacement vector is calculated by substracting
       the point M_i from M_K (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _2DVectorField2DVectorFieldComposition( bal_transformation *tkj,
                                                   bal_transformation *tji,
                                                   bal_transformation *tki )
{
  char *proc = "_2DVectorField2DVectorFieldComposition";

  size_t first = 0;
  size_t last;
  int i;
  bal_transformation tmp, *ptr;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* allocation of auxiliary transformation
   * when tkj = tki
   */

  BAL_InitTransformation( &tmp );
  ptr = (bal_transformation*)NULL;
  if ( tkj == tki ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: output is equal to 1st input\n", proc );
    if ( BAL_AllocTransformation( &tmp, VECTORFIELD_2D, &(tki->vx) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
      return( -1 );
    }
    ptr = &tmp;
  }
  else {
    ptr = tki;
  }



  /* some tests before computation
   */

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tki->vx.ncols != tji->vx.ncols || tki->vx.nrows != tji->vx.nrows
       || tki->vx.nplanes != tji->vx.nplanes ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: different dimensions for T_{J<-I} and T_{K<-I}\n", proc );
    }
    return( -1 );
  }

  if ( tkj->type != VECTORFIELD_2D || tji->type != VECTORFIELD_2D || tki->type != VECTORFIELD_2D ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */

  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = ptr;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);



  /* computation
   */

  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    freeChunks( &chunks );
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    return( -1 );

  case VOXEL_UNIT :

    if ( processChunks( &_2DVectorField2DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      if ( tkj == tki ) BAL_FreeTransformation( &tmp );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_2DVectorField2DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      if ( tkj == tki ) BAL_FreeTransformation( &tmp );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );



  /* finalize result
   */

  ptr->transformation_unit = tji->transformation_unit;

  if ( tkj == tki ) {
    if ( BAL_CopyTransformation( &tmp, tki ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy auxiliary transformation\n", proc );
      BAL_FreeTransformation( &tmp );
      return( -1 );
    }
    BAL_FreeTransformation( &tmp );
  }

  /* copy geometry characteristics of tji into tki
   * (not sure it is necessary)
   */
  if ( BAL_CopyImageGeometry( &(tji->vx), &(tki->vx) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vx'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vy), &(tki->vy) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vy'\n", proc );
    return( -1 );
  }

  return( 1 );
}





static void *_3DVectorField3DVectorFieldCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;
  float ***vjiz = (float ***)(tji->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double xj, yj, zj;
  double xk, yk, zk;


  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */
    xj = i + vjix[k][j][i];
    yj = j + vjiy[k][j][i];
    zj = k + vjiz[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = xj + BAL_GetXYZvalue( &(tkj->vx), xj, yj, zj );
    yk = yj + BAL_GetXYZvalue( &(tkj->vy), xj, yj, zj );
    zk = zj + BAL_GetXYZvalue( &(tkj->vz), xj, yj, zj );

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;
    resz[k][j][i] = zk - (double)k;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_3DVectorField3DVectorFieldCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;
  float ***vjiz = (float ***)(tji->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mi_to_r = tji->vx.to_real.m;
  double *mj_to_z = tkj->vx.to_voxel.m;
  double xir, yir, zir;
  double xjr, yjr, zjr;
  double xjz, yjz, zjz;
  double xkr, ykr, zkr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tji->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      zir =                                   mi_to_r[10] * k;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      zir =                                   mi_to_r[10] * k + mi_to_r[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      zir = mi_to_r[ 8] * i + mi_to_r[ 9] * j + mi_to_r[10] * k + mi_to_r[11];
      break;
    }


    /* (xjr, yjr, zjr) = M_J is the transformed point (xir, yir, zir) by tji
       in real coordinates
    */
    xjr = xir + vjix[k][j][i];
    yjr = yir + vjiy[k][j][i];
    zjr = zir + vjiz[k][j][i];

    /* (xjz, yjz, zjz) = M_J is the same point in voxel coordinates
    */
    switch( tkj->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr;
      yjz =                     mj_to_z[ 5] * yjr;
      zjz =                                         mj_to_z[10] * zjr;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr                                         + mj_to_z[ 3];
      yjz =                     mj_to_z[ 5] * yjr                     + mj_to_z[ 7];
      zjz =                                         mj_to_z[10] * zjr + mj_to_z[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr + mj_to_z[ 1] * yjr + mj_to_z[ 2] * zjr + mj_to_z[ 3];
      yjz = mj_to_z[ 4] * xjr + mj_to_z[ 5] * yjr + mj_to_z[ 6] * zjr + mj_to_z[ 7];
      zjz = mj_to_z[ 8] * xjr + mj_to_z[ 9] * yjr + mj_to_z[10] * zjr + mj_to_z[11];
      break;
    }

    /* (xkr, ykr, zkr) is the transformed point (xjr, yjr, zjr) by tkj
       in real coordinates
    */
    xkr = xjr + BAL_GetXYZvalue( &(tkj->vx), xjz, yjz, zjz );
    ykr = yjr + BAL_GetXYZvalue( &(tkj->vy), xjz, yjz, zjz );
    zkr = zjr + BAL_GetXYZvalue( &(tkj->vz), xjz, yjz, zjz );

    /* the displacement vector is calculated by substracting
       the point M_i from M_K (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;
    resz[k][j][i] = zkr - zir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _3DVectorField3DVectorFieldComposition( bal_transformation *tkj,
                                                   bal_transformation *tji,
                                                   bal_transformation *tki )
{
  char *proc = "_3DVectorField3DVectorFieldComposition";

  size_t first = 0;
  size_t last;
  int i;
  bal_transformation tmp, *ptr;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* allocation of auxiliary transformation
   * when tkj = tki
   */

  BAL_InitTransformation( &tmp );
  ptr = (bal_transformation*)NULL;
  if ( tkj == tki ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: output is equal to 1st input\n", proc );
    if ( BAL_AllocTransformation( &tmp, VECTORFIELD_3D, &(tki->vx) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
      return( -1 );
    }
    ptr = &tmp;
  }
  else {
    ptr = tki;
  }



  /* some tests before computation
   */

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tki->vx.ncols != tji->vx.ncols || tki->vx.nrows != tji->vx.nrows
       || tki->vx.nplanes != tji->vx.nplanes ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: different dimensions for T_{J<-I} and T_{K<-I}\n", proc );
    }
    return( -1 );
  }

  if ( tkj->type != VECTORFIELD_3D || tji->type != VECTORFIELD_3D || tki->type != VECTORFIELD_3D ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types\n", proc );
    return( -1 );
  }

  if ( tkj->vx.nplanes <=1 || tji->vx.nplanes <=1 || tki->vx.nplanes <=1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    return( -1 );
  }



  /* preparing parallelism
   */

  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = ptr;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);



  /* computation
   */

  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    freeChunks( &chunks );
    if ( tkj == tki ) BAL_FreeTransformation( &tmp );
    return( -1 );

  case VOXEL_UNIT :

    if ( processChunks( &_3DVectorField3DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      if ( tkj == tki ) BAL_FreeTransformation( &tmp );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_3DVectorField3DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      if ( tkj == tki ) BAL_FreeTransformation( &tmp );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );



  /* finalize result
   */

  ptr->transformation_unit = tji->transformation_unit;

  if ( tkj == tki ) {
    if ( BAL_CopyTransformation( &tmp, tki ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy auxiliary transformation\n", proc );
      BAL_FreeTransformation( &tmp );
      return( -1 );
    }
    BAL_FreeTransformation( &tmp );
  }

  /* copy geometry characteristics of tji into tki
   * (not sure it is necessary)
   */
  if ( BAL_CopyImageGeometry( &(tji->vx), &(tki->vx) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vx'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vy), &(tki->vy) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vy'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vz), &(tki->vz) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vz'\n", proc );
    return( -1 );
  }
  return( 1 );
}





static void *_2DVectorFieldMatrixCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mji = tji->mat.m;
  double xj, yj;
  double xk, yk;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */

    xj = mji[ 0] * i + mji[ 1] * j + mji[ 3];
    yj = mji[ 4] * i + mji[ 5] * j + mji[ 7];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = xj + BAL_GetXYKvalue( &(tkj->vx), xj, yj, k );
    yk = yj + BAL_GetXYKvalue( &(tkj->vy), xj, yj, k );

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_2DVectorFieldMatrixCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mji = tji->mat.m;
  double *mi_to_r = tki->vx.to_real.m;
  double *mj_to_z = tkj->vx.to_voxel.m;
  double xir, yir, zir;
  double xjr, yjr;
  double xjz, yjz, zjz;
  double xkr, ykr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tki->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      zir =                                   mi_to_r[10] * k;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      zir =                                   mi_to_r[10] * k + mi_to_r[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      zir = mi_to_r[ 8] * i + mi_to_r[ 9] * j + mi_to_r[10] * k + mi_to_r[11];
      break;
    }

    /* (xjr, yjr, zjr) is the transformed point (i,j,k) by tji
       in real coordinates
    */

    xjr = mji[ 0] * xir + mji[ 1] * yir + mji[ 2] * zir + mji[ 3];
    yjr = mji[ 4] * xir + mji[ 5] * yir + mji[ 2] * zir + mji[ 7];
    /*
    zjr = mji[ 8] * xir + mji[ 9] * yir + mji[10] * zir + mji[11];
    */

    /* (xjz, yjz, zjz) = M_J is the same point in voxel coordinates
    */
    switch( tkj->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr;
      yjz =                     mj_to_z[ 5] * yjr;
      zjz =                                         mj_to_z[10];
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr                                       + mj_to_z[ 3];
      yjz =                     mj_to_z[ 5] * yjr                   + mj_to_z[ 7];
      zjz =                                         mj_to_z[10] * k + mj_to_z[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr + mj_to_z[ 1] * yjr + mj_to_z[ 2] * k + mj_to_z[ 3];
      yjz = mj_to_z[ 4] * xjr + mj_to_z[ 5] * yjr + mj_to_z[ 6] * k + mj_to_z[ 7];
      zjz = mj_to_z[ 8] * xjr + mj_to_z[ 9] * yjr + mj_to_z[10] * k + mj_to_z[11];
      break;
    }

    /* (xkr, ykr, zkr) is the transformed point (xjr, yjr, zjr) by tkj
       in real coordinates
    */
    xkr = xjr + BAL_GetXYKvalue( &(tkj->vx), xjz, yjz, (int)(zjz + 0.5) );
    ykr = yjr + BAL_GetXYKvalue( &(tkj->vy), xjz, yjz, (int)(zjz + 0.5) );

    /* the displacement vector is calculated by substracting
       the point M_i from M_K (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _2DVectorFieldMatrixComposition( bal_transformation *tkj,
                                            bal_transformation *tji,
                                            bal_transformation *tki )
{
  char *proc = "_2DVectorFieldMatrixComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* some tests before computation
   */

  if ( tkj == tki ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tkj->type != VECTORFIELD_2D || tki->type != VECTORFIELD_2D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 1st input and/or result\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( tji ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 2nd input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = tki;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);




  /* computation
   */

  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :

    if ( processChunks( &_2DVectorFieldMatrixCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_2DVectorFieldMatrixCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );



  /* finalize result
   */

  tki->transformation_unit = tji->transformation_unit;

  return( 1 );
}





static void *_3DVectorFieldMatrixCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mji = tji->mat.m;
  double xj, yj, zj;
  double xk, yk, zk;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */
    xj = mji[ 0] * i + mji[ 1] * j + mji[ 2] * k + mji[ 3];
    yj = mji[ 4] * i + mji[ 5] * j + mji[ 6] * k + mji[ 7];
    zj = mji[ 8] * i + mji[ 9] * j + mji[10] * k + mji[11];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = xj + BAL_GetXYKvalue( &(tkj->vx), xj, yj, zj );
    yk = yj + BAL_GetXYKvalue( &(tkj->vy), xj, yj, zj );
    zk = zj + BAL_GetXYKvalue( &(tkj->vz), xj, yj, zj );

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;
    resz[k][j][i] = zk - (double)k;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_3DVectorFieldMatrixCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mji = tji->mat.m;
  double *mi_to_r = tki->vx.to_real.m;
  double *mj_to_z = tkj->vx.to_voxel.m;
  double xir, yir, zir;
  double xjr, yjr, zjr;
  double xjz, yjz, zjz;
  double xkr, ykr, zkr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tki->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      zir =                                   mi_to_r[10] * k;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      zir =                                   mi_to_r[10] * k + mi_to_r[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      zir = mi_to_r[ 8] * i + mi_to_r[ 9] * j + mi_to_r[10] * k + mi_to_r[11];
      break;
    }

    /* (xjr, yjr, zjr) is the transformed point (i,j,k) by tji
       in real coordinates
    */

    xjr = mji[ 0] * xir + mji[ 1] * yir + mji[ 2] * zir + mji[ 3];
    yjr = mji[ 4] * xir + mji[ 5] * yir + mji[ 2] * zir + mji[ 7];
    zjr = mji[ 8] * xir + mji[ 9] * yir + mji[10] * zir + mji[11];

    /* (xjz, yjz, zjz) = M_J is the same point in voxel coordinates
    */
    switch( tkj->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr;
      yjz =                     mj_to_z[ 5] * yjr;
      zjz =                                         mj_to_z[10];
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr                                       + mj_to_z[ 3];
      yjz =                     mj_to_z[ 5] * yjr                   + mj_to_z[ 7];
      zjz =                                         mj_to_z[10] * k + mj_to_z[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xjz = mj_to_z[ 0] * xjr + mj_to_z[ 1] * yjr + mj_to_z[ 2] * k + mj_to_z[ 3];
      yjz = mj_to_z[ 4] * xjr + mj_to_z[ 5] * yjr + mj_to_z[ 6] * k + mj_to_z[ 7];
      zjz = mj_to_z[ 8] * xjr + mj_to_z[ 9] * yjr + mj_to_z[10] * k + mj_to_z[11];
      break;
    }

    /* (xkr, ykr, zkr) is the transformed point (xjr, yjr, zjr) by tkj
       in real coordinates
    */
    xkr = xjr + BAL_GetXYZvalue( &(tkj->vx), xjz, yjz, zjz );
    ykr = yjr + BAL_GetXYZvalue( &(tkj->vy), xjz, yjz, zjz );
    zkr = zjr + BAL_GetXYZvalue( &(tkj->vy), xjz, yjz, zjz );

    /* the displacement vector is calculated by substracting
       the point M_i from M_K (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;
    resz[k][j][i] = zkr - zir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _3DVectorFieldMatrixComposition( bal_transformation *tkj,
                                            bal_transformation *tji,
                                            bal_transformation *tki )
{
  char *proc = "_3DVectorFieldMatrixComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* some tests before computation
   */

  if ( tkj == tki ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tkj->type != VECTORFIELD_3D || tki->type != VECTORFIELD_3D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 1st input and/or result\n", proc );
    return( -1 );
  }

  if ( tkj->vx.nplanes <=1 || tki->vx.nplanes <=1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( tji ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 2nd input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = tki;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);



  /* computation
   */

  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :

    if ( processChunks( &_3DVectorFieldMatrixCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_3DVectorFieldMatrixCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );

  /* finalize result
   */

  tki->transformation_unit = tji->transformation_unit;

  return( 1 );
}





static void *_Matrix2DVectorFieldCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mkj = tkj->mat.m;
  double xj, yj;
  double xk, yk;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */
    xj = i + vjix[k][j][i];
    yj = j + vjiy[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = mkj[ 0] * xj + mkj[ 1] * yj + mkj[ 3];
    yk = mkj[ 4] * xj + mkj[ 5] * yj + mkj[ 7];

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_Matrix2DVectorFieldCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mkj = tkj->mat.m;
  double *mi_to_r = tji->vx.to_real.m;
  double xir, yir;
  double xjr, yjr;
  double xkr, ykr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tji->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      break;
    }

    /* (xjr, yjr, zjr) = M_J is the transformed point (xir, yir, zir) by tji
       in real coordinates
    */
    xjr = xir + vjix[k][j][i];
    yjr = yir + vjiy[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in real coordinates
    */
    xkr = mkj[ 0] * xjr + mkj[ 1] * yjr + mkj[ 3];
    ykr = mkj[ 4] * xjr + mkj[ 5] * yjr + mkj[ 7];

    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _Matrix2DVectorFieldComposition( bal_transformation *tkj,
                                            bal_transformation *tji,
                                            bal_transformation *tki )
{
  char *proc = "_Matrix2DVectorFieldComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* some tests before computation
   */

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tji->type != VECTORFIELD_2D || tki->type != VECTORFIELD_2D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 2nd input and/or result\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( tkj ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 1st input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = tki;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);



  /* computation
   */

  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :

    if ( processChunks( &_Matrix2DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_Matrix2DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );

  /* finalize result
   */

  tki->transformation_unit = tji->transformation_unit;

  /* copy geometry characteristics of tji into tki
   * (not sure it is necessary)
   */
  if ( BAL_CopyImageGeometry( &(tji->vx), &(tki->vx) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vx'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vy), &(tki->vy) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vy'\n", proc );
    return( -1 );
  }
  return( 1 );
}





static void *_Matrix3DVectorFieldCompositionVoxelUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;
  float ***vjiz = (float ***)(tji->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mkj = tkj->mat.m;
  double xj, yj, zj;
  double xk, yk, zk;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xj, yj, zj) is the transformed point (i,j,k) by tji
       in voxel coordinates
    */
    xj = i + vjix[k][j][i];
    yj = j + vjiy[k][j][i];
    zj = k + vjiz[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in voxel coordinates
    */
    xk = mkj[ 0] * xj + mkj[ 1] * yj + mkj[ 2] * zj + mkj[ 3];
    yk = mkj[ 4] * xj + mkj[ 5] * yj + mkj[ 6] * zj + mkj[ 7];
    zk = mkj[ 8] * xj + mkj[ 8] * yj + mkj[10] * zj + mkj[11];

    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = xk - (double)i;
    resy[k][j][i] = yk - (double)j;
    resz[k][j][i] = zk - (double)k;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_Matrix3DVectorFieldCompositionRealUnit( void *par )
{
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *tkj = p->tkj;
  bal_transformation *tji = p->tji;
  bal_transformation *tki = p->tki;
  float ***resx = (float ***)(tki->vx).array;
  float ***resy = (float ***)(tki->vy).array;
  float ***resz = (float ***)(tki->vz).array;
  float ***vjix = (float ***)(tji->vx).array;
  float ***vjiy = (float ***)(tji->vy).array;
  float ***vjiz = (float ***)(tji->vz).array;

  size_t ncols = (tki->vx).ncols; /* dimx */
  size_t nrows = (tki->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  size_t i, j, k;
  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;

  double *mkj = tkj->mat.m;
  double *mi_to_r = tji->vx.to_real.m;
  double xir, yir, zir;
  double xjr, yjr, zjr;
  double xkr, ykr, zkr;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (xir, yir, zir) is M_I in real coordinates
     */
    switch( tji->vx.geometry ) {
    default :
    case _BAL_UNKNOWN_GEOMETRY_ :
    case _BAL_HOMOTHETY_GEOMETRY_ :
      xir = mi_to_r[ 0] * i;
      yir =                   mi_to_r[ 5] * j;
      zir =                                   mi_to_r[10] * k;
      break;
    case _BAL_TRANSLATION_GEOMETRY_ :
      xir = mi_to_r[ 0] * i                                   + mi_to_r[ 3];
      yir =                   mi_to_r[ 5] * j                 + mi_to_r[ 7];
      zir =                                   mi_to_r[10] * k + mi_to_r[11];
      break;
    case _BAL_QFORM_GEOMETRY_ :
      xir = mi_to_r[ 0] * i + mi_to_r[ 1] * j + mi_to_r[ 2] * k + mi_to_r[ 3];
      yir = mi_to_r[ 4] * i + mi_to_r[ 5] * j + mi_to_r[ 6] * k + mi_to_r[ 7];
      zir = mi_to_r[ 8] * i + mi_to_r[ 9] * j + mi_to_r[10] * k + mi_to_r[11];
      break;
    }

    /* (xjr, yjr, zjr) = M_J is the transformed point (xir, yir, zir) by tji
       in real coordinates
    */
    xjr = xir + vjix[k][j][i];
    yjr = yir + vjiy[k][j][i];
    zjr = zir + vjiz[k][j][i];

    /* (xk, yk, zk) is the transformed point (xj, yj, zj) by tkj
       in real coordinates
    */
    xkr = mkj[ 0] * xjr + mkj[ 1] * yjr + mkj[ 2] * zjr + mkj[ 3];
    ykr = mkj[ 4] * xjr + mkj[ 5] * yjr + mkj[ 6] * zjr + mkj[ 7];
    zkr = mkj[ 8] * xjr + mkj[ 9] * yjr + mkj[10] * zjr + mkj[11];

    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = xkr - xir;
    resy[k][j][i] = ykr - yir;
    resz[k][j][i] = zkr - zir;

  }
  chunk->ret = 1;
  return( (void*)NULL );
}





/* tkj o tji = tki
 */
static int _Matrix3DVectorFieldComposition( bal_transformation *tkj,
                                            bal_transformation *tji,
                                            bal_transformation *tki )
{
  char *proc = "_Matrix3DVectorFieldComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;



  /* some tests before computation
   */

  if ( tkj->transformation_unit != tji->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( tji->type != VECTORFIELD_3D || tki->type != VECTORFIELD_3D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 2nd input and/or result\n", proc );
    return( -1 );
  }

  if ( tji->vx.nplanes <= 1 || tki->vx.nplanes <= 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( tkj ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 1st input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (tki->vx).nplanes * (tki->vx).nrows * (tki->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.tkj = tkj;
  p.tji = tji;
  p.tki = tki;
  for ( i=0; i<chunks.n_allocated_chunks; i++ )
    chunks.data[i].parameters = (void*)(&p);



  switch ( tkj->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :

    if ( processChunks( &_Matrix3DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_Matrix3DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );

  /* finalize result
   */

  tki->transformation_unit = tji->transformation_unit;

  /* copy geometry characteristics of tji into tki
   * (not sure it is necessary)
   */
  if ( BAL_CopyImageGeometry( &(tji->vx), &(tki->vx) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vx'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vy), &(tki->vy) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vy'\n", proc );
    return( -1 );
  }
  if ( BAL_CopyImageGeometry( &(tji->vz), &(tki->vz) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to copy image geometry to 'tki->vz'\n", proc );
    return( -1 );
  }

  return( 1 );
}













/*--------------------------------------------------
 *
 * TRANSFORMATION COMPOSITION
 *
 --------------------------------------------------*/



int BAL_TransformationComposition( bal_transformation *t_res, /* t_res = t1 o t2 */
                                   bal_transformation *t1,
                                   bal_transformation *t2 )
{
  char *proc = "BAL_TransformationComposition";
  _MATRIX tmp;

  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  switch ( t1->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such first transformation type not handled yet\n", proc );
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

    switch ( t2->type ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (matrix o unknown)\n", proc );
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

      if ( _alloc_mat( &tmp, 4, 4) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary matrix\n", proc );
        return( -1 );
      }
      _mult_mat( &(t1->mat), &(t2->mat), &tmp );
      _copy_mat( &tmp, &(t_res->mat) );
      _free_mat( &tmp );
      break;

    case VECTORFIELD_2D :

      if ( _Matrix2DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing 'matrix' o '2D vector field'\n", proc );
        return( -1 );
      }
      break;

    case VECTORFIELD_3D :

      if ( _Matrix3DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing 'matrix' o '3D vector field'\n", proc );
        return( -1 );
      }
      break;

    }

    break;
    /* end of matrix case
     */

  case VECTORFIELD_2D :

    switch ( t2->type ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (2D vector field o unknown)\n", proc );
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
      if ( _2DVectorFieldMatrixComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing '2D vector field' o 'matrix'\n", proc );
        return( -1 );
      }
      break;

    case VECTORFIELD_2D:
      if ( _2DVectorField2DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing '2D vector field' o '2D vector field'\n", proc );
        return( -1 );
      };
      break;
    }

    break;
    /* end of 2D vector field case
     */

  case VECTORFIELD_3D :

    switch ( t2->type ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (3D vector field o unknown)\n", proc );
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
      if ( _3DVectorFieldMatrixComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing '3D vector field' o 'matrix'\n", proc );
        return( -1 );
      }
      break;

    case VECTORFIELD_3D:
      if ( _3DVectorField3DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when composing '3D vector field' o '3D vector field'\n", proc );
        return( -1 );
      };
      break;
    }

    break;
    /* end of 3D vector field case
     */

  }

  t_res->transformation_unit =  t1->transformation_unit;

  return( 1 );
}





int BAL_TransformationListComposition( bal_transformation *t_res,
                                      bal_transformation **t_array,
                                      int n )
{
    char *proc = "BAL_TransformatioListComposition";
    bal_transformation tmpTrsf;
    int i;

    if ( BAL_CopyTransformation( t_array[n-1], t_res ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy first transformation\n", proc );
        return( -1 );
    }

    /* loop over all other transformations
     */
    for ( i = n-2; i >=0; i-- ) {

        if ( BAL_IsTransformationLinear( t_res ) == 1 ) {
            if ( BAL_IsTransformationVectorField( t_array[i] ) == 1 ) {

                BAL_InitTransformation( &tmpTrsf );

                if ( BAL_AllocTransformation( &tmpTrsf, t_res->type, (bal_image*)NULL ) != 1 ) {
                    BAL_FreeTransformation( &tmpTrsf );
                    if ( _verbose_ )
                        fprintf( stderr, "%s: unable to allocate transformation (temporary backup)\n", proc );
                    return( -1 );
                }

                if ( BAL_CopyTransformation( t_res, &tmpTrsf ) != 1 ) {
                    BAL_FreeTransformation( &tmpTrsf );
                    if ( _verbose_ )
                        fprintf( stderr, "%s: unable to copy transformation (temporary backup)\n", proc );
                    return( -1 );
                }

                BAL_FreeTransformation( t_res );

                /* allocation of a new result transformation
                 */
                if ( BAL_AllocTransformation( t_res, t_array[i]->type, &(t_array[i]->vx) ) != 1 ) {
                  BAL_FreeTransformation( &tmpTrsf );
                  if ( _verbose_ )
                      fprintf( stderr, "%s: unable to allocate transformation (temporary backup)\n", proc );
                  return( -1 );
                }

                if ( BAL_CopyTransformation( &tmpTrsf, t_res ) != 1 ) {
                    BAL_FreeTransformation( &tmpTrsf );
                    if (_verbose_ )
                        fprintf( stderr, "%s: unable to copy transformation\n", proc );
                    return( -1 );
                }

                BAL_FreeTransformation( &tmpTrsf );
            }
        }

        /* compose transformations
         * int BAL_TransformationComposition( bal_transformation *t_res,
         *                                    bal_transformation *t1,
         *                                    bal_transformation *t2 )
         * with t_res = t1 o t2
         */

        if ( BAL_TransformationComposition( t_res, t_array[i], t_res ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to compose transformation #%d with intermediary result\n", proc, i );
          return( -1 ) ;
        }

    }

    return( 1 );
}





static enumTypeTransfo _TypeTransformationComposition( enumTypeTransfo type1,
                                                       enumTypeTransfo type2 )
{
  char *proc = "_TypeTransformationComposition";

  switch ( type1 ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such first transformation type not handled yet\n", proc );
    return( UNDEF_TRANSFORMATION );

  case TRANSLATION_2D :
  case TRANSLATION_SCALING_2D :
  case RIGID_2D :
  case SIMILITUDE_2D :
  case AFFINE_2D :

    /* 2D matrix case */
    switch ( type2 ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (2D matrix o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );

    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
      return( AFFINE_2D );

    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
      return( AFFINE_3D );

    case VECTORFIELD_2D :
      return( VECTORFIELD_2D );

    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 2D matrix case */

  case TRANSLATION_3D :
  case TRANSLATION_SCALING_3D :
  case RIGID_3D :
  case SIMILITUDE_3D :
  case AFFINE_3D :

    /* 3D matrix case */
    switch ( type2 ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (3D matrix o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );

    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
      return( AFFINE_3D );

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 3D matrix case */

  case VECTORFIELD_2D :

    /* 2D vector field case */
    switch ( type2 ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (2D vector field o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );

    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
      return( VECTORFIELD_2D );

    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
      return( VECTORFIELD_3D );

    case VECTORFIELD_2D :
      return( VECTORFIELD_2D );

    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 2D vector field case */

  case VECTORFIELD_3D :

    /* 3D vector field case */
    switch ( type2 ) {

    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such second transformation type not handled yet (3D vector field o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );

    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 3D vector field case */

  }

  if ( _verbose_ )
    fprintf( stderr, "%s: should not reach this ...\n", proc );
  return( UNDEF_TRANSFORMATION );
}





enumTypeTransfo BAL_TypeTransformationComposition( bal_transformation *t1,
                                                   bal_transformation *t2 )
{
    return( _TypeTransformationComposition( t1->type, t2->type ) );
}





/* in this procedure, res is preferably allocated accordingly to
 * first t2, then t1, and last ref
 */
int BAL_AllocTransformationComposition( bal_transformation *res,
                                        bal_transformation *t1,
                                        bal_transformation *t2,
                                        bal_image *ref )
{
  char *proc = "BAL_AllocTransformationComposition";
  enumTypeTransfo resType;

  resType = BAL_TypeTransformationComposition( t1, t2 );

   switch ( resType ) {

    default :
      if ( _verbose_ ) {
        fprintf( stderr, "%s: such type '", proc );
        BAL_PrintTransformationType( stderr, resType );
        fprintf( stderr, "' not handled yet\n" );
      }
      return( -1 );

    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D :
      if ( BAL_AllocTransformation( res, resType, (bal_image *)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate transformation (linear case)\n", proc );
        return( -1 );
      }
      return( 1 );

    case VECTORFIELD_2D :
      if ( t2->type == VECTORFIELD_2D ) {
        if ( BAL_AllocTransformation( res, resType, &(t2->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (from t2))\n", proc );
          return( -1 );
        }
        return( 1 );
      }
      else if ( t1->type == VECTORFIELD_2D ) {
        if ( BAL_AllocTransformation( res, resType, &(t1->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (from t1))\n", proc );
          return( -1 );
        }
        return( 1 );
      }
      else if ( ref != (bal_image *)NULL ) {
        if ( BAL_AllocTransformation( res, resType, ref ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (from template))\n", proc );
          return( -1 );
        }
        return( 1 );
      }

      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate transformation (2D vector field case, weird case)\n", proc );
      return( -1 );

    case VECTORFIELD_3D :
      if ( t2->type == VECTORFIELD_3D ) {
        if ( BAL_AllocTransformation( res, resType, &(t2->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (from t2))\n", proc );
          return( -1 );
        }
        return( 1 );
      }
      else if ( t1->type == VECTORFIELD_3D ) {
        if ( BAL_AllocTransformation( res, resType, &(t1->vx) ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (from t1))\n", proc );
          return( -1 );
        }
        return( 1 );
      }
      else if ( ref != (bal_image *)NULL ) {
        if ( BAL_AllocTransformation( res, resType, ref ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (from template))\n", proc );
          return( -1 );
        }
        return( 1 );
      }

      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate transformation (3D vector field case, weird case)\n", proc );
      return( -1 );

   }

  if ( _verbose_ )
    fprintf( stderr, "%s: should not reach this ...\n", proc );
   return( -1 );
}





enumTypeTransfo BAL_TypeTransformationListComposition( bal_transformation **array, int n )
{
    enumTypeTransfo resType = UNDEF_TRANSFORMATION;
    enumTypeTransfo prevType;
    int i;

    if ( n == 0 )
        return( UNDEF_TRANSFORMATION );
    if ( n == 1 )
        return( array[0]->type );

    prevType = array[n-1]->type;
    for ( i=n-2; i>=0; i-- ) {
        resType = _TypeTransformationComposition( array[i]->type, prevType );
        prevType = resType;
    }
    return( resType );
}





int BAL_AllocTransformationListComposition( bal_transformation *res,
                                            bal_transformation **array,
                                            int n,
                                            bal_image *ref )
{
    char *proc = "BAL_AllocTransformationListComposition";
    enumTypeTransfo resType;
    int i;
    bal_transformation *reftrsf = (bal_transformation*)NULL;

    resType = BAL_TypeTransformationListComposition( array, n );

    switch ( resType ) {

     default :
       if ( _verbose_ ) {
         fprintf( stderr, "%s: such type '", proc );
         BAL_PrintTransformationType( stderr, resType );
         fprintf( stderr, "' not handled yet\n" );
       }
       return( -1 );

     case TRANSLATION_2D :
     case TRANSLATION_SCALING_2D :
     case RIGID_2D :
     case SIMILITUDE_2D :
     case AFFINE_2D :
     case TRANSLATION_3D :
     case TRANSLATION_SCALING_3D :
     case RIGID_3D :
     case SIMILITUDE_3D :
     case AFFINE_3D :
       if ( BAL_AllocTransformation( res, resType, (bal_image *)NULL ) != 1 ) {
         if ( _verbose_ )
           fprintf( stderr, "%s: unable to allocate transformation (linear case)\n", proc );
         return( -1 );
       }
       return( 1 );

     case VECTORFIELD_2D :
     case VECTORFIELD_3D :
       /* 1. check for a template
        * 2. check for a vectorfield transformation in the list
        */
       if ( ref != (bal_image *)NULL ) {
         if ( BAL_AllocTransformation( res, resType, ref ) != 1 ) {
           if ( _verbose_ )
             fprintf( stderr, "%s: unable to allocate transformation (vector field case (from template))\n", proc );
           return( -1 );
         }
         return( 1 );
       }

       for ( i=n-1; i>=0 && reftrsf == (bal_transformation*)NULL; i-- ) {
         if ( array[i]->type == VECTORFIELD_2D
              || array[i]->type == VECTORFIELD_3D )
           reftrsf = array[i];
       }
       if ( reftrsf == (bal_transformation*)NULL ) {
         if ( _verbose_ )
           fprintf( stderr, "%s: unable to find a template to allocate transformation (vector field case, weird case)\n", proc );
         return( -1 );
       }
       if ( BAL_AllocTransformation( res, resType, &(reftrsf->vx) ) != 1 ) {
         if ( _verbose_ )
           fprintf( stderr, "%s: unable to allocate transformation (vector field case (from list))\n", proc );
         return( -1 );
       }
       return( 1 );
    }

   if ( _verbose_ )
     fprintf( stderr, "%s: should not reach this ...\n", proc );
    return( -1 );
}


