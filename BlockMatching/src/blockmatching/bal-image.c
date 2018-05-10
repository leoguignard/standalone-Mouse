/*************************************************************************
 * bal-image.c -
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




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include <ImageIO.h>

#include <convert.h>
#include <linearFiltering.h>
#include <vtmalloc.h>

#include <bal-image.h>


static int _verbose_ = 1;
static int _debug_ = 0;




int BAL_GetVerboseInBalImage(  )
{
  return( _verbose_ );
}

void BAL_SetVerboseInBalImage( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalImage(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalImage(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalImage( int d )
{
  _debug_ = d;
}

void BAL_IncrementDebugInBalImage(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalImage(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}






/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image geometry
 *
 --------------------------------------------------*/





static void BAL_InitImageGeometry( bal_image *image )
{
  image->vx = 1.0;
  image->vy = 1.0;
  image->vz = 1.0;

  image->geometry = _BAL_UNKNOWN_GEOMETRY_;

  /* matrices
   * - from voxel grid (integer) to real space
   * - from real space to voxel grid
   */
  _init_mat( &(image->to_real) );
  _init_mat( &(image->to_voxel) );
  image->qform_code = 0;
  image->sform_code = 0;
  _init_mat( &(image->sform_to_real) );
}





static void BAL_FreeImageGeometry( bal_image *image )
{
  _free_mat( &(image->to_real) );
  _free_mat( &(image->to_voxel) );
  _free_mat( &(image->sform_to_real) );

  BAL_InitImageGeometry( image );
}





int BAL_AllocImageGeometry( bal_image *image )
{
  char *proc = "BAL_AllocImageGeometry";

  /* desallocation of matrices if dimensions are not correct
   */
  if ( (image->to_real.l > 0 && image->to_real.l != 4)
       || (image->to_real.c > 0 && image->to_real.c != 4) )
    _free_mat( &(image->to_real) );

  if ( (image->to_voxel.l > 0 && image->to_voxel.l != 4)
       || (image->to_voxel.c > 0 && image->to_voxel.c != 4) )
    _free_mat( &(image->to_voxel) );

  if ( (image->sform_to_real.l > 0 && image->sform_to_real.l != 4)
       || (image->sform_to_real.c > 0 && image->sform_to_real.c != 4) )
    _free_mat( &(image->sform_to_real) );

  /* allocations
   */
  if ( image->to_real.m == (double*)NULL ) {
    if ( _alloc_mat( &(image->to_real), 4, 4 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: can not allocate 'to_real' matrix\n", proc );
        return( -1 );
    }
  }

  if ( image->to_voxel.m == (double*)NULL ) {
    if ( _alloc_mat( &(image->to_voxel), 4, 4 ) != 1 ) {
      _free_mat( &(image->to_real) );
        if ( _verbose_ )
          fprintf( stderr, "%s: can not allocate 'to_voxel' matrix\n", proc );
        return( -1 );
    }
  }

  if ( image->sform_to_real.m == (double*)NULL ) {
    if ( _alloc_mat( &(image->sform_to_real), 4, 4 ) != 1 ) {
      _free_mat( &(image->to_voxel) );
      _free_mat( &(image->to_real) );
      if ( _verbose_ )
        fprintf( stderr, "%s: can not allocate 'sform_to_real' matrix\n", proc );
      return( -1 );
    }
  }

  /* set matrices to identity
   */
  image->vx = 1.0;
  image->vy = 1.0;
  image->vz = 1.0;
  _identity_mat( &(image->to_real) );
  _identity_mat( &(image->to_voxel) );
  _identity_mat( &(image->sform_to_real) );

  image->geometry = _BAL_HOMOTHETY_GEOMETRY_;

  return( 1 );
}





int BAL_SetImageVoxelSizes( bal_image *image, typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz )
{
  char *proc = "BAL_SetImageVoxelSizes";

  switch ( image->geometry ) {
  case _BAL_UNKNOWN_GEOMETRY_ :
    /* geometry is set to identity
     */
    if ( BAL_AllocImageGeometry( image ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate image geometry\n", proc );
      return( -1 );
    }
  case _BAL_HOMOTHETY_GEOMETRY_ :
    break;
  case _BAL_TRANSLATION_GEOMETRY_ :
  case _BAL_QFORM_GEOMETRY_ :
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image orientation information may be inconsistent\n", proc );
  }

  /* recompute geometry information
   */
  image->vx = vx;
  image->vy = vy;
  image->vz = vz;

  _identity_mat( &(image->to_real) );
  image->to_real.m[ 0] = vx;
  image->to_real.m[ 5] = vy;
  image->to_real.m[10] = vz;
  image->to_real.m[15] = 1.0;

  _identity_mat( &(image->to_voxel) );
  image->to_voxel.m[ 0] = 1.0 / vx;
  image->to_voxel.m[ 5] = 1.0 / vy;
  image->to_voxel.m[10] = 1.0 / vz;
  image->to_voxel.m[15] = 1.0;

  image->geometry = _BAL_HOMOTHETY_GEOMETRY_;

  image->qform_code = 0;
  image->sform_code = 0;

  return( 1 );
}





int BAL_CopyImageGeometry( bal_image *theIm, bal_image *resIm )
{
    resIm->vx = theIm->vx;
    resIm->vy = theIm->vy;
    resIm->vz = theIm->vz;

    resIm->geometry = theIm->geometry;
    _copy_mat( &(theIm->to_real), &(resIm->to_real) );
    _copy_mat( &(theIm->to_voxel), &(resIm->to_voxel) );

    resIm->qform_code = theIm->qform_code;
    resIm->sform_code = theIm->sform_code;
    _copy_mat( &(theIm->sform_to_real), &(resIm->sform_to_real) );

    return( 1 );
}



int BAL_ResizeImageGeometry( bal_image *theIm, bal_image *resIm )
{
  char *proc = "BAL_ResizeImageGeometry";
  _MATRIX theFROMres;
  int i;

  _init_mat( &theFROMres );

  switch( theIm->geometry ) {
  case _BAL_UNKNOWN_GEOMETRY_ :
  case _BAL_HOMOTHETY_GEOMETRY_ :
    if ( BAL_SetImageVoxelSizes( resIm, resIm->vx, resIm->vy, resIm->vz ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to set image voxel sizes\n", proc );
      return(-1);
    }
    resIm->geometry = _BAL_HOMOTHETY_GEOMETRY_;
    break;

  case _BAL_TRANSLATION_GEOMETRY_ :
  case _BAL_QFORM_GEOMETRY_ :
    if ( _alloc_mat( &theFROMres, 4, 4 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: matrix #1 allocation error\n", proc );
      return( -1 );
    }

    for ( i=0; i<16; i++ )
      theFROMres.m[i] = 0.0;

    /* compute voxel to voxel conversion matrix
     */
    theFROMres.m[ 0] = resIm->vx/theIm->vx;
    theFROMres.m[ 5] = resIm->vy/theIm->vy;
    theFROMres.m[10] = resIm->vz/theIm->vz;
    theFROMres.m[15] = 1.0;
    theFROMres.m[ 3] = (resIm->vx - theIm->vx) / (2 * theIm->vx);
    theFROMres.m[ 7] = (resIm->vy - theIm->vy) / (2 * theIm->vy);
    theFROMres.m[11] = (resIm->vz - theIm->vz) / (2 * theIm->vz);

    _mult_mat( &(theIm->to_real), &theFROMres, &(resIm->to_real) );

    if ( InverseMat4x4( resIm->to_real.m, resIm->to_voxel.m ) != 4 ) {
      _free_mat( &theFROMres );
      if ( _verbose_ )
        fprintf( stderr, "%s: matrix inversion error\n", proc );
      return( -1 );
    }
    resIm->geometry = theIm->geometry;
    resIm->qform_code = theIm->qform_code;

    if ( theIm->sform_code ) {
      _mult_mat( &(theIm->sform_to_real), &theFROMres, &(resIm->sform_to_real) );
    }
    resIm->sform_code = theIm->sform_code;

    _free_mat( &theFROMres );
    break;
  }

  return( 1 );
}





/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image structure
 *
 --------------------------------------------------*/





int BAL_InitImage( bal_image *image, char *name,
                   int dimx, int dimy, int dimz, int dimv, bufferType type )
{
  char *proc = "BAL_InitImage";

  if ( image == NULL ) return( -1 );

  /*
  if ( dimx <= 0 || dimy <= 0 || dimz <= 0 || dimv <= 0 )
    return( -1 );
  */

  if ( name != NULL && name[0] != '\0' && strlen(name) > 0 ) {
    image->name = (char*)vtmalloc( (strlen(name)+1)*sizeof( char ),
                                   "image->name", proc );
    if ( image->name == NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: can not allocate name\n", proc );
      return( -1 );
    }
    strcpy( image->name, name );
  }
  else
    image->name = NULL;

  image->type = type;
  
  image->ncols = dimx;
  image->nrows = dimy;
  image->nplanes = dimz;
  image->vdim = dimv;
  
  image->data = NULL;
  image->array = NULL;

  /* image geometry
   */
  BAL_InitImageGeometry( image );

  return( 1 );
}





int BAL_InitFullImage( bal_image *image, char *name,
                   int dimx, int dimy, int dimz, int dimv,
                   typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz,
                   bufferType type )
{
  char *proc = "BAL_InitImage";

  if ( BAL_InitImage( image, name, dimx, dimy, dimz, dimv, type ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to initialize image structure\n", proc );
  }
  if ( BAL_SetImageVoxelSizes( image, vx, vy, vz ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to set image voxel sizes\n", proc );
  }

  return( 1 );
}





int BAL_InitImageFromImage( bal_image *image, char *name,
                            bal_image *from, bufferType type )
{

  if ( BAL_InitImage( image, name, from->ncols, from->nrows, from->nplanes, from->vdim, type ) != 1 )
    return( -1 );

  if ( BAL_AllocImageGeometry( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }

  if ( BAL_CopyImageGeometry( from, image ) != 1 ) {
      BAL_FreeImage( image );
      return( -1 );
  }

  return( 1 );
}





int BAL_InitScalarImageFromImage( bal_image *image, char *name,
                                  bal_image *from, bufferType type )
{

  if ( BAL_InitImage( image, name, from->ncols, from->nrows, from->nplanes, from->vdim, type ) != 1 )
    return( -1 );

  image->vdim = 1;

  if ( BAL_AllocImageGeometry( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }

  if ( BAL_CopyImageGeometry( from, image ) != 1 ) {
      BAL_FreeImage( image );
      return( -1 );
  }


  return( 1 );
}





/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image allocation
 *
 --------------------------------------------------*/





void BAL_FreeImage( bal_image *image )
{
  char *proc = "BAL_FreeImage";

  if ( image->array != NULL ) vtfree( image->array );
  image->array = NULL;
  if ( image->data != NULL ) vtfree( image->data );
  image->data = NULL;
  if ( image->name != NULL ) vtfree( image->name );
  image->name = NULL;

  BAL_FreeImageGeometry( image );

  if ( BAL_InitImage( image, (char *)NULL, 0, 0, 0, 0, TYPE_UNKNOWN ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error\n", proc );
  }
}





static int _AllocArrayImage( bal_image *image )
{
  char *proc = "_AllocArrayImage";
  size_t size=0;
  
  if ( image == NULL || image->data == NULL ) 
    return( -1 );
  if ( image->ncols * image->nrows * image->nplanes * image->vdim <= 0 )
    return( -1 );
  
  switch ( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    size += image->nplanes * sizeof( u8 **);
    size += image->nplanes * image->nrows * sizeof( u8*);
    break;
  case SCHAR :
    size += image->nplanes * sizeof( s8 **);
    size += image->nplanes * image->nrows * sizeof( s8 *);
    break;
  case USHORT :
    size += image->nplanes * sizeof( u16 **);
    size += image->nplanes * image->nrows * sizeof( u16 *);
    break;
  case SSHORT :
    size += image->nplanes * sizeof( s16 **);
    size += image->nplanes * image->nrows * sizeof( s16 *);
    break;
  case UINT :
    size += image->nplanes * sizeof( u32 **);
    size += image->nplanes * image->nrows * sizeof( u32 *);
    break;
  case SINT :
    size += image->nplanes * sizeof( s32 **);
    size += image->nplanes * image->nrows * sizeof( s32 *);
    break;
  case FLOAT :
    size += image->nplanes * sizeof( r32 **);
    size += image->nplanes * image->nrows * sizeof( r32 *);
    break;
  case DOUBLE :
    size += image->nplanes * sizeof( r64 **);
    size += image->nplanes * image->nrows * sizeof( r64 *);
    break;  }
  
  image->array = (void***)vtmalloc( size, "image->array", proc );
  if ( image->array == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  return( 1 );
}



static int _BuildArrayImage( bal_image *image )
{  
  char *proc = "_BuildArrayImage";
  size_t k, j;

  if ( image == NULL || image->data == NULL || image->array == NULL ) 
    return( -1 );
  if ( image->ncols * image->nrows * image->nplanes * image->vdim <= 0 )
    return( -1 );
  
  switch ( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      unsigned char ***z = (unsigned char ***)image->array;
      unsigned char **zy = (unsigned char **)(z + image->nplanes);
      unsigned char *zyx = (unsigned char *)(image->data);
      for ( k=0; k<image->nplanes; k++, z++ ) {
        *z = zy;
        for ( j=0 ; j < image->nrows; j++, zy++, zyx += image->ncols*image->vdim )
          *zy = zyx;
      }
    }
    break;
  case USHORT :
    {
      unsigned short int ***z = (unsigned short int ***)image->array;
      unsigned short int **zy = (unsigned short int **)(z + image->nplanes);
      unsigned short int *zyx = (unsigned short int *)(image->data);
      for ( k=0; k<image->nplanes; k++, z++ ) {
        *z = zy;
        for ( j=0 ; j < image->nrows; j++, zy++, zyx += image->ncols*image->vdim )
          *zy = zyx;
      }
    }
    break;
  case SSHORT :
    {
      short int ***z = (short int ***)image->array;
      short int **zy = (short int **)(z + image->nplanes);
      short int *zyx = (short int *)(image->data);
      for ( k=0; k<image->nplanes; k++, z++ ) {
        *z = zy;
        for ( j=0 ; j < image->nrows; j++, zy++, zyx += image->ncols*image->vdim )
          *zy = zyx;
      }
    }
    break;
  case FLOAT :
    {
      float ***z = (float ***)image->array;
      float **zy = (float **)(z + image->nplanes);
      float *zyx = (float *)(image->data);
      for ( k=0; k<image->nplanes; k++, z++ ) {
        *z = zy;
        for ( j=0 ; j < image->nrows; j++, zy++, zyx += image->ncols*image->vdim )
          *zy = zyx;
      }
    }
    break;
  }
  return( 1 );
}



int BAL_AllocArrayImage( bal_image *image )
{
    char *proc = "BAL_AllocArrayImage";

    if ( _AllocArrayImage( image ) != 1 ) {
        if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to allocate image array\n", proc );
        }
        return( -1 );
    }
    if ( _BuildArrayImage( image ) != 1 ) {
      if ( image->array != NULL ) vtfree( image->array );
      image->array = NULL;
      if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to build image array\n", proc );
      }
      return( -1 );
    }
    return( 1 );
}





int BAL_AllocImage( bal_image *image )
{
  char *proc = "BAL_AllocImage";
  size_t size=0;
  
  if ( image == NULL ) 
    return( -1 );
  if ( image->ncols * image->nrows * image->nplanes * image->vdim <= 0 )
    return( -1 );
  
  if ( image->data != (void*)NULL || image->array != (void***)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: warning, image to be allocated seems already allocated\n", proc );
  }

  image->data = image->array = NULL;

  size = BAL_ImageDataSize( image );
  if ( size <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad image size\n", proc );
    return( -1 );
  }

  if ( _debug_ ) 
    fprintf( stderr, "%s: allocate image of size %lu\n", proc, size );

  image->data = (void*)vtmalloc( size, "image->data", proc );
  if ( image->data == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  memset( image->data, 0, size );
  
  if ( _debug_ ) 
    fprintf( stderr, "%s: allocation done\n", proc );

  if ( BAL_AllocArrayImage( image ) != 1 ) {
      BAL_FreeImage( image );
      if (_verbose_ )
          fprintf( stderr, "%s: unable to build image array\n", proc );
      return( -1 );
  }

  return( 1 );
}





int BAL_AllocFullImage( bal_image *image, char *name,
                        int dimx, int dimy, int dimz, int dimv,
                        typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz,
                        bufferType type )
{
  if ( BAL_InitFullImage( image, name, dimx, dimy, dimz, dimv, vx, vy, vz, type ) != 1 )
    return( -1 );
  if ( BAL_AllocImage( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }
  return( 1 );
}





int BAL_AllocImageFromImage( bal_image *image, char *name,
                                 bal_image *from, bufferType type )
{

  if ( BAL_InitImageFromImage( image, name, from, type ) != 1 ) {
    return( -1 );
  }

  if ( BAL_AllocImage( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }
  return( 1 );
}





int BAL_AllocScalarImageFromImage( bal_image *image, char *name,
                                 bal_image *from, bufferType type )
{

  if ( BAL_InitImageFromImage( image, name, from, type ) != 1 ) {
    return( -1 );
  }

  image->vdim = 1;

  if ( BAL_AllocImage( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }
  return( 1 );
}





/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT:
 *
 --------------------------------------------------*/





int BAL_CopyImage( bal_image *theIm, bal_image *resIm )
{
  bufferType theType = theIm->type;
  bufferType resType = resIm->type;
 
  if ( theIm->ncols != resIm->ncols 
       || theIm->nrows != resIm->nrows 
       || theIm->nplanes != resIm->nplanes 
       || theIm->vdim != resIm->vdim ) {
    if ( _verbose_ ) 
      fprintf( stderr, "BAL_CopyImage: image have different dimensions\n" );
    return( -1 );
  }
  if ( theType == TYPE_UNKNOWN || resType == TYPE_UNKNOWN ) {
    if ( _verbose_ )
      fprintf( stderr, "BAL_CopyImage: unable to deal with such image types\n" );
    return( -1 );
  }

  if ( ConvertBuffer( theIm->data, theType, resIm->data, resType, 
                      theIm->ncols * theIm->nrows * theIm->nplanes * theIm->vdim ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "BAL_CopyImage: unable to convert such image type\n" );
    return( -1 ); 
  }

  if ( BAL_CopyImageGeometry( theIm, resIm ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "BAL_CopyImage: unable to copy image geometry\n" );
      return( -1 );
  }

  return( 1 );
       
}



/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/

int BAL_FillImage( bal_image *theIm, float fval )
{
  char *proc = "BAL_FillImage";
  int ival;
  size_t i, v;

  v = theIm->ncols * theIm->nrows * theIm->nplanes * theIm->vdim;

  switch( theIm->type ) {
  default :
    if (_verbose_)
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 *theBuf = (u8*)theIm->data;
      if ( fval < 0.0 ) ival = 0;
      else if ( fval > 255.0 ) ival = 255;
      else ival = (int)(fval + 0.5);
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = (u8)ival;
    }
    break;
  case SSHORT :
    {
      s16 *theBuf = (s16*)theIm->data;
      if ( fval < -32768.0 ) ival = -32768;
      else if ( fval > 32767.0 ) ival = 32767;
      else ival = (fval > 0 ) ? (int)(fval + 0.5) : (int)(fval - 0.5);
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = (s16)ival;
    }
    break;
  case USHORT :
    {
      u16 *theBuf = (u16*)theIm->data;
      if ( fval < 0.0 ) ival = 0;
      else if ( fval > 65535.0 ) ival = 65535;
      else ival = (int)(fval + 0.5);
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = (u16)ival;
    }
    break;
  case FLOAT :
    {
      r32 *theBuf = (r32*)theIm->data;
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = fval;
    }
    break;
  }
  return( 1 );
}



int BAL_NormaImage( bal_image *theIm, bal_image *resIm )
{
  bufferType theType = theIm->type;
  bufferType resType = resIm->type;
  
  double a, v, theMin, theMax;
  int vol;

  if ( theIm->ncols != resIm->ncols 
       || theIm->nrows != resIm->nrows 
       || theIm->nplanes != resIm->nplanes 
       || theIm->vdim != resIm->vdim ) {
    if ( _verbose_ ) 
      fprintf( stderr, "BAL_NormaImage: image have different dimensions\n" );
    return( -1 );
  }

  if ( theType == resType ) 
    return( 1 );

  if ( theType == TYPE_UNKNOWN || resType == TYPE_UNKNOWN ) {
    if ( _verbose_ ) 
      fprintf( stderr, "BAL_NormaImage: unable to deal with such image types\n" );
    return( -1 );
  }

  /* computation of minimum and maximum values
   */
  vol = theIm->ncols * theIm->nrows * theIm->nplanes;
  switch ( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "BAL_NormaImage: input image type not handled yet\n" );
    return( -1 );
  case UCHAR :
    {
      unsigned char *theBuf = (unsigned char *)theIm->data;
      theMin = theMax = *theBuf; vol--;
      while ( vol-- > 0 ) {
        v = (double)(*theBuf++);
        if ( theMin > v ) theMin = v;
        if ( theMax < v ) theMax = v;
      }
    }
    break;
  case USHORT :
    {
      unsigned short int *theBuf = (unsigned short int *)theIm->data;
      theMin = theMax = *theBuf; vol--;
      while ( vol-- > 0 ) {
        v = (double)(*theBuf++);
        if ( theMin > v ) theMin = v;
        if ( theMax < v ) theMax = v;
      }
    }
    break;
  case SSHORT :
    {
      short int *theBuf = (short int *)theIm->data;
      theMin = theMax = *theBuf; vol--;
      while ( vol-- > 0 ) {
        v = (double)(*theBuf++);
        if ( theMin > v ) theMin = v;
        if ( theMax < v ) theMax = v;
      }
    }
    break;
  }

  if ( theMax <= theMin ) {
    if ( _verbose_ )
      fprintf( stderr, "BAL_NormaImage: constant input image\n" );
    return( -1 );
  }
  

  /* normalisation
   */
  vol = theIm->ncols * theIm->nrows * theIm->nplanes;
  switch ( resIm->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "BAL_NormaImage: output image type not handled yet in normalisation step\n" );
    return( -1 );
  case UCHAR :
    {
      unsigned char *resBuf = (unsigned char *)resIm->data;
      a = ((double)255 - (double)0)/(theMax-theMin);
      switch( theIm->type ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "BAL_NormaImage: input image type not handled yet in normalisation step\n" );
      return( -1 );
      case USHORT :
        {
          unsigned short int *theBuf = (unsigned short int *)theIm->data;
          while ( vol-- > 0 )
            *resBuf++ = (unsigned char)( a*(((double)(*theBuf++))-theMin) + (double)0 + 0.5);
        }
        break;
      case SSHORT :
        {
          short int *theBuf = (short int *)theIm->data;
          while ( vol-- > 0 )
            *resBuf++ = (unsigned char)( a*(((double)(*theBuf++))-theMin) + (double)0 + 0.5);
        }
        break;
      }
    }
    /* end of case UCHAR */
    break;
  }

  return( 1 );
       
}





size_t BAL_ImageDataSize( bal_image *image )
{
  size_t size=0;
  
  if ( image == NULL ) 
    return( -1 );
  if ( image->ncols * image->nrows * image->nplanes * image->vdim <= 0 )
    return( -1 );
  
  switch ( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "BAL_ImageDataSize: image type not handled yet\n" );
    return( -1 );
  case UCHAR :
    size = sizeof( u8 );
    break;
  case SCHAR :
    size = sizeof( s8 );
    break;
  case USHORT :
    size = sizeof( u16 );
    break;
  case SSHORT :
    size = sizeof( s16 );
    break;
  case UINT :
    size = sizeof( u32 );
    break;
  case SINT :
    size = sizeof( s32 );
    break;
  case FLOAT :
    size = sizeof( r32 );
    break;
  case DOUBLE :
    size = sizeof( r64 );
    break;
  }
  size *= image->ncols * image->nrows * image->nplanes * image->vdim;
  return( size );
}






/*--------------------------------------------------
 *
 * Image value
 *
 --------------------------------------------------*/

double BAL_GetXYKvalue( bal_image *image, double x, double y, int k ) {
  char *proc = "BAL_GetXYKvalue";
  int ix, iy;
  int dimx;
  int dimy;
  int dimz;
  double res = 0.0;
  double dx, dy;
  
  if ( image == (bal_image*)NULL ) return( 0.0 );
  if ( image->data == NULL || image->array == NULL ) return( 0.0 );

  dimx = image->ncols;
  dimy = image->nrows;
  dimz = image->nplanes;

  if ( dimx == 0 || dimy == 0 || dimz == 0  ) return( 0.0 );

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

  if ( ix >= dimx - 1 ) {
    ix = dimx-2;
    dx = 1;
  }
  if ( iy >= dimy - 1 ) {
    iy = dimy-2;
    dy = 1;
  }

  switch( image->type ) {
    
  case FLOAT :
    {
      float ***buf = (float***)image->array;
      res = (1-dx) * ( (1-dy) * buf[k][iy][ix] 
                       + (dy) * buf[k][iy+1][ix] )
            + (dx) * ( (1-dy) * buf[k][iy][ix+1]
                       + (dy) * buf[k][iy+1][ix+1] );
    }
    break;
    
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such image type not handled in the switch\n", proc );
    return( (double)0.0 );
  }

  return( res );

}



double BAL_GetXYZvalue( bal_image *image, double x, double y, double z ) {
  char *proc = "BAL_GetXYZvalue";
  int ix, iy, iz;
  int dimx;
  int dimy;
  int dimz;
  double res = 0.0;
  double dx, dy, dz;
  
  if ( image == (bal_image*)NULL ) return( 0.0 );
  if ( image->data == NULL || image->array == NULL ) return( 0.0 );

  dimx = image->ncols;
  dimy = image->nrows;
  dimz = image->nplanes;

  if ( dimx == 0 || dimy == 0 || dimz == 0  ) return( 0.0 );

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

  if ( dimx == 1 ) {
      ix = 0;
      dx = 0;
  }
  else if ( ix >= dimx - 1 ) {
    ix = dimx-2;
    dx = 1;
  }
  if ( dimy == 1 ) {
      iy = 0;
      dy = 0;
  }
  else if ( iy >= dimy - 1 ) {
    iy = dimy-2;
    dy = 1;
  }
  if ( dimz == 1 ) {
      iz = 0;
      dz = 0;
  }
  else if ( iz >= dimz - 1 ) {
    iz = dimz-2;
    dz = 1;
  }

  switch( image->type ) {
    
  case UCHAR :
    {
      u8 ***buf = (u8***)image->array;
      if ( dimz == 1 ) {
          res = (1-dx) * ( (1-dy) * buf[iz][iy][ix]
                           + (dy) * buf[iz][iy+1][ix] )
                + (dx) * ( (1-dy) * buf[iz][iy][ix+1]
                           + (dy) * buf[iz][iy+1][ix+1] );
      }
      else {
          res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                      + (dz)*buf[iz+1][iy][ix] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                      + (dz)*buf[iz+1][iy+1][ix] ) )
                + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                      + (dz)*buf[iz+1][iy][ix+1] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                      + (dz)*buf[iz+1][iy+1][ix+1] ) );
      }
    }
    break;

  case SCHAR :
    {
      s8 ***buf = (s8***)image->array;
      if ( dimz == 1 ) {
          res = (1-dx) * ( (1-dy) * buf[iz][iy][ix]
                           + (dy) * buf[iz][iy+1][ix] )
                + (dx) * ( (1-dy) * buf[iz][iy][ix+1]
                           + (dy) * buf[iz][iy+1][ix+1] );
      }
      else {
          res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                      + (dz)*buf[iz+1][iy][ix] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                      + (dz)*buf[iz+1][iy+1][ix] ) )
                + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                      + (dz)*buf[iz+1][iy][ix+1] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                      + (dz)*buf[iz+1][iy+1][ix+1] ) );
      }
    }
    break;

  case USHORT :
    {
      u16 ***buf = (u16***)image->array;
      if ( dimz == 1 ) {
          res = (1-dx) * ( (1-dy) * buf[iz][iy][ix]
                           + (dy) * buf[iz][iy+1][ix] )
                + (dx) * ( (1-dy) * buf[iz][iy][ix+1]
                           + (dy) * buf[iz][iy+1][ix+1] );
      }
      else {
          res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                      + (dz)*buf[iz+1][iy][ix] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                      + (dz)*buf[iz+1][iy+1][ix] ) )
                + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                      + (dz)*buf[iz+1][iy][ix+1] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                      + (dz)*buf[iz+1][iy+1][ix+1] ) );
      }
    }
    break;

  case SSHORT :
    {
      s16 ***buf = (s16***)image->array;
      if ( dimz == 1 ) {
          res = (1-dx) * ( (1-dy) * buf[iz][iy][ix]
                           + (dy) * buf[iz][iy+1][ix] )
                + (dx) * ( (1-dy) * buf[iz][iy][ix+1]
                           + (dy) * buf[iz][iy+1][ix+1] );
      }
      else {
          res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                      + (dz)*buf[iz+1][iy][ix] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                      + (dz)*buf[iz+1][iy+1][ix] ) )
                + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                      + (dz)*buf[iz+1][iy][ix+1] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                      + (dz)*buf[iz+1][iy+1][ix+1] ) );
      }
    }
    break;

  case FLOAT :
    {
      float ***buf = (float***)image->array;
      if ( dimz == 1 ) {
          res = (1-dx) * ( (1-dy) * buf[iz][iy][ix]
                           + (dy) * buf[iz][iy+1][ix] )
                + (dx) * ( (1-dy) * buf[iz][iy][ix+1]
                           + (dy) * buf[iz][iy+1][ix+1] );
      }
      else {
          res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix]
                                      + (dz)*buf[iz+1][iy][ix] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix]
                                      + (dz)*buf[iz+1][iy+1][ix] ) )
                + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1]
                                      + (dz)*buf[iz+1][iy][ix+1] )
                           + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
                                      + (dz)*buf[iz+1][iy+1][ix+1] ) );
      }
    }
    break;
    
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such image type not handled in the switch\n", proc );
    return( (double)0.0 );
  }

  return( res );

}



/*--------------------------------------------------
 *
 * IMAGE I/O
 *
 --------------------------------------------------*/

static void _PrintImageStatistics( FILE *f, bal_image *image )
{
  size_t i;
  size_t v = image->ncols * image->nrows * image->nplanes * image->vdim;
  double m;

  if ( image->data == NULL || v <= 0 ) return;

  switch ( image->type ) {
  default :                
    fprintf( f, "TYPE_UNKNOWN\n" );
    break;
  case UCHAR :  
    {
      unsigned char *buf = (unsigned char *)image->data;
      unsigned char min, max;
      min = max = buf[0];
      m = buf[0];
      for ( i=1; i<v; i++ ) {
        if ( min > buf[i] ) min = buf[i];
        if ( max < buf[i] ) max = buf[i];
        m += buf[i];
      }
      fprintf( f, "  - min = %d , max = %d , mean = %f\n", min, max, m/(double)v );
    }
    break;
  case USHORT : 
    {
      unsigned short int *buf = (unsigned short int *)image->data;
      unsigned short int min, max;
      min = max = buf[0];
      m = buf[0];
      for ( i=1; i<v; i++ ) {
        if ( min > buf[i] ) min = buf[i];
        if ( max < buf[i] ) max = buf[i];
        m += buf[i];
      }
      fprintf( f, "  - min = %d , max = %d , mean = %f\n", min, max, m/(double)v );
    }
    break;
  case SSHORT :   
    {
      short int *buf = (short int *)image->data;
      short int min, max;
      min = max = buf[0];
      m = buf[0];
      for ( i=1; i<v; i++ ) {
        if ( min > buf[i] ) min = buf[i];
        if ( max < buf[i] ) max = buf[i];
        m += buf[i];
      }
      fprintf( f, "  - min = %d , max = %d , mean = %f\n", min, max, m/(double)v );
    }
    break;
  case FLOAT :
    {
      float *buf = (float *)image->data;
      float min, max;
      min = max = buf[0];
      m = buf[0];
      for ( i=1; i<v; i++ ) {
        if ( min > buf[i] ) min = buf[i];
        if ( max < buf[i] ) max = buf[i];
        m += buf[i];
      }
      fprintf( f, "  - min = %f , max = %f , mean = %f\n", min, max, m/(double)v );
    }
    break;
  }  
}





void BAL_PrintImage( FILE *f, bal_image *image, char *s )
{
  if ( s != (char *)NULL )
    fprintf( f, "'%s' information:\n", s );
  else {
    if ( image->name != (char *)NULL )
      fprintf( f, "'%s' information:\n", image->name );
    else
      fprintf( f, "image information:\n" );
  }

  fprintf( f, "  - dimensions [x y z] = %lu %lu %lu\n", 
           image->ncols, image->nrows, image->nplanes );
  if ( image->vdim > 1 )
    fprintf( f, "  - vectorial dimension [v] = %lu\n", image->vdim );
  fprintf( f, "  - voxel size [x y z] = %f %f %f\n",
           image->vx, image->vy, image->vz );
  fprintf( f, "  - image type is: " );
  switch ( image->type ) {
  default :           fprintf( f, "not handled\n" ); break;
  case TYPE_UNKNOWN : fprintf( f, "TYPE_UNKNOWN\n" ); break;
  case UCHAR :        fprintf( f, "UCHAR\n" ); break;
  case USHORT :       fprintf( f, "USHORT\n" ); break;
  case UINT :         fprintf( f, "UINT\n" ); break;
  case SCHAR :        fprintf( f, "SCHAR\n" ); break;
  case SSHORT :       fprintf( f, "SSHORT\n" ); break;
  case SINT :         fprintf( f, "SINT\n" ); break;
  case ULINT :        fprintf( f, "ULINT\n" ); break;
  case FLOAT :        fprintf( f, "FLOAT\n" ); break;
  case DOUBLE :       fprintf( f, "DOUBLE\n" ); break;
  }  

  fprintf( f, "  - image geometry is: " );
  switch ( image->geometry ) {
  default :           fprintf( f, "not handled\n" ); break;
  case _BAL_UNKNOWN_GEOMETRY_ : fprintf( f, "_BAL_UNKNOWN_GEOMETRY_\n" ); break;
  case _BAL_HOMOTHETY_GEOMETRY_ : fprintf( f, "_BAL_HOMOTHETY_GEOMETRY_\n" ); break;
  case _BAL_TRANSLATION_GEOMETRY_ : fprintf( f, "_BAL_TRANSLATION_GEOMETRY_\n" ); break;
  case _BAL_QFORM_GEOMETRY_ :   fprintf( f, "_BAL_QFORM_GEOMETRY_\n" ); break;
  }
  fprintf( f, "  - conversion matrix voxel to real:\n" );
  _print_mat ( f, &(image->to_real), (char *)NULL );
  fprintf( f, "  - conversion matrix real to voxel:\n" );
  _print_mat ( f, &(image->to_voxel), (char *)NULL );

  fprintf( f, "  - qform_code: %d\n", image->qform_code );
  fprintf( f, "  - sform_code: %d\n", image->sform_code );
  if ( image->sform_code ) {
    fprintf( f, "  - SForm conversion matrix voxel to real:\n" );
    _print_mat ( f, &(image->sform_to_real), (char *)NULL );
  }

  if ( _verbose_ >= 4 || _debug_ ) {
    fprintf( f, "  - image buffer address is %p\n", image->data );
    fprintf( f, "      \"   array     \"    is %p\n", image->array );
  }

  if ( image->data != NULL ) {
    _PrintImageStatistics( f, image );
  }
  else {
    fprintf( f, "- empty image\n" );
  }
}





void BAL_PrintParImage( FILE *f, bal_image *image, char *s )
{

  if ( s != (char *)NULL )
    fprintf( f, "%s: ", s );

  if ( image->vdim > 1 )
    fprintf( f, "-v %lu ", image->vdim );

  fprintf( f, "-x %lu ", image->ncols);
  fprintf( f, "-y %lu ", image->nrows );
  fprintf( f, "-z %lu ", image->nplanes);

  switch ( image->type ) {
  case SCHAR :
    fprintf( f, "-f -o 1 -s ");
    break;
  case UCHAR :
    fprintf( f, "-f -o 1 ");
    break;
  case SSHORT :
    fprintf( f, "-f -o 2 -s ");
    break;
  case USHORT :
    fprintf( f, "-f -o 2 ");
    break;
  case SINT :
    fprintf( f, "-f -o 4 -s ");
    break;
  case UINT :
    fprintf( f, "-f -o 4 ");
    break;
  case ULINT :
    fprintf( f, "-r -o 8 ");
    break;
  case FLOAT :
    fprintf( f, "-r -o 4 ");
    break;
  case DOUBLE :
    fprintf( f, "-r -o 8 ");
    break;
  case TYPE_UNKNOWN :
  default :
    fprintf( f, "unknown type ");
    break;
  }

  fprintf( f, "-vx %f ", image->vx );
  fprintf( f, "-vy %f ", image->vy );
  fprintf( f, "-vz %f ", image->vz );

  fprintf( f, "\n" );
}







int BAL_ReadImage( bal_image *image, char *name, int normalisation )
{
  char *proc = "BAL_ReadImage";
  _image *theIm;
  bufferType type=UCHAR;
  size_t size;
  int i, j, k;
  
  if ( _verbose_ >= 2 ) 
    fprintf( stderr, "%s: will read '%s'\n", proc, name );


  if ( name == (char*) NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no image name\n", proc );
    return( -1 );
  }

  theIm = _readImage( name );
  if ( theIm == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read '%s'\n", proc, name );
    return( -1 );
  }

  
  switch( theIm->wordKind ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such word kind not handled yet\n", proc );
    _freeImage( theIm );
    return( -1 );

  case WK_FLOAT :
    switch( theIm->wdim ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such float word dim not handled yet\n", proc );
      _freeImage( theIm );
      return( -1 );
    case 4 :
      type = FLOAT;
      break;
    case 8 :
      type = DOUBLE;
      break;
    }
    break;

  case WK_FIXED :
    switch( theIm->sign ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such sign not handled yet\n", proc );
      _freeImage( theIm );
      return( -1 );
    case SGN_SIGNED :
      switch( theIm->wdim ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such signed word dim not handled yet\n", proc );
        _freeImage( theIm );
        return( -1 );
      case 1 :
        type = SCHAR;
        break;
      case 2 :
        type = SSHORT;
        break;
      case 4 :
        type = SINT;
        break;
      }
      break;
    case SGN_UNSIGNED :
      switch( theIm->wdim ) {
      default :
        if ( _verbose_ )
          fprintf( stderr, "%s: such unsigned word dim not handled yet\n", proc );
        _freeImage( theIm );
        return( -1 );
      case 1 :
        type = UCHAR;
        break;
      case 2 :
        type = USHORT;
        break;
      case 4 :
        type = UINT;
        break;
      }
      break;
    } /* switch( theIm->sign ) */
    break;
  }


  
  if ( normalisation && type != UCHAR ) {

    if ( _verbose_ ) {
      fprintf( stderr, "%s: normalization of image '%s' into unsigned char\n", proc, name );
    }

    if ( BAL_InitImage( image, name, theIm->xdim, theIm->ydim,
                        theIm->zdim, theIm->vdim, UCHAR ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize image '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }

    if ( BAL_AllocImage( image ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate image '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }
    
    if ( ConvertScaleBuffer( theIm->data, type, 
                             image->data, UCHAR,
                             theIm->xdim * theIm->ydim * theIm->zdim * theIm->vdim ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to convert image '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }

  }
  else {

    if ( BAL_InitImage( image, name, theIm->xdim, theIm->ydim,
                        theIm->zdim, theIm->vdim, type ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize image '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }

    if ( BAL_AllocImage( image ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate image '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }

    size = BAL_ImageDataSize( image );
    if ( size <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: can not copy raw data of '%s'\n", proc, name );
      _freeImage( theIm );
      return( -1 );
    }

    (void)memcpy( image->data, theIm->data, size );

  }

  if ( BAL_AllocImageGeometry( image ) != 1 ) {
    BAL_FreeImage( image );
    _freeImage( theIm );
    if ( _verbose_ )
      fprintf( stderr, "%s: can not set initialize geometry of '%s'\n", proc, name );
    return( -1 );
  }

  image->vx = theIm->vx;
  image->vy = theIm->vy;
  image->vz = theIm->vz;

  if ( theIm->qform_code ) {


    for ( k=0, i=0; i<4; i++ )
    for ( j=0; j<4; j++, k++ )
      image->to_real.m[k] = theIm->qform_toreal[i][j];
    image->geometry = _BAL_QFORM_GEOMETRY_;
    if ( InverseMat4x4( image->to_real.m, image->to_voxel.m ) != 4 ) {
      BAL_FreeImage( image );
      _freeImage( theIm );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert qform matrix of '%s'\n", proc, name );
      return( -1 );
    }
  }
  else {
    /*  0  1  2  3
     *  4  5  6  7
     *  8  9 10 11
     * 12 13 14 15
     */
    BAL_SetImageVoxelSizes( image, image->vx, image->vy, image->vz );
    if ( theIm->t_is_set ) {
      image->to_real.m[ 3] = theIm->tx;
      image->to_real.m[ 7] = theIm->ty;
      image->to_real.m[11] = theIm->tz;
      image->to_voxel.m[ 3] = - theIm->tx / theIm->vx;
      image->to_voxel.m[ 7] = - theIm->ty / theIm->vy;
      image->to_voxel.m[11] = - theIm->tz / theIm->vz;
      image->geometry = _BAL_TRANSLATION_GEOMETRY_;
    }
  }

  image->qform_code = theIm->qform_code;
  image->sform_code = theIm->sform_code;

  for ( k=0, i=0; i<4; i++ )
  for ( j=0; j<4; j++, k++ )
    image->sform_to_real.m[k] = theIm->sform_toreal[i][j];


  _freeImage( theIm );
  
  return( 1 );
}



int BAL_WriteImage( bal_image *image, char *name )
{
  char *proc = "BAL_WriteImage";
  _image *theIm = _initImage();
  int i, j, k;
  
  if ( _verbose_ >= 2 ) 
    fprintf( stderr, "%s: will write '%s'\n", proc, name );

  if ( theIm == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate structure for '%s'\n", proc, name );
    return( -1 );
  }

  theIm->xdim = image->ncols;
  theIm->ydim = image->nrows;
  theIm->zdim = image->nplanes;
  theIm->vdim = image->vdim;

  theIm->vx = image->vx;
  theIm->vy = image->vy;
  theIm->vz = image->vz;

  theIm->qform_code = 0;
  if ( image->to_real.l == 4 && image->to_real.c == 4 && image->to_real.m != (double*)NULL ) {
    if ( image->geometry == _BAL_QFORM_GEOMETRY_ ) {
      for ( k=0, i=0; i<4; i++ )
      for ( j=0; j<4; j++, k++ )
        theIm->qform_toreal[i][j] = image->to_real.m[k];
      theIm->qform_code = image->qform_code;
    }
    else if ( image->geometry == _BAL_TRANSLATION_GEOMETRY_ ) {
      theIm->tx = image->to_real.m[ 3];
      theIm->ty = image->to_real.m[ 4];
      theIm->tz = image->to_real.m[11];
      theIm->t_is_set = 1;
    }
  }
  else {
    if ( image->qform_code != 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird, qform_code is %d but qform matrix is not well formed\n", proc, image->qform_code );
    }
  }

  theIm->sform_code = image->sform_code;
  if ( image->sform_to_real.l == 4 && image->sform_to_real.c == 4 && image->sform_to_real.m != (double*)NULL ) {
    for ( k=0, i=0; i<4; i++ )
    for ( j=0; j<4; j++, k++ )
      theIm->sform_toreal[i][j] = image->sform_to_real.m[k];
  }
  else {
    if ( image->sform_code != 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird, sform_code is %d but sform matrix is not well formed\n", proc, image->sform_code );
    }
  }

  switch( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    _freeImage( theIm );
    return( -1 );
  case UCHAR :
    theIm->wdim = 1;
    theIm->wordKind = WK_FIXED;
    theIm->sign = SGN_UNSIGNED;
    break;
  case USHORT :
    theIm->wdim = 2;
    theIm->wordKind = WK_FIXED;
    theIm->sign = SGN_UNSIGNED;
    break;
  case SSHORT :
    theIm->wdim = 2;
    theIm->wordKind = WK_FIXED;
    theIm->sign = SGN_SIGNED;
    break;
  case FLOAT :
    theIm->wdim = sizeof( float );
    theIm->wordKind = WK_FLOAT;
    theIm->sign = SGN_UNKNOWN;
    break;
  }
  
  theIm->data = image->data;

  if ( _writeImage( theIm, name ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write '%s'\n", proc, name );
    theIm->data = NULL;
    _freeImage( theIm );
    return( -1 );
  }
  
  theIm->data =NULL;
  _freeImage( theIm );
  return( 1 );
}











/*--------------------------------------------------
 *
 * IMAGE FILTERING
 *
 --------------------------------------------------*/
static filterType theFilter = GAUSSIAN_YOUNG_2002;

void BAL_SetFilterType( filterType filter )
{
  theFilter = filter;
}



int BAL_SmoothImage( bal_image *theIm,
                     bal_doublePoint *theSigma )
{
  return( BAL_SmoothImageIntoImage( theIm, theIm, theSigma ) );
}



int BAL_SmoothImageIntoImage( bal_image *theIm, bal_image *resIm,
                              bal_doublePoint *theSigma )
{
  char *proc = "BAL_SmoothImageIntoImage";
  bufferType theType = TYPE_UNKNOWN;
  bufferType resType = TYPE_UNKNOWN;
  int theDim[3];
  int borders[3] = {0, 0, 0};
  typeFilteringCoefficients filter[3];

  

  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );

  theType = theIm->type;
  resType = resIm->type;
  if ( theType == TYPE_UNKNOWN || resType == TYPE_UNKNOWN ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with such image type\n", proc );
    return( -1 );
  }
  
  if ( theIm->ncols != resIm->ncols || theIm->nrows != resIm->nrows
       || theIm->nplanes != resIm->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: images should have the same dimensions\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;

  if ( theIm->ncols > 1 && theSigma->x > 0.0 ) {
    filter[0].type = theFilter;
    filter[0].derivative = SMOOTHING;
    filter[0].coefficient = theSigma->x;
    borders[0] = (int)floor( theSigma->x );
  }
  else {
    filter[0].derivative = NODERIVATIVE;
  }

  if ( theIm->nrows > 1 && theSigma->y > 0.0 ) {
    filter[1].type = theFilter;
    filter[1].derivative = SMOOTHING;
    filter[1].coefficient = theSigma->y;
    borders[1] = (int)floor( theSigma->y );
  }
  else {
    filter[1].derivative = NODERIVATIVE;
  }

  if ( theIm->nplanes > 1 && theSigma->z > 0.0 ) {
    filter[2].type = theFilter;
    filter[2].derivative = SMOOTHING;
    filter[2].coefficient = theSigma->z;
    borders[2] = (int)floor( theSigma->z );
  }
  else {
    filter[2].derivative = NODERIVATIVE;
  }

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)resIm->data, resType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to smooth image\n", proc );
    return( -1 );
  }

  /*
  if ( RecursiveFilterOnBuffer( (void*)theIm->data, theType,
                                (void*)resIm->data, resType,
                                theDim, borders, derivatives,
                                coeffs, theFilter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to smooth image\n", proc );
    return( -1 );
  }
  */

  return( 1 );
}



int BAL_2DDerivativesOfImage( bal_image *theIm, 
                              bal_image *theDx, bal_image *theDy,
                              bal_doublePoint *theSigma )
{
  char *proc = "BAL_2DDerivativesOfImage";
  bufferType theType = TYPE_UNKNOWN;
  bufferType dxType = TYPE_UNKNOWN;
  bufferType dyType = TYPE_UNKNOWN;
  int theDim[3];
  int borders[3] = {0, 0, 0};
  typeFilteringCoefficients filter[3];

  

  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );

  theType = theIm->type;
  dxType = theDx->type;
  dyType = theDy->type;
  if ( theType == TYPE_UNKNOWN || dxType == TYPE_UNKNOWN || dyType == TYPE_UNKNOWN || dyType == TYPE_UNKNOWN) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with such image type\n", proc );
    return( -1 );
  }
  
  if ( theIm->ncols != theDx->ncols || theIm->nrows != theDx->nrows
       || theIm->nplanes != theDx->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: input image and X derivative should have the same dimensions\n", proc );
    return( -1 );
  }

  if ( theIm->ncols != theDy->ncols || theIm->nrows != theDy->nrows
       || theIm->nplanes != theDy->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: input image and Y derivative should have the same dimensions\n", proc );
    return( -1 );
  }



  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;


  
  filter[0].type = theFilter;
  filter[0].coefficient = theSigma->x;
  borders[0] = (int)floor( theSigma->x );
  
  filter[1].type = theFilter;
  filter[1].coefficient = theSigma->y;
  borders[1] = (int)floor( theSigma->y );

  filter[2].type = theFilter;
  filter[2].coefficient = theSigma->z;
  borders[2] = (int)floor( theSigma->z );
  filter[2].derivative = NODERIVATIVE;


  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = SMOOTHING;

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)theDx->data, dxType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute X derivative\n", proc );
    return( -1 );
  }

  filter[0].derivative = SMOOTHING;
  filter[1].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)theDy->data, dyType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute Y derivative\n", proc );
    return( -1 );
  }

  return( 1 );
}



int BAL_3DDerivativesOfImage( bal_image *theIm, 
                              bal_image *theDx, bal_image *theDy, bal_image *theDz,
                              bal_doublePoint *theSigma )
{
  char *proc = "BAL_3DDerivativesOfImage";
  bufferType theType = TYPE_UNKNOWN;
  bufferType dxType = TYPE_UNKNOWN;
  bufferType dyType = TYPE_UNKNOWN;
  bufferType dzType = TYPE_UNKNOWN;
  int theDim[3];
  int borders[3] = {0, 0, 0};
  typeFilteringCoefficients filter[3];

  

  initFilteringCoefficients( &(filter[0]) );
  initFilteringCoefficients( &(filter[1]) );
  initFilteringCoefficients( &(filter[2]) );

  theType = theIm->type;
  dxType = theDx->type;
  dyType = theDy->type;
  dzType = theDz->type;
  if ( theType == TYPE_UNKNOWN || dxType == TYPE_UNKNOWN || dyType == TYPE_UNKNOWN || dyType == TYPE_UNKNOWN) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to deal with such image type\n", proc );
    return( -1 );
  }
  
  if ( theIm->ncols != theDx->ncols || theIm->nrows != theDx->nrows
       || theIm->nplanes != theDx->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: input image and X derivative should have the same dimensions\n", proc );
    return( -1 );
  }

  if ( theIm->ncols != theDy->ncols || theIm->nrows != theDy->nrows
       || theIm->nplanes != theDy->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: input image and Y derivative should have the same dimensions\n", proc );
    return( -1 );
  }

  if ( theIm->ncols != theDz->ncols || theIm->nrows != theDz->nrows
       || theIm->nplanes != theDz->nplanes ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: input image and Z derivative should have the same dimensions\n", proc );
    return( -1 );
  }



  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;


  
  filter[0].type = theFilter;
  filter[0].coefficient = theSigma->x;
  borders[0] = (int)floor( theSigma->x );
  
  filter[1].type = theFilter;
  filter[1].coefficient = theSigma->y;
  borders[1] = (int)floor( theSigma->y );

  filter[2].type = theFilter;
  filter[2].coefficient = theSigma->z;
  borders[2] = (int)floor( theSigma->z );



  filter[0].derivative = DERIVATIVE_1;
  filter[1].derivative = SMOOTHING;
  filter[2].derivative = SMOOTHING;

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)theDx->data, dxType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute X derivative\n", proc );
    return( -1 );
  }

  filter[0].derivative = SMOOTHING;
  filter[1].derivative = DERIVATIVE_1;
  filter[2].derivative = SMOOTHING;

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)theDy->data, dyType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute Y derivative\n", proc );
    return( -1 );
  }

  filter[0].derivative = SMOOTHING;
  filter[1].derivative = SMOOTHING;
  filter[2].derivative = DERIVATIVE_1;

  if ( separableLinearFiltering( (void*)theIm->data, theType,
                                 (void*)theDz->data, dzType,
                                 theDim, borders, filter )  != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to compute Z derivative\n", proc );
    return( -1 );
  }

  return( 1 );
}




