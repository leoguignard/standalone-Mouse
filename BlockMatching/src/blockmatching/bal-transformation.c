/*************************************************************************
 * bal-transformation.c -
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
#include <string.h>
#include <sys/stat.h>

#include <vtmalloc.h>

#include <bal-transformation.h>



static int _debug_ = 0;
static int _verbose_ = 1;



int BAL_GetVerboseInBalTransformation( )
{
  return( _verbose_ );
}

void BAL_SetVerboseInBalTransformation( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformation(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformation(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}



static int _fileSize( char *name );










/*--------------------------------------------------
 *
 * TRANSFORMATION MANAGEMENT
 *
 --------------------------------------------------*/

void BAL_InitTransformation( bal_transformation *t ) 
{
  t->type = UNDEF_TRANSFORMATION; 

  /* matrix dimensions are set to 0
     pointer (ie t->mat.m) is set to NULL
  */
  _init_mat( &(t->mat) );
  
  t->transformation_unit = REAL_UNIT;
  
  BAL_InitImage( &(t->vx), NULL, 0, 0, 0, 0, FLOAT );
  BAL_InitImage( &(t->vy), NULL, 0, 0, 0, 0, FLOAT );
  BAL_InitImage( &(t->vz), NULL, 0, 0, 0, 0, FLOAT );

  t->weight = 1.0;
  t->error = 0.0;
  
}




/* Allocate a transformation
   
   Default is REAL_UNIT
*/



static int _InitAllocImageFromImage( bal_image *image, char *name,
                                     bal_image *from, bufferType type )
{
  /* this is mostly BAL_InitAllocImageFromImage() (bal-image.c),
   * however, we have to pay attention that the reference image
   * may be vectorial
   */

  if ( BAL_InitImageFromImage( image, name, from, type ) != 1 ) {
    return( -1 );
  }
  /* in case of vectorial template image
   */
  image->vdim = 1;

  if ( BAL_AllocImage( image ) != 1 ) {
    BAL_FreeImage( image );
    return( -1 );
  }
  return( 1 );

}



int BAL_AllocTransformation( bal_transformation *t, enumTypeTransfo type, bal_image *ref ) 
{
  char *proc = "BAL_AllocTransformation";

  BAL_InitTransformation( t );

  t->type = type;

  if ( ref != NULL ) {
    if ( ref->nplanes == 1 ) {
      switch ( t->type ) {
      default :
        break;
      case TRANSLATION_3D :
        t->type = TRANSLATION_2D;
        break;
      case TRANSLATION_SCALING_3D :
        t->type = TRANSLATION_SCALING_2D;
        break;
      case RIGID_3D :
        t->type = RIGID_2D;
        break;
      case SIMILITUDE_3D :
        t->type = SIMILITUDE_2D;
        break;
      case AFFINE_3D :
        t->type = AFFINE_2D;
        break;
      case VECTORFIELD_3D :
        t->type = VECTORFIELD_2D;
        break;
      }
    }
  }

  if ( _verbose_ >= 5 ) {
      fprintf( stderr, "%s: allocate transformation of type ", proc );
      BAL_PrintTransformationType( stderr, t->type );
  }

  switch ( t->type ) {

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
    
    if ( _alloc_mat(&(t->mat), 4, 4) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate matrix\n", proc );
      return( -1 );
    }
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    if ( ref == (bal_image *)NULL ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: template image was NULL\n", proc );
      return( -1 );
    }

    /* x component
     */
    if ( _InitAllocImageFromImage (&(t->vx), "component_x.inr", ref, FLOAT) == -1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate vector field #1\n", proc );
      return ( -1 );
    }

    /* y component
     */
    if ( _InitAllocImageFromImage (&(t->vy), "component_y.inr", ref, FLOAT) == -1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate vector field #2\n", proc );
      BAL_FreeImage( &(t->vx) );
      return ( -1 );
    }
    t->vy.vx = ref->vx;
    t->vy.vy = ref->vy;
    t->vy.vz = ref->vz;


    if ( t->type == VECTORFIELD_2D ) {

      /* z component
         initialized but not allocated
      */
      
      if ( BAL_InitImageFromImage( &(t->vz), (char*)NULL, ref, FLOAT ) == -1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to initialize vector field #3\n", proc );
        BAL_FreeImage( &(t->vy) );
        BAL_FreeImage( &(t->vx) );
        return ( -1 );
      }
      t->vz.vdim = 1;

    }
    else {

      /* z component
       */
      if ( _InitAllocImageFromImage (&(t->vz), "component_z.inr", ref, FLOAT) == -1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate vector field #3\n", proc );
        BAL_FreeImage( &(t->vx) );
        BAL_FreeImage( &(t->vy) );
        return ( -1 );
      }
      t->vz.vx = ref->vx;
      t->vz.vy = ref->vy;
      t->vz.vz = ref->vz;
      
    }
    break;

  }
  
  BAL_SetTransformationToIdentity( t );

  t->transformation_unit = REAL_UNIT;

  return ( 1 );
}





void BAL_FreeTransformation( bal_transformation *t )
{
  char *proc = "BAL_FreeTransformation";
  
  switch ( t->type ) {

  default :

    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return;

  case UNDEF_TRANSFORMATION :
    break;

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
    
    _free_mat( &(t->mat) );
    break;

  case VECTORFIELD_3D :

    BAL_FreeImage( &(t->vz) );

  case VECTORFIELD_2D :

    BAL_FreeImage( &(t->vx) );
    BAL_FreeImage( &(t->vy) );
    break;

  }
  BAL_InitTransformation( t );

}











void BAL_SetTransformationToIdentity( bal_transformation *t ) 
{
  char *proc = "BAL_SetTransformationToIdentity";
  size_t i, v;
  float *buf;
  
  switch ( t->type ) {

  default :

    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return;
    
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
    
    for ( i=0; i<16; i++ ) t->mat.m[i] = 0.0;
    t->mat.m[0] = t->mat.m[5] = t->mat.m[10] = t->mat.m[15] = 1.0;
    break;

  case VECTORFIELD_3D :

    v = t->vz.ncols * t->vz.nrows * t->vz.nplanes;
    buf = (float*)(t->vz.data);
    for ( i=0; i<v; i++ ) buf[i] = 0.0;
    
  case VECTORFIELD_2D :

    v = t->vx.ncols * t->vx.nrows * t->vx.nplanes;
    buf = (float*)(t->vx.data);
    for ( i=0; i<v; i++ ) buf[i] = 0.0;

    v = t->vy.ncols * t->vy.nrows * t->vy.nplanes;
    buf = (float*)(t->vy.data);
    for ( i=0; i<v; i++ ) buf[i] = 0.0;

    break;

  }
  
}










/*--------------------------------------------------
 *
 * TRANSFORMATION LIST MANAGEMENT
 *
 --------------------------------------------------*/

void BAL_InitTransformationList( bal_transformationList *l )
{
  l->data = (bal_transformation *)NULL;
  l->pointer = (bal_transformation **)NULL;
  l->n_selected_trsfs = 0;
  l->n_trsfs = 0;
  l->n_allocated_trsfs = 0;
}



int BAL_AllocTransformationList( bal_transformationList *l, 
                                 int n_transformations )
{
  char *proc = "BAL_AllocTransformationList";
  int i;

  if ( n_transformations <= 0 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: negative number of transformations (%d)\n", proc, n_transformations );
    return( -1 );
  }
  l->data = (bal_transformation *)vtmalloc( n_transformations*sizeof(bal_transformation),
                                            "l->data", proc );
  if ( l->data == (bal_transformation *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate transformation list (data)\n", proc );
    return( -1 );
  }
  
  l->pointer = (bal_transformation **)vtmalloc( n_transformations*sizeof(bal_transformation*),
                                                "l->pointer", proc );
  if ( l->pointer == (bal_transformation **)NULL ) {
    vtfree( l->data );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate transformation list (pointer)\n", proc );
    return( -1 );
  }

  for ( i=0; i<n_transformations; i++ ) {
    BAL_InitTransformation( &(l->data[i]) );
    l->pointer[i] = &(l->data[i]);
  }

  l->n_allocated_trsfs = n_transformations;

  return( 1 );
}





int BAL_FullAllocTransformationList( bal_transformationList *l, 
                                     int n_transformations,
                                     enumTypeTransfo type, bal_image *ref )
{
  char *proc = "BAL_FullAllocTransformationList";
  int i, j;

  if ( BAL_AllocTransformationList( l, n_transformations ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate transformation list\n", proc );
    return( -1 );
  }

  for ( i=0; i<n_transformations; i++ ) {
    
    if ( BAL_AllocTransformation( &(l->data[i]), type, ref ) != 1 ) {
      for ( j=0; j<i; j++ )
        BAL_FreeTransformation( &(l->data[j]) );
      if ( _verbose_ ) 
        fprintf( stderr, "%s: allocation fails for transformation #%d\n", proc, i );
      return( -1 );
    }
  }

  l->n_trsfs = n_transformations;

  return( 1 );
}





void BAL_FreeTransformationList( bal_transformationList *l )
{
  int i;
  
  if ( l->data != (bal_transformation *)NULL ) {
    for ( i=0; i<l->n_allocated_trsfs; i++ ) {
      BAL_FreeTransformation( &(l->data[i]) );
    }
    vtfree( l->data );
  }

  if ( l->pointer != (bal_transformation **)NULL ) {
    vtfree( l->pointer );
  }

  BAL_InitTransformationList( l );
}













/*--------------------------------------------------
 *
 * TRANSFORMATION ARRAY MANAGEMENT
 *
 --------------------------------------------------*/

void BAL_InitTransformationArray( bal_transformationArray *l )
{
  l->data = (bal_transformation *)NULL;
  l->array = (bal_transformation **)NULL;
  l->n_allocated_trsfs = 0;
  l->ncols = 0;
  l->nrows = 0;
}



int BAL_AllocTransformationArray( bal_transformationArray *l, 
                                  int ncols, int nrows )
{
  char *proc = "BAL_AllocTransformationArray";
  int i, j;
  int n_transformations;
  
  if ( ncols <= 0 || nrows <= 0 ) {
    fprintf( stderr, "%s: negative dimensions\n", proc );
    return( -1 );
  }

  /* allocate data
   */
  n_transformations = ncols*nrows;
  l->data = (bal_transformation *)vtmalloc( n_transformations*sizeof(bal_transformation),
                                            "l->data", proc );
  if ( l->data == (bal_transformation *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate transformation list\n", proc );
    return( -1 );
  }

  for ( i=0; i<n_transformations; i++ )
    BAL_InitTransformation( &(l->data[i]) );

  /* allocate array
   */
  l->array = (bal_transformation **)vtmalloc( nrows * sizeof( bal_transformation* ),
                                              "l->array", proc );
  if ( l->array == (bal_transformation **)NULL ) {
    vtfree( l->data );
    fprintf( stderr, "%s: unable to allocate array\n", proc );
    return( -1 );
  }
  
  /* build array
   */
  for ( j=0, i=0 ; j<nrows; j++, i+=ncols )
    l->array[j] = &(l->data[i]);
  l->nrows = nrows;
  l->ncols = ncols;

  l->n_allocated_trsfs = n_transformations;
  
  return( 1 );
}



void BAL_FreeTransformationArray( bal_transformationArray *l )
{
  int i;
  
  if ( l->data != (bal_transformation *)NULL ) {
    for ( i=0; i<l->n_allocated_trsfs; i++ )
      BAL_FreeTransformation( &(l->data[i]) );
    vtfree( l->data );
  }
  if ( l->array != (bal_transformation **)NULL ) {
    vtfree( l->array );
  }
  BAL_InitTransformationArray( l );
}







/*--------------------------------------------------
 *
 * TRANSFORMATION I/O
 *
 --------------------------------------------------*/

static void _PrintVectorFieldStatistics( FILE *f, bal_image *vx, bal_image *vy, bal_image *vz )
{
  size_t i;
  size_t v = vx->ncols * vx->nrows * vx->nplanes;
  double n, nmin, nmax;

  fprintf( f, "Vector field statistics:\n" );

  switch ( vx->type ) {
  default :                
    fprintf( f, "TYPE_UNKNOWN\n" );
    break;

  case FLOAT :
    {
      float *xbuf = (float *)vx->data;
      float *ybuf = (float *)vy->data;
      float *zbuf;
      float vmin[3] = {0,0,0};
      float vmax[3] = {0,0,0};
      float vm[3] = {0,0,0};

      if ( vz != (bal_image*)NULL ) {

        zbuf = (float *)vz->data;

        vm[0] = vmin[0] = vmax[0] = xbuf[0];
        vm[1] = vmin[1] = vmax[1] = ybuf[0];
        vm[2] = vmin[2] = vmax[2] = zbuf[0];

        nmin = nmax = (xbuf[0] * xbuf[0] + ybuf[0] * ybuf[0] + zbuf[0] * zbuf[0]);

        for ( i=1; i<v; i++ ) {
          vm[0] += xbuf[i];
          vm[1] += ybuf[i];
          vm[2] += zbuf[i];
          n = (xbuf[i] * xbuf[i] + ybuf[i] * ybuf[i] + zbuf[i] * zbuf[i]);
          if ( nmin > n ) {
            vmin[0] = xbuf[i];
            vmin[1] = ybuf[i];
            vmin[2] = zbuf[i];
            nmin = n;
          }
          if ( nmax < n ) {
            vmax[0] = xbuf[i];
            vmax[1] = ybuf[i];
            vmax[2] = zbuf[i];
            nmax = n;
          }
        }
        fprintf( f, "  - mean vector = " );
        fprintf( f, "[%f,%f,%f]\n", vm[0]/v, vm[1]/v, vm[2]/v );
        fprintf( f, "  - modulus:    " );
        fprintf( f, "min = [%f,%f,%f], max = [%f,%f,%f]\n",
                 vmin[0], vmin[1], vmin[2],
                 vmax[0], vmax[1], vmax[2] );
      }
      else {

        vm[0] = vmin[0] = vmax[0] = xbuf[0];
        vm[1] = vmin[1] = vmax[1] = ybuf[0];

        nmin = nmax = (xbuf[0] * xbuf[0] + ybuf[0] * ybuf[0]);

        for ( i=1; i<v; i++ ) {
          vm[0] += xbuf[i];
          vm[1] += ybuf[i];
          n = (xbuf[i] * xbuf[i] + ybuf[i] * ybuf[i]);
          if ( nmin > n ) {
            vmin[0] = xbuf[i];
            vmin[1] = ybuf[i];
            nmin = n;
          }
          if ( nmax < n ) {
            vmax[0] = xbuf[i];
            vmax[1] = ybuf[i];
            nmax = n;
          }
        }

        fprintf( f, "  - mean vector = " );
        fprintf( f, "[%f,%f]\n", vm[0]/v, vm[1]/v );
        fprintf( f, "  - modulus:    " );
        fprintf( f, "min = [%f,%f], max = [%f,%f]\n",
                 vmin[0], vmin[1],
                 vmax[0], vmax[1] );
      }
      
    }
    break;
  }  
}





void BAL_PrintTransformationType( FILE *f, enumTypeTransfo t )
{
  char *proc = "BAL_PrintTransformationType"; 

  switch ( t ) {
  default : fprintf( f, "not handled yet in %s\n", proc );  break;
  case UNDEF_TRANSFORMATION :   fprintf( f, "UNDEF_TRANSFORMATION\n" ); break;
  case TRANSLATION_2D :         fprintf( f, "TRANSLATION_2D\n" ); break;
  case TRANSLATION_3D :         fprintf( f, "TRANSLATION_3D\n" ); break;
  case TRANSLATION_SCALING_2D : fprintf( f, "TRANSLATION_SCALING_2D\n" ); break;
  case TRANSLATION_SCALING_3D : fprintf( f, "TRANSLATION_SCALING_3D\n" ); break;
  case RIGID_2D :               fprintf( f, "RIGID_2D\n" ); break;
  case RIGID_3D :               fprintf( f, "RIGID_3D\n" ); break;
  case SIMILITUDE_2D :          fprintf( f, "SIMILITUDE_2D\n" ); break;
  case SIMILITUDE_3D :          fprintf( f, "SIMILITUDE_3D\n" ); break;
  case AFFINE_2D :              fprintf( f, "AFFINE_2D\n" ); break;
  case AFFINE_3D :              fprintf( f, "AFFINE_3D\n" ); break;
  case VECTORFIELD_2D :         fprintf( f, "VECTORFIELD_2D\n" ); break;
  case VECTORFIELD_3D :         fprintf( f, "VECTORFIELD_3D\n" ); break;
  }
}



void BAL_PrintTransformationUnit( FILE *f, enumUnitTransfo t )
{
  char *proc = "BAL_PrintTransformationUnit";

  switch( t ) {
  default : fprintf( f, "not handled yet in %s\n", proc );  break;
  case UNDEF_UNIT : fprintf( f, "UNDEF_UNIT\n" ); break;
  case VOXEL_UNIT : fprintf( f, "VOXEL_UNIT\n" ); break;
  case REAL_UNIT :  fprintf( f, "REAL_UNIT\n" ); break;
  }
}



void BAL_PrintTransformation( FILE *f, bal_transformation *t, char *s )
{
  char *proc = "BAL_PrintTransformation"; 
  float det, scale = 1.0;
  float trace = 10.0, theta = 0.0;

  if ( t == (bal_transformation *)NULL ) {
    if ( s != (char *)NULL )
      fprintf( f, "'%s' is NULL\n", s );
    else
      fprintf( f, "transformation is NULL\n" );
    return;
  }

  if ( s != (char *)NULL )
    fprintf( f, "type of '%s' is ", s );
  else
    fprintf( f, "transformation type is " );
  BAL_PrintTransformationType( f, t->type );


  if ( s != (char *)NULL )
    fprintf( f, "unit of '%s' is ", s );
  else
    fprintf( f, "transformation unit is " );
  BAL_PrintTransformationUnit( f, t->transformation_unit );

  
  switch ( t->type ) {
  default :
    fprintf( f, "%s: such transformation type not handled yet\n", proc );
    return;
    
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
    /*  0  1  2  3
     *  4  5  6  7
     *  8  9 10 11
     * 12 13 14 15
     */
    switch( t->type ) {
    default :
        break;
    case SIMILITUDE_2D :
        det = t->mat.m[ 0] * t->mat.m[ 5]
                - t->mat.m[ 1] * t->mat.m[ 4];
        break;
    case SIMILITUDE_3D :
        det = t->mat.m[ 0] * t->mat.m[ 5] * t->mat.m[10]
                + t->mat.m[ 4] * t->mat.m[ 9] * t->mat.m[ 2]
                + t->mat.m[ 8] * t->mat.m[ 1] * t->mat.m[ 6]
                - t->mat.m[ 0] * t->mat.m[ 6] * t->mat.m[ 9]
                - t->mat.m[ 1] * t->mat.m[ 4] * t->mat.m[10]
                - t->mat.m[ 2] * t->mat.m[ 5] * t->mat.m[ 8];
        break;
    }
    scale = sqrt(det);

  case RIGID_2D :
  case RIGID_3D :

    switch( t->type ) {
    default :
        break;
    case RIGID_2D :
        scale = 1.0;
    case SIMILITUDE_2D :
        theta = atan2( t->mat.m[ 4]/scale, t->mat.m[ 0]/scale );
        break;
    case RIGID_3D :
        scale = 1.0;
    case SIMILITUDE_3D :
        trace = t->mat.m[ 0] + t->mat.m[ 5] + t->mat.m[10];
        trace /= scale;
        trace -= 1.0;
        trace /= 2.0;
        theta = acos( trace );
        break;
    }

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    _print_mat ( f, &(t->mat), (char *)NULL );
    if ( t->type == SIMILITUDE_2D || t->type == SIMILITUDE_3D )
      fprintf( f, "   scale is %f\n", scale );
    if ( t->type == SIMILITUDE_2D || t->type == SIMILITUDE_3D
         || t->type == RIGID_2D || t->type == RIGID_3D ) {
        fprintf( f, "   angle is %f (%f degrees)\n", theta, theta * 180.0 / 3.14159265359 );
    }
    break;

  case VECTORFIELD_2D :

    BAL_PrintImage( f, &(t->vx), "X component of vector field" );
    BAL_PrintImage( f, &(t->vy), "Y component of vector field" );
    _PrintVectorFieldStatistics( f, &(t->vx), &(t->vy), (bal_image *)NULL );
    break;

  case VECTORFIELD_3D :

    BAL_PrintImage( f, &(t->vx), "X component of vector field" );
    BAL_PrintImage( f, &(t->vy), "Y component of vector field" );
    BAL_PrintImage( f, &(t->vz), "Z component of vector field" );
    _PrintVectorFieldStatistics( f, &(t->vx), &(t->vy), &(t->vz) );
    break;

  }
}







/* Deformation norm calculation
 */
int BAL_TransformationModulus( bal_transformation *theTrsf, bal_image *image )
{
  char *proc = "BAL_TransformationModulus";
  size_t i, j, k;
  float u, v, w;
  float x, y, z;
  double *mat = (double*)NULL;
  double *to_r = image->to_real.m;
  float ***vx = (float ***)NULL;
  float ***vy = (float ***)NULL;
  float ***vz = (float ***)NULL;


  switch( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case FLOAT :
    {
      float ***res= (float***)image->array;
      switch( theTrsf->type ) {
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
        mat = theTrsf->mat.m;
        switch( theTrsf->transformation_unit ) {
        default :
          if ( _verbose_ )
            fprintf( stderr, "%s: such transformation unit not handled yet\n", proc );
          return( -1 );
        case VOXEL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            u = mat[ 0] * i + mat[ 1] * j + mat[ 2] * k + mat[ 3] - i;
            v = mat[ 4] * i + mat[ 5] * j + mat[ 6] * k + mat[ 7] - j;
            w = mat[ 8] * i + mat[ 8] * j + mat[10] * k + mat[11] - k;
            res[k][j][i] = sqrt( u*u + v*v + w*w );
          }
          break;
        case REAL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            x = to_r[ 0] * i + to_r[ 1] * j + to_r[ 2] * k + to_r[ 3];
            y = to_r[ 4] * i + to_r[ 5] * j + to_r[ 6] * k + to_r[ 7];
            z = to_r[ 8] * i + to_r[ 8] * j + to_r[10] * k + to_r[11];
            u = mat[ 0] * x + mat[ 1] * y + mat[ 2] * z + mat[ 3] - x;
            v = mat[ 4] * x + mat[ 5] * y + mat[ 6] * z + mat[ 7] - y;
            w = mat[ 8] * x + mat[ 8] * y + mat[10] * z + mat[11] - z;
            res[k][j][i] = sqrt( u*u + v*v + w*w );
          }
          break;
        }
        break;
      case VECTORFIELD_2D :
        vx = (float***)(theTrsf->vx).array;
        vy = (float***)(theTrsf->vy).array;
        switch( theTrsf->transformation_unit ) {
        default :
          if ( _verbose_ )
            fprintf( stderr, "%s: such transformation unit not handled yet\n", proc );
          return( -1 );
        case VOXEL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            u = vx[k][j][i];
            v = vy[k][j][i];
            res[k][j][i] = sqrt( u*u + v*v );
          }
          break;
        case REAL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            x = to_r[ 0] * i + to_r[ 1] * j + to_r[ 2] * k + to_r[ 3];
            y = to_r[ 4] * i + to_r[ 5] * j + to_r[ 6] * k + to_r[ 7];
            u = vx[k][j][i] - x;
            v = vy[k][j][i] - y;
            res[k][j][i] = sqrt( u*u + v*v );
          }
          break;
        }
        break; /* theTrsf->type == VECTORFIELD_2D */
      case VECTORFIELD_3D :
        vx = (float***)(theTrsf->vx).array;
        vy = (float***)(theTrsf->vy).array;
        vz = (float***)(theTrsf->vz).array;
        switch( theTrsf->transformation_unit ) {
        default :
          if ( _verbose_ )
            fprintf( stderr, "%s: such transformation unit not handled yet\n", proc );
          return( -1 );
        case VOXEL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            u = vx[k][j][i];
            v = vy[k][j][i];
            w = vz[k][j][i];
            res[k][j][i] = sqrt( u*u + v*v + w*w );
          }
          break;
        case REAL_UNIT :
          for ( k=0; k<image->nplanes; k++ )
          for ( j=0; j<image->nrows; j++ )
          for ( i=0; i<image->ncols; i++ ) {
            x = to_r[ 0] * i + to_r[ 1] * j + to_r[ 2] * k + to_r[ 3];
            y = to_r[ 4] * i + to_r[ 5] * j + to_r[ 6] * k + to_r[ 7];
            z = to_r[ 8] * i + to_r[ 8] * j + to_r[10] * k + to_r[11];
            u = vx[k][j][i] - x;
            v = vy[k][j][i] - y;
            w = vz[k][j][i] - z;
            res[k][j][i] = sqrt( u*u + v*v + w*w );
          }
          break;
        }
        break; /* theTrsf->type == VECTORFIELD_3D */
      } /* switch( theTrsf->type ) */
    }
    break; /* image->type == FLOAT */
  }

  return( 1 );
}







/* Read a transformation
   
   Transformations are supposed to be in REAL_UNIT
*/

int BAL_ReadTransformation( bal_transformation *theTrsf, char *name )
{
  char *proc = "BAL_ReadTransformation";
  bal_image theIm;
  int size;
  int distant_verbose;
  size_t i, j ,k;

  if ( name == (char*)NULL || name[0] == '\0' ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: invalid file name '%s'\n", proc, name );
      return( -1 );
  }

  size = _fileSize( name );
  if ( size == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: an error occurs when accessing to '%s'\n", proc, name );
    return( -1 );
  }
  
  /* on a bench of matrix files, the size was between 153 and 163
     if less than 200, should be a matrix
   */
  
  if ( size <= 200 ) {

    if ( BAL_AllocTransformation( theTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to allocate transformation (matrix)\n", proc );
      return( -1 );
    }

    distant_verbose = BAL_GetVerboseInBalMatrix();
    BAL_SetVerboseInBalMatrix( 0 );

    if ( _read_mat( name, &(theTrsf->mat) ) == 1 ) {
      BAL_SetVerboseInBalMatrix( distant_verbose );
      return( 1 );
    }    

    BAL_SetVerboseInBalMatrix( distant_verbose );
    BAL_FreeTransformation( theTrsf );
    /* this is not a matrix
     */
  }

  

  /* unset errors when reading an image
     (in case it's a matrix)
     unless it's a debugging mode
   */
  distant_verbose = BAL_GetVerboseInBalImage();
  if ( _debug_  ) {
    if ( distant_verbose >= 1 ) BAL_SetVerboseInBalImage( distant_verbose );
    else BAL_SetVerboseInBalImage( 1 );
  }
  else {
    BAL_SetVerboseInBalImage( 0 );
  }



  /* try to read a vector field
   */

  if ( BAL_ReadImage( &theIm, name, 0 ) == 1 ) {
    
    if ( _debug_ ) {
      fprintf( stderr, "\n" );
      fprintf( stderr, "%s: reading image '%s'", proc, name );
      BAL_PrintImage( stderr, &theIm, "read image" );
      fprintf( stderr, "\n" );
    }

    BAL_SetVerboseInBalImage( distant_verbose );

    switch( theIm.vdim ) {
    default : 

      if ( _verbose_ ) 
        fprintf( stderr, "%s: '%s' is not a vectorial image (with dimv = 2 or 3)\n", proc, name );
      BAL_FreeImage( &theIm );
      return( -1 );

    case 2 :

      if ( BAL_AllocTransformation( theTrsf, VECTORFIELD_2D, &theIm ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate transformation (2D vector field)\n", proc );
        BAL_FreeImage( &theIm );
        return( -1 );
      }

      switch ( theIm.type ) {
      default : 
        if ( _verbose_ )
          fprintf( stderr, "%s: image type of '%s' not handled yet (2D vector field)\n", proc, name );
        BAL_FreeTransformation( theTrsf );
        BAL_FreeImage( &theIm );
        return( -1 );

      case FLOAT :
        {
          float ***vx = (float***)(theTrsf->vx.array);
          float ***vy = (float***)(theTrsf->vy.array);
          float ***vec = (float***)(theIm.array);

          for ( k=0; k<theIm.nplanes; k++ )
          for ( j=0; j<theIm.nrows; j++ )
          for ( i=0; i<theIm.ncols; i++ ) {
            vx[k][j][i] = vec[k][j][2*i];
            vy[k][j][i] = vec[k][j][2*i+1];
          }
        }
        break;
      }
      break;
      /* end of 2D vector field
       */
      
    case 3 :

      if ( BAL_AllocTransformation( theTrsf, VECTORFIELD_3D, &theIm ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate transformation (3D vector field)\n", proc );
        BAL_FreeImage( &theIm );
        return( -1 );
      }

      switch ( theIm.type ) {
      default : 
        if ( _verbose_ )
          fprintf( stderr, "%s: image type of '%s' not handled yet (3D vector field)\n", proc, name );
        BAL_FreeTransformation( theTrsf );
        BAL_FreeImage( &theIm );
        return( -1 );

      case FLOAT :
        {
          float ***vx = (float***)(theTrsf->vx.array);
          float ***vy = (float***)(theTrsf->vy.array);
          float ***vz = (float***)(theTrsf->vz.array);
          float ***vec = (float***)(theIm.array);

          for ( k=0; k<theIm.nplanes; k++ )
          for ( j=0; j<theIm.nrows; j++ )
          for ( i=0; i<theIm.ncols; i++ ) {
            vx[k][j][i] = vec[k][j][3*i];
            vy[k][j][i] = vec[k][j][3*i+1];
            vz[k][j][i] = vec[k][j][3*i+2];
          }
        }
        break;
      }
      break;
      /* end of 3D vector field
       */
      
    }
    BAL_FreeImage( &theIm );
    return( 1 );
  }

  /* this is not a vector field
   */

  else {

    /* try to read a matrix 
     */

    if ( BAL_AllocTransformation( theTrsf, AFFINE_3D, (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to allocate transformation (matrix)\n", proc );
      return( -1 );
    }

    distant_verbose = BAL_GetVerboseInBalMatrix();
    BAL_SetVerboseInBalMatrix( 0 );

    if ( _read_mat( name, &(theTrsf->mat) ) == 1 ) {
      BAL_SetVerboseInBalMatrix( distant_verbose );
      return( 1 );
    }    

    BAL_SetVerboseInBalMatrix( distant_verbose );
    BAL_FreeTransformation( theTrsf );
    /* this is not a matrix
     */
  }


  BAL_SetVerboseInBalImage( distant_verbose );


  if ( _verbose_ ) 
    fprintf( stderr, "%s: transformation type of '%s' was not recognized\n", proc, name );
  
  return( -1 );
}





int BAL_WriteTransformation( bal_transformation *theTrsf, char *name )
{
  char *proc = "BAL_WriteTransformation";
  bal_image theIm;
  size_t i, j, k;

  if ( name == (char*)NULL || name[0] == '\0' ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: invalid file name '%s'\n", proc, name );
      return( -1 );
  }

  switch ( theTrsf->type ) {
  default : 
    if ( _verbose_ )
      fprintf( stderr, "%s: transformation type not handled yet\n", proc );  
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
    
    if ( _write_mat( name, &(theTrsf->mat) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write matrix in '%s'\n", proc, name );
      return( -1 );
    }
    return( 1 );

  case VECTORFIELD_2D :

    if ( BAL_InitImageFromImage( &theIm, (char*)NULL, &(theTrsf->vx), theTrsf->vx.type ) == -1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to initialize 2D vector field\n", proc );
        return( -1 );
    }
    theIm.vdim = 2;
    if ( BAL_AllocImage( &theIm ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate 2D vector field\n", proc );
      return( -1 );
    }
    
    theIm.vx = theTrsf->vx.vx;
    theIm.vy = theTrsf->vx.vy;
    theIm.vz = theTrsf->vx.vz;

    switch ( theTrsf->vx.type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled yet (2D vector field)\n", proc );
      BAL_FreeImage( &theIm );
      return( -1 );
    case FLOAT :
      {
        float ***vx = (float***)(theTrsf->vx.array);
        float ***vy = (float***)(theTrsf->vy.array);
        float ***vec = (float***)(theIm.array);
        for ( k=0; k<theIm.nplanes; k++ )
        for ( j=0; j<theIm.nrows; j++ )
        for ( i=0; i<theIm.ncols; i++ ) {
          vec[k][j][2*i]   = vx[k][j][i];
          vec[k][j][2*i+1] = vy[k][j][i];
        }
      }
      break;
    }
    
    if ( BAL_WriteImage( &theIm, name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write 2D vector field in '%s'\n", proc, name );
      BAL_FreeImage( &theIm );
      return( -1 );
    }
    
    BAL_FreeImage( &theIm );
    return( 1 );
    
  case VECTORFIELD_3D :
    
    if ( BAL_InitImageFromImage( &theIm, (char*)NULL, &(theTrsf->vx), theTrsf->vx.type ) == -1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to initialize 3D vector field\n", proc );
        return( -1 );
    }
    theIm.vdim = 3;
    if ( BAL_AllocImage( &theIm ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate 3D vector field\n", proc );
      return( -1 );
    }

    theIm.vx = theTrsf->vx.vx;
    theIm.vy = theTrsf->vx.vy;
    theIm.vz = theTrsf->vx.vz;

    switch ( theTrsf->vx.type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled yet (3D vector field)\n", proc );
      BAL_FreeImage( &theIm );
      return( -1 );
    case FLOAT :
      {
        float ***vx = (float***)(theTrsf->vx.array);
        float ***vy = (float***)(theTrsf->vy.array);
        float ***vz = (float***)(theTrsf->vz.array);
        float ***vec = (float***)(theIm.array);
        for ( k=0; k<theIm.nplanes; k++ )
        for ( j=0; j<theIm.nrows; j++ )
        for ( i=0; i<theIm.ncols; i++ ) {
          vec[k][j][3*i]   = vx[k][j][i];
          vec[k][j][3*i+1] = vy[k][j][i];
          vec[k][j][3*i+2] = vz[k][j][i];
        }
      }
      break;
    }
    
    if ( BAL_WriteImage( &theIm, name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write 3D vector field in '%s'\n", proc, name );
      BAL_FreeImage( &theIm );
      return( -1 );
    }
    
    BAL_FreeImage( &theIm );
    return( 1 );

  }

  return( 1 );
}









/*--------------------------------------------------
 *
 * TRANSFORMATION LIST I/O
 *
 --------------------------------------------------*/

void BAL_PrintTransformationList( FILE *f, bal_transformationList *l, char *s )
{
  if ( l == (bal_transformationList *)NULL ) {
    if ( s != (char *)NULL )
      fprintf( f, "'%s' is NULL\n", s );
    else
      fprintf( f, "transformation list is NULL\n" );
    return;
  }

  if ( s != (char *)NULL )
    fprintf( f, "information on '%s':", s );
  else
    fprintf( f, "information:" );

  fprintf( f, "   n_selected_trsfs  = %4d\n", l->n_selected_trsfs );
  fprintf( f, "   n_trsfs           = %4d\n", l->n_trsfs );
  fprintf( f, "   n_allocated_trsfs = %4d\n", l->n_allocated_trsfs );
}



int BAL_ReadTransformationList( bal_transformationList *theList, 
                                stringList *trsfFileList )
{
  char *proc = "BAL_ReadTransformationList";
  int i;
  int n_trsfs = theList->n_trsfs;
  int _local_verbose_= _verbose_;

  _verbose_ = 0;

  if ( 0 ) {
    fprintf( stderr, "%s: trsfFileList->n_data = %d\n", proc, trsfFileList->n_data );
    fprintf( stderr, "\t theList->n_trsfs = %d \n", theList->n_trsfs );
    fprintf( stderr, "\t theList->n_allocated_trsfs = %d\n", theList->n_allocated_trsfs );
  }

  for ( i=0; i<trsfFileList->n_data && theList->n_trsfs<theList->n_allocated_trsfs; i++ ) {
    if ( 0 ) {
      fprintf(stderr, "%s: read transformation #%d '%s' -> [%d]\n", 
              proc, i, trsfFileList->data[i], theList->n_trsfs );
    }
    if ( BAL_ReadTransformation( &(theList->data[theList->n_trsfs]), 
                                 trsfFileList->data[i] ) != 1 ) {
      if ( _local_verbose_ )
        fprintf( stderr, "%s: unable to read transformation #%3d '%s'\n",
                 proc, i, trsfFileList->data[i] );
    }
    else {
      theList->n_trsfs ++;
    }
  }

  _verbose_ = _local_verbose_;

  if ( i < trsfFileList->n_data ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read last transformations [%d-%d]\n", 
               proc, i, trsfFileList->n_data-1 );
  }

  if ( n_trsfs == theList->n_trsfs ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no transformations read\n", proc );
    return( -1 );
  }

  theList->n_selected_trsfs = theList->n_trsfs;

  return( 1 );
}





int BAL_WriteTransformationList( bal_transformationList *theList, 
                                 stringList *trsfFileList )
{
  char *proc = "BAL_WriteTransformationList";
  int i;

  if ( theList->n_trsfs > trsfFileList->n_data ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: there are more transformations (%d) than file names (%d)\n",
               proc, theList->n_trsfs, trsfFileList->n_data );
  }
  if ( theList->n_trsfs < trsfFileList->n_data ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: there are less transformations (%d) than file names (%d)\n",
               proc, theList->n_trsfs, trsfFileList->n_data );
  }

  

  for ( i=0; i<theList->n_trsfs && i<trsfFileList->n_data; i++ ) {
    if ( 0 ) {
      fprintf(stderr, "%s: write transformation #%d '%s'\n", 
              proc, i, trsfFileList->data[i] );
    }
    if ( BAL_WriteTransformation( &(theList->data[i]), 
                                 trsfFileList->data[i] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write transformation #%3d '%s'\n",
                 proc, i, trsfFileList->data[i] );
    }
  }

  return( 1 );
}





/*--------------------------------------------------
 *
 * TRANSFORMATION ARRAY I/O
 *
 --------------------------------------------------*/

int BAL_ReadTransformationArray( bal_transformationArray *theArray, 
                                stringArray *trsfFileArray )
{
  char *proc = "BAL_ReadTransformationArray";
  int i, j, n;
  int _local_verbose_ = _verbose_;

  _verbose_ = 0;

  if ( theArray->ncols != trsfFileArray->ncols || theArray->nrows != trsfFileArray->nrows ) {
    _verbose_ = _local_verbose_;
    if ( _verbose_ ) {
      fprintf( stderr, "%s: transformation array and string array have different dimensions\n",
               proc );
    }
    return( -1 );
  }

 for ( j=0; j<theArray->nrows; j++ ) {
   for ( i=0; i<theArray->ncols; i++ ) {
     if ( 0 ) {
       fprintf(stderr, "%s: read transformation [%d][%d] '%s'\n", 
               proc, j, i, trsfFileArray->array[j][i] );
     }
     if ( BAL_ReadTransformation( &(theArray->array[j][i]), 
                                  trsfFileArray->array[j][i] ) != 1 ) {
       if ( _local_verbose_ >= 3 )
         fprintf( stderr, "%s: unable to read transformation [%3d][%3d] '%s'\n",
                  proc, j, i, trsfFileArray->array[j][i] );
     }
   }
 }

  _verbose_ = _local_verbose_;

 for ( n=0, j=0; j<theArray->nrows; j++ ) {
   for ( i=0; i<theArray->ncols; i++ ) {
     if ( BAL_DoesTransformationExist( &(theArray->array[j][i]) ) == 1 )
       n ++;
   }
 }

 if ( _verbose_ >=2 ) {
   fprintf( stderr, "%s: read %d/[%d=%dx%d] transformations\n", proc,
            n, theArray->nrows*theArray->ncols, theArray->ncols, theArray->nrows );
 }

 if ( n == 0 ) {
   if ( _verbose_ )
     fprintf( stderr, "%s: no transformation read, error in names?\n", proc );
   return( -1 );
 }

  return( 1 );
}







/*--------------------------------------------------
 *
 * TRANSFORMATION TESTS
 *
 --------------------------------------------------*/


int BAL_DoesTransformationExist( bal_transformation *t ) 
{
  char *proc = "BAL_DoesTransformationExist";

  if ( t == (bal_transformation*)NULL )
    return( 0 );

  switch ( t->type ) {

  default :
  case UNDEF_TRANSFORMATION :
    return( 0 );
    
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
  
    if ( t->mat.l <= 0 || t->mat.c <= 0 || t->mat.m == (double*)NULL )
      return( 0 );
    return( 1 );

  case VECTORFIELD_3D :
    if ( t->vz.ncols <= 0 
         || t->vz.nrows <= 0
         || t->vz.nplanes <= 0
         || t->vz.data == (void*)NULL
         || t->vz.array == (void***)NULL )
      return( 0 );
    
  case VECTORFIELD_2D :
    if ( t->vy.ncols <= 0 
         || t->vy.nrows <= 0
         || t->vy.nplanes <= 0
         || t->vy.data == (void*)NULL
         || t->vy.array == (void***)NULL )
      return( 0 );
    if ( t->vx.ncols <= 0 
         || t->vx.nrows <= 0
         || t->vx.nplanes <= 0
         || t->vx.data == (void*)NULL
         || t->vx.array == (void***)NULL )
      return( 0 );
    return( 1 );

  case SPLINE :
    return( 0 );
  }

  if ( _verbose_ )
    fprintf( stderr, "%s: weird, this should not be reachable\n", proc );

  return( 1 );
}



int BAL_IsTransformationLinear( bal_transformation *t ) 
{
    if ( t == (bal_transformation *)NULL )
        return( -1 );
    return( BAL_IsTransformationTypeLinear( t->type) );
}

int BAL_IsTransformationVectorField( bal_transformation *t ) 
{
    if ( t == (bal_transformation *)NULL )
        return( -1 );
    return( BAL_IsTransformationTypeVectorField( t->type) );
}

int BAL_IsTransformationTypeLinear( enumTypeTransfo type ) 
{
  switch ( type ) {
  default :
    break;
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
    return( 1 );
  }
  return( 0 );
}


int BAL_IsTransformationTypeVectorField( enumTypeTransfo type ) 
{
  switch ( type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    return( 1 );
  }
  return( 0 );
}










/*--------------------------------------------------
 *
 * TRANSFORMATION ARRAY TESTS
 *
 --------------------------------------------------*/

void BAL_TestTransformationArray( bal_transformationArray *theArray, 
                                 stringArray *trsfFileArray )
{
  int i, j;
  for ( j=0; j<theArray->nrows; j++ ) {
    for ( i=0; i<theArray->ncols; i++ ) {
      if ( BAL_DoesTransformationExist( &(theArray->array[j][i]) ) == 1 ) {
        fprintf( stdout, "transformation(%3d,%3d) exists", i, j );
        if ( trsfFileArray != (stringArray*)NULL )
          fprintf( stdout, ", read from '%s'", trsfFileArray->array[j][i] );
        fprintf( stdout, "\n" );
      }
    }
  }

}










/************************************************************
 *
 * MISC
 *
 ************************************************************/

static int _fileSize( char *name )
{
  char *proc = "_fileSize";
  struct stat buf;
  if ( stat( name, &buf ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: can not access to '%s'\n", proc, name );
    return( -1 );
  }
  return( buf.st_size );
}





int BAL_RotationAngle( bal_transformation *theTrsf, double *radian_angle )
{
  char *proc = "BAL_RotationAngle";
  double s;

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
    /* the rotation matrix trace is 2 cos(theta) + 1 
     */
    s = theTrsf->mat.m[ 0] + theTrsf->mat.m[ 5] + theTrsf->mat.m[10];
    s -= 1.0;
    s /= 2.0;

    if ( s < -1.0 || s > 1.0 ) {
      fprintf( stderr, "%s: cosinus value is out of range\n", proc );
      if ( s < -1.0 ) s = -0.99;
      else if ( s > 1.0 ) s = 0.99;
    }
    *radian_angle = acos( s );
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    fprintf( stderr, "%s: vector field transformation type not handled yet\n", proc );
    return( -1 );
    break;

  }

  return( 1 );

}







/*--------------------------------------------------
 *
 * 
 *
 --------------------------------------------------*/










