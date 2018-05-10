/*************************************************************************
 * pointmatching.c - template for executable creation
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 14 sep 2016 18:28:04 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <vtmalloc.h>

#include <bal-point.h>
#include <bal-transformation-copy.h>

#include <api-pointmatching.h>







static int _verbose_ = 1;





/* static function definitions
 */

static char *_Array2Str( int argc, char *argv[] );
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamPointmatching par;
  bal_image imread, imtemplate;
  bal_image *ptrtemplate = (bal_image*)NULL;
  bal_doublePointList theFloatingPoints;
  bal_doublePointList theReferencePoints;
  bal_transformation *resTransformation;
  bal_image imfloat, imreference;
  char *lineoptions;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_pointmatching( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_pointmatching( 1, argc, argv, &par );
  

  /* parameter writing
   */
  if ( par.command_line_file != (char*)NULL && par.command_line_file[0]!= '\0' ) {
    FILE *f;
    f = fopen( par.command_line_file, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: warning, unable to open '%s'\n", _BaseName( argv[0] ), _BaseName( argv[0] ) );
    }
    else {
      int i;
      for ( i=0; i<argc; i++ )
        fprintf( f, "%s ", argv[i] );
      fprintf( f, "\n" );
      fclose( f );
    }
  }


  /* template image
   * for vector field definition
   */
  (void)BAL_InitImage( &imtemplate, NULL, 0, 0, 0, 0, UCHAR );

  switch ( par.transformation_type ) {
  default :
    API_ErrorParse_pointmatching( _BaseName( argv[0] ), "such transformation type non handled yet ...\n", 0 );
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
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    if ( par.template_name != (char*)NULL && par.template_name[0]!= '\0' ) {
      if ( BAL_ReadImage( &imread, par.template_name, 1 ) != 1 ) {
        fprintf( stderr, "%s: unable to read '%s'\n", _BaseName( argv[0] ), par.template_name );
        API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to read template image ...\n", 0 );
      }
      if ( BAL_InitScalarImageFromImage( &imtemplate, NULL, &imread, UCHAR ) != 1 ) {
        BAL_FreeImage( &imread );
        API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
      }
      BAL_FreeImage( &imread );
    }
    else if ( par.dim.x > 0.0&& par.dim.y > 0.0 ) {
      if ( par.dim.z > 0.0 ) {
        if ( BAL_InitImage( &imtemplate, (char*)NULL, par.dim.x, par.dim.y, par.dim.z, 1, UCHAR ) != 1 ) {
          API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
        }
      }
      else {
        if ( BAL_InitImage( &imtemplate, (char*)NULL, par.dim.x, par.dim.y, 1, 1, UCHAR ) != 1 ) {
          API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
        }
      }
      imtemplate.vx = (par.voxel.x > 0) ? par.voxel.x : ((par.reference_voxel.x > 0) ? par.reference_voxel.x : 1.0);
      imtemplate.vy = (par.voxel.y > 0) ? par.voxel.y : ((par.reference_voxel.y > 0) ? par.reference_voxel.y : 1.0);
      imtemplate.vz = (par.voxel.z > 0) ? par.voxel.z : ((par.reference_voxel.z > 0) ? par.reference_voxel.z : 1.0);
      if ( BAL_SetImageVoxelSizes( &imtemplate, imtemplate.vx, imtemplate.vy, imtemplate.vz ) != 1 ) {
        API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to initialize template image voxel sizes\n", 0 );
      }
    }
    else {
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "no information to build template image ...\n", 0 );
    }
    ptrtemplate = &imtemplate;
    break;
  }



  /* reading points
   */
  BAL_InitDoublePointList( &theFloatingPoints );
  BAL_InitDoublePointList( &theReferencePoints );


  if ( par.floating_points != (char*)NULL && par.floating_points[0]!= '\0' ) {
    if ( BAL_ReadDoublePointList( &theFloatingPoints, par.floating_points ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening/reading '%s'\n",
                 _BaseName( argv[0] ),  par.floating_points );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "error when opening/reading floating point list...\n", 0 );
    }
  }
  else {
    API_ErrorParse_pointmatching( _BaseName( argv[0] ), "no floating point list...\n", 0 );
  }

  theFloatingPoints.unit =  par.points_unit;
  if ( par.floating_voxel.x > 0.0 ) {
    theFloatingPoints.vx = par.floating_voxel.x;
  } else if ( par.voxel.x > 0.0 ) {
    theFloatingPoints.vx = par.voxel.x;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theFloatingPoints.vx = ptrtemplate->vx;
  }
  if ( par.floating_voxel.y > 0.0 ) {
    theFloatingPoints.vy = par.floating_voxel.y;
  } else if ( par.voxel.y > 0.0 ) {
    theFloatingPoints.vy = par.voxel.y;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theFloatingPoints.vy = ptrtemplate->vy;
  }
  if ( par.floating_voxel.z > 0.0 ) {
    theFloatingPoints.vz = par.floating_voxel.z;
  } else if ( par.voxel.z > 0.0 ) {
    theFloatingPoints.vz = par.voxel.z;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theFloatingPoints.vz = ptrtemplate->vz;
  }



  if ( par.reference_points != (char*)NULL && par.reference_points[0]!= '\0' ) {
    if ( BAL_ReadDoublePointList( &theReferencePoints, par.reference_points ) != 1 ) {
      BAL_FreeDoublePointList( &theFloatingPoints );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening/reading '%s'\n",
                 _BaseName( argv[0] ),  par.reference_points );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "error when opening/reading reference point list...\n", 0 );
    }
  }
  else {
    BAL_FreeDoublePointList( &theFloatingPoints );
    API_ErrorParse_pointmatching( _BaseName( argv[0] ), "no reference point list...\n", 0 );
  }

  theReferencePoints.unit =  par.points_unit;
  if ( par.reference_voxel.x > 0.0 ) {
    theReferencePoints.vx = par.reference_voxel.x;
  } else if ( par.voxel.x > 0.0 ) {
    theReferencePoints.vx = par.voxel.x;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theReferencePoints.vx = ptrtemplate->vx;
  }
  if ( par.reference_voxel.y > 0.0 ) {
    theReferencePoints.vy = par.reference_voxel.y;
  } else if ( par.voxel.y > 0.0 ) {
    theReferencePoints.vy = par.voxel.y;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theReferencePoints.vy = ptrtemplate->vy;
  }
  if ( par.reference_voxel.z > 0.0 ) {
    theReferencePoints.vz = par.reference_voxel.z;
  } else if ( par.voxel.z > 0.0 ) {
    theReferencePoints.vz = par.voxel.z;
  } else if ( ptrtemplate != (bal_image*)NULL ) {
    theReferencePoints.vz = ptrtemplate->vz;
  }





  /* API call
   */

  lineoptions = _Array2Str( argc, argv );
  if ( lineoptions == (char*)NULL )
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "unable to translate command line options ...\n", 0 );

  resTransformation = API_pointmatching( &theFloatingPoints,
                                         &theReferencePoints,
                                         ptrtemplate,
                                         par.result_residual,
                                         lineoptions, (char*)NULL );
   if ( resTransformation == (bal_transformation*)NULL ) {
      vtfree( lineoptions );
      BAL_FreeDoublePointList( &theReferencePoints );
      BAL_FreeDoublePointList( &theFloatingPoints );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "some error occurs during processing ...\n", -1 );

  }

  vtfree( lineoptions );


  /* writing result transformation
   */

  (void)BAL_InitFullImage( &imfloat, NULL, 1, 1, 1, 1,
                       theFloatingPoints.vx, theFloatingPoints.vy, theFloatingPoints.vy, UCHAR );
  (void)BAL_InitFullImage( &imreference, NULL, 1, 1, 1, 1,
                       theReferencePoints.vx, theReferencePoints.vy, theReferencePoints.vz, UCHAR );

  BAL_FreeDoublePointList( &theReferencePoints );
  BAL_FreeDoublePointList( &theFloatingPoints );

  if ( par.result_real_transformation != (char*)NULL && par.result_real_transformation[0] != '\0' ) {
    if ( BAL_WriteTransformation( resTransformation, par.result_real_transformation ) != 1 ) {
      BAL_FreeImage( &imreference );
      BAL_FreeImage( &imfloat );
      BAL_FreeTransformation( resTransformation );
      vtfree( resTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write real result transformation '%s'\n",
                 _BaseName( argv[0] ), par.result_real_transformation  );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "error when writing transformation ...\n", 0 );
    }
  }

  if ( par.result_voxel_transformation != (char*)NULL && par.result_voxel_transformation[0] != '\0' ) {
    if ( BAL_ChangeTransformationToVoxelUnit( &imfloat, &imreference,
                                              resTransformation, resTransformation ) != 1 ) {
      BAL_FreeImage( &imreference );
      BAL_FreeImage( &imfloat );
      BAL_FreeTransformation( resTransformation );
      vtfree( resTransformation );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "error when converting to voxel unit ...\n", 0 );
    }
    if ( BAL_WriteTransformation( resTransformation, par.result_voxel_transformation ) != 1 ) {
      BAL_FreeImage( &imreference );
      BAL_FreeImage( &imfloat );
      BAL_FreeTransformation( resTransformation );
      vtfree( resTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write voxel result transformation '%s'\n",
                 _BaseName( argv[0] ), par.result_voxel_transformation  );
      API_ErrorParse_pointmatching( _BaseName( argv[0] ), "error when writing transformation ...\n", 0 );
    }
  }

  BAL_FreeImage( &imreference );
  BAL_FreeImage( &imfloat );
  BAL_FreeTransformation( resTransformation );
  vtfree( resTransformation );

  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char *_Array2Str( int argc, char *argv[] )
{
  char *proc = "_Array2Str";
  int i, l;
  char *s, *t;

  if ( argc <= 1 || argv == (char**)NULL ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: no options in argv[]\n", proc );
    return( (char*)NULL );
  }

  /* there are argc-1 strings
   * compute the sum of string lengths from 1 to argc-1
   * + number of interval between successive strings (argc-2)
   * + 1 to add a trailing '\0'
   */
  for ( l=argc-1, i=1; i<argc; i++ ) {
    l += strlen( argv[i] );
  }

  s = (char*)vtmalloc( l * sizeof( char ), "s", proc );
  if ( s == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (char*)NULL );
  }

  for ( t=s, i=1; i<argc; i++ ) {
    (void)strncpy( t, argv[i], strlen( argv[i] ) );
    t += strlen( argv[i] );
    if ( i < argc-1 ) {
      *t = ' ';
      t++;
    }
    else {
      *t = '\0';
    }
  }

  return( s );
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}



static double _GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}



static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
