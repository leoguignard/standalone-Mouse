/*************************************************************************
 * printTrsf.c -
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
#include <string.h>

#include <vtmalloc.h>

#include <bal-stddef.h>
#include <bal-transformation.h>









static char *program = NULL;

static char *usage = "%s ... %s\n\
[-transformation-type|-transformation|-trsf-type %s]\n\
[-module %s]";

static char *detail = "";



static int _verbose_ = 1;


static void _ErrorParse( char *str, int flag );
static char *_BaseName( char *p );





int main(int argc, char *argv[])
{
  int i;
  int *isatrsf = (int*)NULL;
  bal_transformation theTrsf;
  bal_image theModule;
  char *nameModule = (char*)NULL;
  enumTypeTransfo type = UNDEF_TRANSFORMATION;



  BAL_InitImage( &theModule, NULL, 0, 0, 0, 0, UCHAR );


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  isatrsf = (int*)vtmalloc( argc * sizeof(int), "isatrsf", argv[0] );
  if ( isatrsf == (int*)NULL ) {
    _ErrorParse( "allocation failed\n", 1 );
  }
  for ( i=0; i<argc; i++ )  isatrsf[i] = 1;
  isatrsf[0] = 0;

  
  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );

  
  /* displaying help is required
   */
  i = 1;
  while ( i < argc ) {
    if ( ( strcmp ( argv[i], "-help") == 0 ) 
        || ( strcmp ( argv[i], "-h") == 0 && argv[i][2] == '\0' )
        || ( strcmp ( argv[i], "--help") == 0 )
        || ( strcmp ( argv[i], "--h") == 0 && argv[i][3] == '\0' ) ) {
      _ErrorParse( NULL, 1 );
    }
    else if ( strcmp ( argv[i], "-verbose" ) == 0
              || ( strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0' ) ) {
      isatrsf[i] = 0;
      BAL_IncrementVerboseInBalTransformation();
    }
    else if ( strcmp ( argv[i], "-no-verbose" ) == 0
              || strcmp ( argv[i], "-noverbose" ) == 0
              || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
      isatrsf[i] = 0;
      BAL_SetVerboseInBalImage( 0 );
    }

    else if ( strcmp ( argv[i], "-transformation-type" ) == 0
              || strcmp ( argv[i], "-transformation" ) == 0
              || (strcmp ( argv[i], "-trsf-type" ) == 0 && argv[i][10] == '\0') ) {
       isatrsf[i] = 0;
       i ++;
       if ( i >= argc)
           _ErrorParse( "parsing -transformation-type", 0 );
       if ( strcmp ( argv[i], "translation2D" ) == 0 ) {
         type = TRANSLATION_2D;
       }
       else if ( strcmp ( argv[i], "translation3D" ) == 0 ) {
         type = TRANSLATION_3D;
       }
       else if ( strcmp ( argv[i], "translation" ) == 0 && argv[i][11] == '\0') {
         type = TRANSLATION_3D;
       }
       else if ( strcmp ( argv[i], "translation-scaling2D" ) == 0 ) {
         type = TRANSLATION_SCALING_2D;
       }
       else if ( strcmp ( argv[i], "translation-scaling3D" ) == 0 ) {
         type = TRANSLATION_SCALING_3D;
       }
       else if ( strcmp ( argv[i], "rigid2D" ) == 0 ) {
         type = RIGID_2D;
       }
       else if ( strcmp ( argv[i], "rigid3D" ) == 0 ) {
         type = RIGID_3D;
       }
       else if ( (strcmp ( argv[i], "rigid" ) == 0 && argv[i][5] == '\0') ) {
         type = RIGID_3D;
       }
       else if ( strcmp ( argv[i], "similitude2D" ) == 0 ) {
         type = SIMILITUDE_2D;
       }
       else if ( strcmp ( argv[i], "similitude3D" ) == 0 ) {
         type = SIMILITUDE_3D;
       }
       else if ( strcmp ( argv[i], "similitude" ) == 0 ) {
         type = SIMILITUDE_3D;
       }
       else if ( strcmp ( argv[i], "affine2D" ) == 0 ) {
         type = AFFINE_2D;
       }
       else if ( strcmp ( argv[i], "affine3D" ) == 0 ) {
         type = AFFINE_3D;
       }
       else if ( strcmp ( argv[i], "affine" ) == 0 ) {
         type = AFFINE_3D;
       }
       /*
         else if ( strcmp ( argv[i], "spline" ) == 0 ) {
         type = SPLINE;
         }
       */
       else if ( strcmp ( argv[i], "vectorfield" ) == 0
                 || strcmp ( argv[i], "vector" ) == 0 ) {
         type = VECTORFIELD_3D;
       }
       else if ( strcmp ( argv[i], "vectorfield3D" ) == 0
                 || strcmp ( argv[i], "vector3D" ) == 0 ) {
         type = VECTORFIELD_3D;
       }
       else if ( strcmp ( argv[i], "vectorfield2D" ) == 0
                 || strcmp ( argv[i], "vector2D" ) == 0 ) {
         type = VECTORFIELD_2D;
       }
       else {
         fprintf( stderr, "unknown transformation type: '%s'\n", argv[i] );
         _ErrorParse( "parsing -transformation-type", 0 );
       }
       isatrsf[i] = 0;
    }
    else if ( strcmp ( argv[i], "-module" ) == 0 ) {
       isatrsf[i] = 0;
       i ++;
       if ( i >= argc )
           _ErrorParse( "parsing -module", 0 );
       nameModule = argv[i];
       isatrsf[i] = 0;
    }

    i++;
  }



  /***************************************************
   *
   * 
   *
   ***************************************************/
  BAL_InitTransformation( &theTrsf );

  for ( i=1; i<argc; i++ ) {

    if ( isatrsf[i] == 0 ) continue;

    if ( BAL_ReadTransformation( &theTrsf, argv[i] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read '%s'\n", program, argv[i] );
      fprintf( stderr, "\n" );
    }

    if ( type != UNDEF_TRANSFORMATION )
        theTrsf.type = type;
    
    fprintf( stdout, "***** Information about '%s'\n\n", argv[i] );
    BAL_PrintTransformation( stdout, &theTrsf, argv[i] );
    fprintf( stdout, "\n" );

    if ( nameModule != (char*)NULL ) {
      switch ( theTrsf.type ) {
      default:
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to compute transformation displacements\n", program );
          fprintf( stderr, "\t such transformation type not handled yet\n" );
        }
        break;
      case VECTORFIELD_2D :
      case VECTORFIELD_3D :
        if ( BAL_AllocImageFromImage( &theModule, (char *)NULL, &(theTrsf.vx), FLOAT ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: can not allocate modulus image\n", program );
          break;
        }
        if ( BAL_TransformationModulus(  &theTrsf, &theModule ) != 1 ) {
          BAL_FreeImage( &theModule );
          if ( _verbose_ )
            fprintf( stderr, "%s: can not compute modulus image\n", program );
          break;
        }
        if ( BAL_WriteImage( &theModule, nameModule ) != 1 ) {
          BAL_FreeImage( &theModule );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to write module image '%s'\n", program, nameModule );
          break;
        }
        BAL_FreeImage( &theModule );
        break;
      }
    } /* if ( nameModule != (char*)NULL ) */

    BAL_FreeTransformation( &theTrsf );
  }

  return( 0 );
}




/***************************************************
 *
 * 
 *
 ***************************************************/





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage: %s %s\n",_BaseName(program), usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s\n",detail);
  if ( str != NULL ) (void)fprintf(stderr,"Error: %s\n",str);
  exit( 1 );
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

