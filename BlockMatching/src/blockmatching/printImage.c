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
#include <bal-image.h>









static char *program = NULL;

static char *usage = "%s ... %s\n";

static char *detail = "";



static int _verbose_ = 1;


static void _ErrorParse( char *str, int flag );
static char *_BaseName( char *p );


typedef enum _enumFormatPrint {
  _DEFAULT_FORMAT_PRINT_,
  _ZPAR_FORMAT_PRINT_
} _enumFormatPrint;


int main(int argc, char *argv[])
{
  int i;
  int *isanimage = (int*)NULL;
  bal_image theIm;
  _enumFormatPrint formatPrint = _DEFAULT_FORMAT_PRINT_;


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  isanimage  = (int*)vtmalloc( argc * sizeof(int), "isanimage ", argv[0] );
  if ( isanimage  == (int*)NULL ) {
    _ErrorParse( "allocation failed\n", 1 );
  }
  for ( i=0; i<argc; i++ )  isanimage [i] = 1;
  isanimage [0] = 0;

  
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
      isanimage [i] = 0;
      BAL_IncrementVerboseInBalImage();
    }
    else if ( strcmp ( argv[i], "-no-verbose" ) == 0
              || strcmp ( argv[i], "-noverbose" ) == 0
              || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
      isanimage [i] = 0;
      BAL_SetVerboseInBalImage( 0 );
    }
    else if ( (strcmp ( argv[i], "-par" ) == 0 && argv[i][4] == '\0')
              || (strcmp ( argv[i], "-zpar" ) == 0 && argv[i][5] == '\0') ) {
      isanimage [i] = 0;
      formatPrint = _ZPAR_FORMAT_PRINT_;
    }
    i++;
  }



  /***************************************************
   *
   * 
   *
   ***************************************************/

  for ( i=1; i<argc; i++ ) {

    if ( isanimage [i] == 0 ) continue;

    if ( BAL_ReadImage( &theIm, argv[i], 0 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to read '%s'\n", program, argv[i] );
      fprintf( stderr, "\n" );
    }
    
    switch( formatPrint ) {
    default :
    case _DEFAULT_FORMAT_PRINT_ :
        fprintf( stdout, "***** Information about '%s'\n\n", argv[i] );
        BAL_PrintImage( stdout, &theIm, argv[i] );
        fprintf( stdout, "\n" );
        break;
    case _ZPAR_FORMAT_PRINT_ :
        BAL_PrintParImage( stdout, &theIm, argv[i] );
        break;
    }

    BAL_FreeImage( &theIm );
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

