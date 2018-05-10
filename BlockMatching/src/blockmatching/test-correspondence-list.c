/*************************************************************************
 * test-pyramid-image.c -
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

#include <bal-interpolation.h>



int main(int argc, char *argv[])
{
  typeCorrespondenceList list;
  FILE *f;

  BAL_InitCorrespondenceList( &list );
  if ( argc <= 2 ) {
    fprintf( stderr, "must provide two args\n" );
    exit(-1);
  }

  f = fopen( argv[1], "r" );
  (void)BAL_ReadCorrespondenceList( f, &list );
  fclose( f );

  f = fopen( argv[2], "w" );
  BAL_PrintCorrespondenceList( f, &list );
  fclose( f );



  BAL_FreeCorrespondenceList( &list );
  return( 0 );
}
