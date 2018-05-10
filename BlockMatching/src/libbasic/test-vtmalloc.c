/*************************************************************************
 * test-vtmalloc.c -
 *
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 16 oct 2016 22:27:01 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdio.h>

#include <vtmalloc.h>

#define TEST 10

int main (int argc, char *argv[] )
{
  int i;
  void *p[TEST];
  char s[20];

  setTraceInVtMalloc( 1 );

  for ( i=0; i<TEST; i++ ) {
    sprintf( s, "p[%d]", i );
    p[i] = vtmalloc( 100000, s, "main" );
  }

  fprintfVtMallocTrace( stdout );

  for ( i=0; i<TEST; i++ ) {
    vtfree( p[i] );
  }

  fprintfVtMallocTrace( stdout );

  clearVtMalloc();

  return( 0 );
}


