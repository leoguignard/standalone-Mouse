/*************************************************************************
 * test-sizeof.c - 
 *
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 11:08:39 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>


int main (int argc, char *argv[] )
{
 
  fprintf( stderr, "sizeof( size_t ) = %lu\n", sizeof( size_t ) );

  return( 0 );
}


