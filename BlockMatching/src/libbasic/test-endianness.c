/****************************************************
 * test-endianness.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 14:26:30 CET 2012
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>


int _getEndianness()
{
  union {
    unsigned char uc[2];
    unsigned short us;
  } twobytes;
  twobytes.us = 255;
  /* little endian
   */
  if ( twobytes.uc[1] == 0 ) return( 1 );
  /* big endian
   */
  return( 2);
}

int get_byteorder()
{
        union {
                long l;
                char c[4];
        } test;
        test.l = 1;
        if( test.c[3] && !test.c[2] && !test.c[1] && !test.c[0] )
            /* big endian
             */
                return 2;

        if( !test.c[3] && !test.c[2] && !test.c[1] && test.c[0] )
            /* little endian
             */
                return 1;

        return 0;
}



int main (int argc, char const *argv[])
{
    fprintf( stderr, "test #1:" );
    switch( _getEndianness() ) {
    default : fprintf( stderr, "default\n" ); break;
    case 1 : fprintf( stderr, "little endian\n" ); break;
    case 2 : fprintf( stderr, "big endian\n" ); break;
    }
    fprintf( stderr, "test #2:" );
    switch( get_byteorder() ) {
    default : fprintf( stderr, "default\n" ); break;
    case 0 : fprintf( stderr, "little endian\n" ); break;
    case 1 : fprintf( stderr, "little endian\n" ); break;
    case 2 : fprintf( stderr, "big endian\n" ); break;
    }

}


