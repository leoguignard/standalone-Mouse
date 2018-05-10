/*************************************************************************
 * threshold.c - image thresholding
 *
 * $$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 17 nov 2016 09:17:49 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <threshold.h>

static int _verbose_ = 1;



static int _thresholdBuffer( void *bufferIn,
                     bufferType typeIn,
                     void *bufferOut,
                     bufferType typeOut,
                     size_t bufferLength,
                     float threshold )
{
  char *proc = "_thresholdBuffer";
  size_t i;
  int ithreshold;

#define _THRESHOLD( itype, otype, t, v ) {               \
  itype *theBuf = (itype*)bufferIn;                      \
  otype *resBuf = (otype*)bufferOut;                     \
  for ( i=0; i<bufferLength; i++, theBuf++, resBuf++ ) { \
    *resBuf = (otype)(( *theBuf >= (itype)t ) ? v : 0);  \
  }                                                      \
}

#define _SWITCH_THRESHOLD( otype, v ) {                                 \
  switch( typeIn ) {                                                    \
  default :                                                             \
    if ( _verbose_ )                                                    \
      fprintf( stderr, "%s: such input type not handled yet\n", proc ); \
    return( -1 );                                                       \
  case UCHAR :                                                          \
    _THRESHOLD( u8, otype, ithreshold, v );                             \
    break;                                                              \
  case SCHAR :                                                          \
    _THRESHOLD( s8, otype, ithreshold, v );                             \
    break;                                                              \
  case USHORT :                                                         \
    _THRESHOLD( u16, otype, ithreshold, v );                            \
    break;                                                              \
  case SSHORT :                                                         \
    _THRESHOLD( s16, otype, ithreshold, v );                            \
    break;                                                              \
  case UINT :                                                           \
    _THRESHOLD( u32, otype, ithreshold, v );                            \
    break;                                                              \
  case SINT :                                                           \
    _THRESHOLD( s32, otype, ithreshold, v );                            \
    break;                                                              \
  case FLOAT :                                                          \
    _THRESHOLD( r32, otype, threshold, v );                             \
    break;                                                              \
  case DOUBLE :                                                         \
    _THRESHOLD( r64, otype, threshold, v );                             \
    break;                                                              \
  }                                                                     \
}

  if ( threshold >= 0.0 ) ithreshold = (int)(threshold + 0.5);
  else ithreshold = (int)(threshold - 0.5);

  switch( typeOut ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
    _SWITCH_THRESHOLD( u8, 255 );
    break;
  case SCHAR :
    _SWITCH_THRESHOLD( s8, 127 );
    break;
  case USHORT :
    _SWITCH_THRESHOLD( u16, 65535 );
    break;
  case SSHORT :
    _SWITCH_THRESHOLD( s16, 32767 );
    break;
  case UINT :
    _SWITCH_THRESHOLD( u32, 1 );
    break;
  case SINT :
    _SWITCH_THRESHOLD( s32, 1 );
    break;
  case FLOAT :
    _SWITCH_THRESHOLD( r32, 1.0 );
    break;
  case DOUBLE :
    _SWITCH_THRESHOLD( r64, 1.0 );
    break;

  }

  return( 1 );

}





int thresholdBuffer( void *bufferIn,
                     bufferType typeIn,
                     void *bufferOut,
                     bufferType typeOut,
                     int *bufferDims,
                     float threshold )
{
  size_t bufferLength = 0;
  bufferLength = (size_t)bufferDims[0] * (size_t)bufferDims[1] * (size_t)bufferDims[2];
  return( _thresholdBuffer( bufferIn, typeIn, bufferOut, typeOut, bufferLength, threshold) );
}
