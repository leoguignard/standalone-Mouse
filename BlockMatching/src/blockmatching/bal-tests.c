/*************************************************************************
 * bal-tests.c -
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



/* random(), srandom(): Feature Test Macro Requirements for glibc
 * _SVID_SOURCE || _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED
 *
 * compilation with [gcc (GCC) 5.3.1 20151207 (Red Hat 5.3.1-2)] yields
 * "_BSD_SOURCE and _SVID_SOURCE are deprecated, use _DEFAULT_SOURCE"
 */
#define _DEFAULT_SOURCE



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vtmalloc.h>

#include <bal-transformation-tools.h>
#include <bal-pyramid.h>

/*--------------------------------------------------*
 *
 * TESTS
 *
 *--------------------------------------------------*/

void BAL_TestPyramidLevel( int n ) 
{
  int x, y, z, i;
  int m = 1000;
  bal_pyramid_level *pyramid_level = NULL;
  int l, nl, levels;
  bal_image theIm;
  bal_blockmatching_pyramidal_param p;

  srandom( time(0) );

  levels = (int)( log((double)(m)) / log(2.0) ) + 2;
  pyramid_level = (bal_pyramid_level *)vtmalloc( (levels) * sizeof(bal_pyramid_level),
                                                 "pyramid_level", "BAL_TestPyramidLevel" );

  BAL_InitBlockMatchingPyramidalParameters( &p );
  p.pyramid_lowest_level = 0;;
  p.pyramid_highest_level = levels-1;
  p.blocks_fraction.lowest = 0.25;
  p.blocks_fraction.highest = 0.5;

  for (i=0; i<n; i++ ) {
    
    x = (int)( m * (double)random() / (double)RAND_MAX );
    y = (int)( m * (double)random() / (double)RAND_MAX );
    z = (int)( m * (double)random() / (double)RAND_MAX );
    (void)BAL_AllocFullImage( &theIm, "test.inr.gz", x, y, z, 1, 1.0, 1.0, 1.0, UCHAR );
    nl = _ComputePyramidLevel( pyramid_level, levels-1, &theIm, &p );

    fprintf( stdout, "\n" );
    fprintf( stdout, "TEST #%d : x=%4d y=%4d z=%4d\n", i, x, y, z );
    if ( nl < 0 ) {
      fprintf( stdout, "- no pyramid level\n" );
    }
    else { 
      for ( l=0; l<=nl; l++ ) 
        _PrintPyramidLevel( stdout, &(pyramid_level[l]), 1 );
    }
    BAL_FreeImage( &theIm );
  }
  
}




/* 
   - pyramid_gaussian_filtering
     the pyramid is built with smooting
   - test_inverse_trsf
     test the consistance between direct and inverse transformation
 */
void BAL_TestPyramidImage( bal_image *theIm, 
                           int pyramid_gaussian_filtering )
{
  int m;
  bal_pyramid_level *pyramid_level = NULL;
  int l, nl, levels;
  bal_image resIm;
  bal_blockmatching_pyramidal_param p;
  char imagename[256];


  srandom( time(0) );

  if ( theIm->nplanes == 1 )
    m = max ( theIm->ncols, theIm->nrows );
  else 
    m = max ( theIm->ncols, max( theIm->nrows, theIm->nplanes ) );

  levels = (int)( log((double)(m)) / log(2.0) ) + 2;
  pyramid_level = (bal_pyramid_level *)vtmalloc( (levels) * sizeof(bal_pyramid_level),
                                                 "pyramid_level", "BAL_TestPyramidImage" );

  BAL_InitBlockMatchingPyramidalParameters( &p );
  p.pyramid_lowest_level = 0;
  p.pyramid_highest_level = levels-1;
  p.blocks_fraction.lowest = 0.25;
  p.blocks_fraction.highest = 0.5;

  if ( theIm->nplanes == 1 ) {
    p.block_dim.z = 1;
  }

  nl = _ComputePyramidLevel( pyramid_level, levels-1, theIm, &p );
  
  if ( nl < 0 ) {
    fprintf( stdout, "- no pyramid level\n" );
  }
  else { 

    BAL_InitImage( &resIm, NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
    for ( l=0; l<=nl; l++ ) {
      _PrintPyramidLevel( stdout, &(pyramid_level[l]), 1 ); 

      sprintf( imagename, "imageatlevel%02d.hdr", l );

      (void)BAL_AllocComputeSubsampledImage( &resIm, pyramid_level[l].ncols,
                                       pyramid_level[l].nrows, pyramid_level[l].nplanes,
                                       theIm,
                                       &(pyramid_level[l]), pyramid_gaussian_filtering );


      (void)BAL_WriteImage( &resIm, imagename );
      BAL_FreeImage( &resIm );
    }

  }
}
