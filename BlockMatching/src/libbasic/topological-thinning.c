/****************************************************
 * topological-thinning.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mer  2 mar 2016 22:02:15 CET
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

#include <time.h>
#include <string.h>

#include <chamferdistance.h>
#include <t04t08.h>
#include <t06t26.h>

#include <vtmalloc.h>

#include <topological-operations-common.h>
#include <topological-thinning.h>





static void _ComputeTimeFrom2Clocks( float c1, float c2,
                                     int *hours, int *minutes, float *seconds );
static void _PrintTimeFrom2Clocks( char *str, float c1, float c2 );





/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;
static int _time_ = 0;

void setVerboseInTopologicalThinning( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalThinning(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalThinning(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalThinning( int v )
{
  _debug_ = v;
}

void incrementDebugInTopologicalThinning(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalThinning(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}

void setTimeInTopologicalThinning( int v )
{
  _time_ = v;
}








/**************************************************
 *
 * parameter management
 *
 **************************************************/



void initTypeThinningParameters( typeThinningParameters *p )
{
  /* list construction: bin size
   */
  p->binLength = 1;

  /* distance-based thinning conditions
   */
  p->typeOrdering = _INCREASING_ORDER_;
  p->typeChanging = _FOREGROUND_TO_BACKGROUND_;

  /* points to be changed
   */
  p->typeThickness = _26THICKNESS_;
  p->typeEndPoint = _SURFACE_;

  /* test conditions before end condition
   */
  p->cyclesBeforeEnding = -1;
  p->valueBeforeEnding = -1;

  /* end condition
   */
  p->maxIteration = -1;
}

static void _fprintfTypeThinningParameters( FILE *f, typeThinningParameters *par )
{
    fprintf( f, " ... bin length           = %d\n", par->binLength );
    fprintf( f, " ... value spanning order = " );
    switch( par->typeOrdering ) {
    default : fprintf( f, "unknown" );
    case _INCREASING_ORDER_ : fprintf( f, "_INCREASING_ORDER_" ); break;
    case _DECREASING_ORDER_ : fprintf( f, "_DECREASING_ORDER_" ); break;
    }
    fprintf( f, "\n" );
    fprintf( f, " ... changing type        = " );
    switch( par->typeChanging ) {
    default : fprintf( f, "unknown" );
    case _FOREGROUND_TO_BACKGROUND_ : fprintf( f, "_FOREGROUND_TO_BACKGROUND_" ); break;
    case _BACKGROUND_TO_FOREGROUND_ : fprintf( f, "_BACKGROUND_TO_FOREGROUND_" ); break;
    }
    fprintf( f, "\n" );
    fprintf( f, " ... typeThickness        = " );
    switch ( par->typeThickness ) {
    default : fprintf( f, "unknown" );
    case _04THICKNESS_ : fprintf( f, "04THICKNESS" ); break;
    case _08THICKNESS_ : fprintf( f, "08THICKNESS" ); break;
    case _06THICKNESS_ : fprintf( f, "06THICKNESS" ); break;
    case _18THICKNESS_ : fprintf( f, "18THICKNESS" ); break;
    case _26THICKNESS_ : fprintf( f, "26THICKNESS" ); break;
    }
    fprintf( f, "\n" );
    fprintf( f, " ... typeEndPoint         = " );
    switch ( par->typeEndPoint ) {
    default : fprintf( f, "unknown" );
    case _SURFACE_      : fprintf( f, "SURFACE" ); break;
    case _PURE_SURFACE_ : fprintf( f, "SURFACE" ); break;
    case _CURVE_        : fprintf( f, "CURVE" ); break;
    case _PURE_CURVE_   : fprintf( f, "CURVE" ); break;
    case _NO_END_POINT_ : fprintf( f, "NO_END_POINT" ); break;
    }
    fprintf( f, "\n" );
    fprintf( f, " ... cyclesBeforeEnding   = %d\n", par->cyclesBeforeEnding );
    fprintf( f, " ... valueBeforeEnding    = %d\n", par->valueBeforeEnding );
    fprintf( f, " ... max. iterations      = %d\n", par->maxIteration );

}












/**************************************************
 *
 * interiority tests
 *
 **************************************************/



typedef int (*typeIsInsideFunction)( int, int, int, int * );



/* __attribute__ ((unused))
 * allows to suppress the "warning: unused parameter ‘z’"
 * when compiled with gcc
 */
static int _is2DInside( int x, int y, int z __attribute__ ((unused)), int *theDim )
{
  if ( 0 < x && x < theDim[0]-1 &&
       0 < y && y < theDim[1]-1 )
    return( 1 );
  return( 0 );
}



static int _is3DInside( int x, int y, int z, int *theDim )
{
  if ( 0 < x && x < theDim[0]-1 &&
       0 < y && y < theDim[1]-1 &&
       0 < z && z < theDim[2]-1 )
    return( 1 );
  return( 0 );
}










/**************************************************
 *
 * thickness tests
 *
 **************************************************/



typedef int (*typeCheckThickness)( int *,
                                   enumThickness,
                                   enumTypeChange,
                                   int );



static int _check2DThickness( int *neighb,
                              enumThickness thickness,
                              enumTypeChange typeChanging,
                              int direction )
{
  int t;
  /* neighborhood
     0 1 2
     3 4 5
     6 7 8
  */
  /* check the thickness ?
   * FOREGROUND TO BACKGROUND (thinnig)
   * 1. foreground = (points > 0)
   *    (ANCHOR=FOREGROUND=255, CANBECHANGED, TOBECHANGED)
   *    background = (points = 0)
   *    (BACKGROUND=0)
   *    => FOREGROUND > _BACKGROUND_VALUE_
   *       BACKGROUND <= _BACKGROUND_VALUE_
   *
   * BACKGROUND TO FOREGROUND (thickening)
   * 1. foreground = (ANCHOR=FOREGROUND)
   *    background = (BACKGROUND, CANBECHANGED, TOBECHANGED)
   *    => FOREGROUND > _ENDPOINT_VALUE_
   *       BACKGROUND <= _ENDPOINT_VALUE_
   *
   */
  /* est-ce une epaisseur ?
     0 = . 0 .   1 = . 1 .   2 = . . .  3 = . . .
         . x .       . x .       1 x 0      0 x 1
         . 1 .       . 0 .       . . .      . . .
  */

  switch( typeChanging ) {
  default :
  case _FOREGROUND_TO_BACKGROUND_ :
      t = _BACKGROUND_VALUE_;
      break;
  case _BACKGROUND_TO_FOREGROUND_ :
      t = _ENDPOINT_VALUE_;
      break;
}

  switch ( thickness ) {
  case _04THICKNESS_ :
  case _06THICKNESS_ :
    switch( direction ) {
    default :
    case 0 :
      if ( neighb[1]  > t )  return( 0 );
      if ( neighb[7] <= t )  return( 0 );
      break;
    case 1 :
      if ( neighb[7]  > t )  return( 0 );
      if ( neighb[1] <= t )  return( 0 );
      break;
    case 2 :
      if ( neighb[5]  > t )  return( 0 );
      if ( neighb[3] <= t )  return( 0 );
      break;
    case 3 :
      if ( neighb[3]  > t )  return( 0 );
      if ( neighb[5] <= t )  return( 0 );
      break;
    }
  case _08THICKNESS_ :
  case _18THICKNESS_ :
  case _26THICKNESS_ :
  default :
    switch( direction ) {
    default :
    case 0 :
      if ( neighb[1]  > t )  return( 0 );
      if ( neighb[6] <= t && neighb[7] <= t && neighb[8] <= t )  return( 0 );
      break;
    case 1 :
      if ( neighb[7]  > t )  return( 0 );
      if ( neighb[0] <= t && neighb[1] <= t && neighb[2] <= t )  return( 0 );
      break;
    case 2 :
      if ( neighb[5]  > t )  return( 0 );
      if ( neighb[0] <= t && neighb[3] <= t && neighb[6] <= t )  return( 0 );
      break;
    case 3 :
      if ( neighb[3]  > t )  return( 0 );
      if ( neighb[2] <= t && neighb[5] <= t && neighb[8] <= t )  return( 0 );
      break;
    }
    break;
  }
  return( 1 );
}



static int _check3DThickness( int *neighb,
                              enumThickness thickness,
                              enumTypeChange typeChanging,
                              int direction )
{
  int t;
  /* neighborhood
     0  1  2     9 10 11    18 19 20
     3  4  5    12 13 14    21 22 23
     6  7  8    15 16 17    24 25 26
  */
    /* check the thickness ?
     * FOREGROUND TO BACKGROUND (thinnig)
     * 1. foreground = (points > 0)
     *    (ANCHOR=FOREGROUND=255, CANBECHANGED, TOBECHANGED)
     *    background = (points = 0)
     *    (BACKGROUND=0)
     *    => FOREGROUND > _BACKGROUND_VALUE_
     *       BACKGROUND <= _BACKGROUND_VALUE_
     *
     * BACKGROUND TO FOREGROUND (thickening)
     * 1. foreground = (ANCHOR=FOREGROUND)
     *    background = (BACKGROUND, CANBECHANGED, TOBECHANGED)
     *    => FOREGROUND > _ENDPOINT_VALUE_
     *       BACKGROUND <= _ENDPOINT_VALUE_
     *
     */

  switch( typeChanging ) {
  default :
  case _FOREGROUND_TO_BACKGROUND_ :
      t = _BACKGROUND_VALUE_;
      break;
  case _BACKGROUND_TO_FOREGROUND_ :
      t = _ENDPOINT_VALUE_;
      break;
  }

  switch ( thickness ) {
  case _04THICKNESS_ :
  case _06THICKNESS_ :

    switch( direction ) {
    default :
    case 0 :
      if ( neighb[ 4]  > t )  return( 0 );
      if ( neighb[22] <= t )  return( 0 );
              break;
    case 1 :
      if ( neighb[22]  > t )  return( 0 );
      if ( neighb[ 4] <= t )  return( 0 );
      break;
    case 2 :
      if ( neighb[10]  > t )  return( 0 );
      if ( neighb[16] <= t )  return( 0 );
      break;
    case 3 :
      if ( neighb[16]  > t )  return( 0 );
      if ( neighb[10] <= t )  return( 0 );
      break;
    case 4 :
      if ( neighb[12]  > t )  return( 0 );
      if ( neighb[14] <= t )  return( 0 );
      break;
    case 5 :
      if ( neighb[14]  > t )  return( 0 );
      if ( neighb[12] <= t )  return( 0 );
      break;
    }
    break;

  case _18THICKNESS_ :

    switch( direction ) {
    default :
    case 0 :
      if ( neighb[ 4]  > t )  return( 0 );
      if ( neighb[22] <= t &&
           neighb[19] <= t && neighb[21] <= t &&
           neighb[23] <= t && neighb[25] <= t )  return( 0 );
      break;
    case 1 :
      if ( neighb[22]  > t )  return( 0 );
      if ( neighb[ 4] <= t &&
           neighb[ 1] <= t && neighb[ 3] <= t &&
           neighb[ 5] <= t && neighb[ 7] <= t )  return( 0 );
      break;
    case 2 :
      if ( neighb[10]  > t )  return( 0 );
      if ( neighb[16] <= t &&
           neighb[ 7] <= t && neighb[15] <= t &&
           neighb[17] <= t && neighb[25] <= t )  return( 0 );
      break;
    case 3 :
      if ( neighb[16]  > t )  return( 0 );
      if ( neighb[10] <= t &&
           neighb[ 1] <= t && neighb[ 9] <= t &&
           neighb[11] <= t && neighb[19] <= t )  return( 0 );
      break;
    case 4 :
      if ( neighb[12]  > t )  return( 0 );
      if ( neighb[14] <= t &&
           neighb[ 5] <= t && neighb[11] <= t &&
           neighb[17] <= t && neighb[23] <= t )  return( 0 );
      break;
    case 5 :
      if ( neighb[14]  > t )  return( 0 );
      if ( neighb[12] <= t &&
           neighb[ 5] <= t && neighb[ 9] <= t &&
           neighb[15] <= t && neighb[21] <= t )  return( 0 );
      break;
    }
    break;

  case _08THICKNESS_ :
  case _26THICKNESS_ :
  default :

    switch( direction ) {
    default :
    case 0 :
      if ( neighb[ 4]  > t )  return( 0 );
      if ( neighb[22] <= t &&
           neighb[19] <= t && neighb[21] <= t &&
           neighb[23] <= t && neighb[25] <= t &&
           neighb[18] <= t && neighb[20] <= t &&
           neighb[24] <= t && neighb[26] <= t )  return( 0 );
      break;
    case 1 :
      if ( neighb[22]  > t )  return( 0 );
      if ( neighb[ 4] <= t &&
           neighb[ 1] <= t && neighb[ 3] <= t &&
           neighb[ 5] <= t && neighb[ 7] <= t &&
           neighb[ 0] <= t && neighb[ 2] <= t &&
           neighb[ 6] <= t && neighb[ 8] <= t )  return( 0 );
      break;
    case 2 :
      if ( neighb[10]  > t )  return( 0 );
      if ( neighb[16] <= t &&
           neighb[ 7] <= t && neighb[15] <= t &&
           neighb[17] <= t && neighb[25] <= t &&
           neighb[ 6] <= t && neighb[ 8] <= t &&
           neighb[24] <= t && neighb[26] <= t )  return( 0 );
      break;
    case 3 :
      if ( neighb[16]  > t )  return( 0 );
      if ( neighb[10] <= t &&
           neighb[ 1] <= t && neighb[ 9] <= t &&
           neighb[11] <= t && neighb[19] <= t &&
           neighb[ 0] <= t && neighb[ 2] <= t &&
           neighb[18] <= t && neighb[20] <= t )  return( 0 );
      break;
    case 4 :
      if ( neighb[12]  > t )  return( 0 );
      if ( neighb[14] <= t &&
           neighb[ 5] <= t && neighb[11] <= t &&
           neighb[17] <= t && neighb[23] <= t &&
           neighb[ 2] <= t && neighb[ 8] <= t &&
           neighb[20] <= t && neighb[26] <= t )  return( 0 );
      break;
    case 5 :
      if ( neighb[14]  > t )  return( 0 );
      if ( neighb[12] <= t &&
           neighb[ 5] <= t && neighb[ 9] <= t &&
           neighb[15] <= t && neighb[21] <= t &&
           neighb[ 0] <= t && neighb[ 6] <= t &&
           neighb[18] <= t && neighb[24] <= t )  return( 0 );
      break;
    }
    break;

  }
  return( 1 );
}










/**************************************************
 *
 * simplicity tests
 *
 * FOREGROUND TO BACKGROUND (thinnig)
 * 1. foreground = (points > 0)
 *    (ANCHOR=FOREGROUND=255, CANBECHANGED, TOBECHANGED)
 *    background = (points = 0)
 *    (BACKGROUND=0)
 * 2. with changes
 *    foreground = (ANCHOR, CANBECHANGED)
 *    background = (BACKGROUND, TOBECHANGED)
 *
 * BACKGROUND TO FOREGROUND (thickening)
 * 2. foreground = (ANCHOR=FOREGROUND)
 *    background = (BACKGROUND, CANBECHANGED, TOBECHANGED)
 * 1. foreground = (ANCHOR=FOREGROUND, TOBECHANGED)
 *    background = (BACKGROUND, CANBECHANGED)
 *
 *
 **************************************************/



typedef int (*typeIsPointSimple)( int *,
                                  int *,
                                  int *,
                                  enumTypeChange typeChanging );



static int _is2DPointSimple( int *neighb, int *t04, int *t08, enumTypeChange typeChanging )
{
  int checkT04, checkT08;
  int n;

  switch( typeChanging ) {
  default :
  case _FOREGROUND_TO_BACKGROUND_ :

      Compute_T04_and_T08( neighb, t04, t08 );
      if ( _debug_ ) {
          fprintf( stderr, " t04=%d, t08=%d\n", *t04, *t08 );
          fprintfNeighborhood( stderr, neighb, 9 );
      }
      if ( *t04 != 1 || *t08 != 1 ) return( 0 );

      for ( n=0; n<9; n++ )
        if ( neighb[n] == _TOBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T04_and_T08( neighb, &checkT04, &checkT08 );
      if ( _debug_ ) {
          fprintf( stderr, " t04=%d, t08=%d\n", checkT04, checkT08 );
          fprintfNeighborhood( stderr, neighb, 9 );
      }
      if ( checkT04 != 1 || checkT08 != 1 ) return( 0 );

      break;

  case _BACKGROUND_TO_FOREGROUND_ :

      for ( n=0; n<9; n++ )
        if ( neighb[n] == _CANBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T04_and_T08( neighb, &checkT04, &checkT08 );
      if ( _debug_ ) {
          fprintf( stderr, " t04=%d, t08=%d\n", checkT04, checkT08 );
          fprintfNeighborhood( stderr, neighb, 9 );
      }
      if ( checkT04 != 1 || checkT08 != 1 ) return( 0 );

      for ( n=0; n<9; n++ )
        if ( neighb[n] == _TOBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T04_and_T08( neighb, t04, t08 );
      if ( _debug_ ) {
          fprintf( stderr, " t04=%d, t08=%d\n", *t04, *t08 );
          fprintfNeighborhood( stderr, neighb, 9 );
      }
      if ( *t04 != 1 || *t08 != 1 ) return( 0 );

      break;
  }

  return( 1 );
}



static int _is3DPointSimple( int *neighb, int *t06, int *t26, enumTypeChange typeChanging )
{
  int checkT06, checkT26;
  int n;


  switch( typeChanging ) {
  default :
  case _FOREGROUND_TO_BACKGROUND_ :

      Compute_T06_and_T26( neighb, t06, t26 );
      if ( *t06 != 1 || *t26 != 1 ) return( 0 );

      for ( n=0; n<27; n++ )
        if ( neighb[n] == _TOBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T06_and_T26( neighb, &checkT06, &checkT26 );
      if ( checkT06 != 1 || checkT26 != 1 ) return( 0 );

      break;

  case _BACKGROUND_TO_FOREGROUND_ :

      for ( n=0; n<27; n++ )
        if ( neighb[n] == _CANBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T06_and_T26( neighb, t06, t26 );
      if ( *t06 != 1 || *t26 != 1 ) return( 0 );

      for ( n=0; n<27; n++ )
        if ( neighb[n] == _TOBECHANGED_ )
          neighb[n] = _BACKGROUND_VALUE_;

      Compute_T06_and_T26( neighb, &checkT06, &checkT26 );
      if ( checkT06 != 1 || checkT26 != 1 ) return( 0 );

      break;
  }



  return( 1 );
}










/**************************************************
 *
 * end point tests
 *
 **************************************************/


typedef int (*typeEndConditionSimplePoint)( int * );



static int _defaultEndConditionCurveSimplePoint( int *neighb __attribute__ ((unused)) )
{
  return( 0 );
}

static int _endConditionCurve2DSimplePoint( int *neighb )
{
  int i, n;
  for ( i=0, n=0; n<9; n++ ) {
    if ( neighb[n] == _CANBECHANGED_ ||
         neighb[n] == _ENDPOINT_ ||
         neighb[n] == _ANCHOR_ )
      i ++;
  }
  if ( i == 2 ) return( 1 );
  return( 0 );
}

static int _endConditionCurve3DSimplePoint( int *neighb )
{
  int i, n;
  for ( i=0, n=0; n<27; n++ ) {
    if ( neighb[n] == _CANBECHANGED_ ||
         neighb[n] == _ENDPOINT_ ||
         neighb[n] == _ANCHOR_ )
      i ++;
  }
  if ( i == 2 ) return( 1 );
  return( 0 );
}




typedef int (*typeEndConditionNonSimplePoint)( int *,
                                               int,
                                               int,
                                               enumTypeEndPoint );



static int _endCondition2DNonSimplePoint( int *neighb,
                                          int t04, int t08,
                                          enumTypeEndPoint typeEndPoint )
{
  int i, n;

  switch( typeEndPoint ) {
  default :
  case _SURFACE_ :
    if ( t04 >= 2 )
      return( 1 );
    break;
  case _PURE_SURFACE_ :
    if ( t04 == 2 )
      return( 1 );
    break;
  case _CURVE_ :
    if ( t08 >= 2 )
      return( 1 );
    break;
  case _PURE_CURVE_ :
    if ( t08 > 2 ) break;
    for ( i=0, n=0; n<9; n++ ) {
      if ( neighb[n] == _CANBECHANGED_ ||
           neighb[n] == _ENDPOINT_ ||
           neighb[n] == _ANCHOR_ )
        i ++;
    }
    if ( i == 3 )
      return( 1 );
    break;
  case _NO_END_POINT_ :
    break;
  }
  return( 0 );
}



static int _endCondition3DNonSimplePoint( int *neighb,
                                          int t06, int t26,
                                          enumTypeEndPoint typeEndPoint )
{
  int i, n;

  switch( typeEndPoint ) {
  default :
  case _SURFACE_ :
    if ( t06 >= 2 )
      return( 1 );
    break;
  case _PURE_SURFACE_ :
    if ( t06 == 2 )
      return( 1 );
    break;
  case _CURVE_ :
    if ( t26 >= 2 )
      return( 1 );
    break;
  case _PURE_CURVE_ :
    if ( t26 > 2 ) break;
    for ( i=0, n=0; n<27; n++ ) {
      if ( neighb[n] == _CANBECHANGED_ ||
           neighb[n] == _ENDPOINT_ ||
           neighb[n] == _ANCHOR_ )
        i ++;
    }
    if ( i == 3 )
      return( 1 );
    break;
  case _NO_END_POINT_ :
    break;
  }
  return( 0 );
}











/**************************************************
 *
 * thinning tools
 *
 **************************************************/


static int _getMaxValue( void *bufferIn, bufferType typeIn,
                         int *bufferDims,
                         int *valmin,
                         int *valmax )
{
    char *proc = "_getMaxValue";
    size_t v;
    size_t i;
    int min, max;

    v = bufferDims[0] * bufferDims[1] * bufferDims[2];
    switch( typeIn ) {
    default :
        if ( _verbose_ )
            fprintf( stderr, "%s: such image type not handled yet\n", proc );
        return( -1 );
    case UCHAR :
        {
            unsigned char *theBuf = (unsigned char*)bufferIn;
            for ( i=0, min=max=*theBuf; i<v; i++, theBuf++ ) {
                if ( max < *theBuf ) max = *theBuf;
                else if ( min > *theBuf ) min = *theBuf;
            }
        }
        /* end of case UCHAR
         */
        break;
    case SSHORT :
        {
            short int *theBuf = (short int*)bufferIn;
            for ( i=0, min=max=*theBuf; i<v; i++, theBuf++ ) {
                if ( max < *theBuf ) max = *theBuf;
                else if ( min > *theBuf ) min = *theBuf;
            }
         }
        /* end of case USHORT
         */
        break;
    case USHORT :
        {
            unsigned short int *theBuf = (unsigned short int*)bufferIn;
            for ( i=0, min=max=*theBuf; i<v; i++, theBuf++ ) {
                if ( max < *theBuf ) max = *theBuf;
                else if ( min > *theBuf ) min = *theBuf;
            }
        }
        /* end of case USHORT
         */
        break;
    }
    *valmin = min;
    *valmax = max;
    return( 1 );
}





static int _buildPointList( topologicalPointListList *pointListList,
                            unsigned char *theBuf,
                            void *theValues,
                            bufferType valuesType,
                            int *theDim )
{
    char *proc = "_buildPointList";
    int i, x, y, z, v;
    topologicalPoint point;
    typeIsInsideFunction _isInside = (typeIsInsideFunction)NULL;
    int npts = 0;

    if ( theDim[2] == 1 )
        _isInside = &_is2DInside;
    else
        _isInside = &_is3DInside;


    /* here the lists are allocated thanks to
     * addTopologicalPointToTopologicalPointList()
     * hence an over-allocation occurs
     *
     * one could first count the number of points to be
     * allocated per distance value, and then
     * make only the required allocation
     */

#define _BUILDPOINTLIST( TYPE ) {            \
  TYPE *theVal = (TYPE*)theValues;           \
  for ( npts=0, i=0, z=0; z<theDim[2]; z++ ) \
  for ( y=0; y<theDim[1]; y++ )              \
  for ( x=0; x<theDim[0]; x++, i++ ) {       \
    if ( theBuf[i] != _CANBECHANGED_ )       \
      continue;                              \
    initTopologicalPoint( &point );          \
    point.x = x;                             \
    point.y = y;                             \
    point.z = z;                             \
    point.i = i;                             \
    point.iteration = 0;                     \
    point.value = theVal[i];                 \
    point.type = _CANBECHANGED_;             \
    point.isinside = (*_isInside)( point.x, point.y, point.z, theDim ); \
    v = theVal[i] - pointListList->firstindex; \
    if ( pointListList->binlength > 1 )      \
      v /= pointListList->binlength;         \
    if ( v < 0 || v >= pointListList->n_data ) { \
      if ( _verbose_ )                       \
        fprintf( stderr, "%s: unable to add point, list overflow\n", proc );   \
      return( -1 );                          \
    }                                        \
    if ( addTopologicalPointToTopologicalPointList( &(pointListList->data[v]), \
                                                    &point ) != 1 ) {   \
      if ( _verbose_ )                       \
        fprintf( stderr, "%s: unable to add point to list\n", proc );   \
      return( -1 );                          \
    }                                        \
    npts ++;                                 \
  }                                          \
}

    switch( valuesType ) {
    default :
        if ( _verbose_ )
            fprintf( stderr, "%s: such image type not handled yet\n", proc );
        return( -1 );
    case UCHAR :
        _BUILDPOINTLIST( u8 );
        break;
    case SSHORT :
        _BUILDPOINTLIST( s16 );
        break;
    case USHORT :
        _BUILDPOINTLIST( u16 );
        break;
    }

    return( npts );
}










/**************************************************
 *
 * thinning procedures
 *
 **************************************************/



int thinningThreshold( void *bufferIn, bufferType typeIn,
                       unsigned char *bufferOut, int *bufferDims,
                       float lowThreshold,
                       float highThreshold )
{
    char *proc = "thinningThreshold";
    size_t v;
    size_t i;
    unsigned char *resBuf = bufferOut;
    int iLowThreshold;
    int iHighThreshold;

    iLowThreshold = ( lowThreshold >= 0.0 ) ? (int)(lowThreshold+0.5) : (int)(lowThreshold-0.5);
    iHighThreshold = ( highThreshold >= 0.0 ) ? (int)(highThreshold+0.5) : (int)(highThreshold-0.5);

    v = bufferDims[0] * bufferDims[1] * bufferDims[2];
    switch( typeIn ) {
    default :
        if ( _verbose_ )
            fprintf( stderr, "%s: such image type not handled yet\n", proc );
        return( -1 );
    case UCHAR :
        {
            unsigned char *theBuf = (unsigned char*)bufferIn;
            if ( iLowThreshold < iHighThreshold ) {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else if ( *theBuf < iHighThreshold ) {
                        *resBuf = _CANBECHANGED_VALUE_;
                    } else {
                        *resBuf = _ANCHOR_VALUE_;
                    }
                }
            }
            else {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else {
                        *resBuf = _CANBECHANGED_VALUE_;
                    }
                }
            }
        }
        /* end of case UCHAR
         */
        break;
    case SSHORT :
        {
            short int *theBuf = (short int*)bufferIn;
            if ( iLowThreshold < iHighThreshold ) {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else if ( *theBuf < iHighThreshold ) {
                        *resBuf = _CANBECHANGED_VALUE_;
                    } else {
                        *resBuf = _ANCHOR_VALUE_;
                    }
                }
            }
            else {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else {
                        *resBuf = _CANBECHANGED_VALUE_;
                    }
                }
            }
        }
        /* end of case USHORT
         */
        break;
    case USHORT :
        {
            unsigned short int *theBuf = (unsigned short int*)bufferIn;
            if ( iLowThreshold < iHighThreshold ) {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else if ( *theBuf < iHighThreshold ) {
                        *resBuf = _CANBECHANGED_VALUE_;
                    } else {
                        *resBuf = _ANCHOR_VALUE_;
                    }
                }
            }
            else {
                for ( i=0; i<v; i++, theBuf++, resBuf++ ) {
                    if ( *theBuf < iLowThreshold ) {
                        *resBuf = _BACKGROUND_VALUE_;
                    } else {
                        *resBuf = _CANBECHANGED_VALUE_;
                    }
                }
            }
        }
        /* end of case USHORT
         */
        break;
    }
    return( 1 );
}










/**************************************************
 *
 * thinning procedures
 *
 **************************************************/



/* used in vt-exec/axeMedianOfComponents.c
 */
/* theBuf is assumed to be pre-processed/threholded
 * it contains 3 values :
 * - _BACKGROUND_VALUE_
 * - _CANBECHANGED_VALUE_
 * - _ANCHOR_VALUE_
 */

int chamferBasedThinning( unsigned char *theBuf,
                          int *theDim,
                          const int chamfer,
                          typeThinningParameters *p )
{
  char *proc = "chamferBasedThinning";
  unsigned short int *theDistance = (unsigned short int *)NULL;
  int v = theDim[0]*theDim[1]*theDim[2];
  int i;
  float exectime[4];
  int tmeas=0;
  int ret;

  exectime[tmeas++] = (float)clock();

  /* distance map computation
     1. initialisation
     2. computation
  */
  theDistance = (unsigned short *)vtmalloc( v * sizeof(unsigned short), "theDistance", proc );
  if ( theDistance == (unsigned short *)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate distance map\n", proc );
    return( -1 );
  }

  for ( i=0; i<v; i++ ) {
    theDistance[i] = (unsigned short)( ( theBuf[i] > 0 ) ? 0 : 65535 );
  }
  

  if ( _verbose_ >= 1 ) 
    fprintf( stderr, " ... distance computation" );



  switch( chamfer ) {
  default :
  case 3 :
    if ( Compute3DNormalizedChamfer3x3x3( (void*)theDistance, USHORT,
                                (void*)theDistance, USHORT,
                                theDim ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to compute distance map\n", proc );
      vtfree( theDistance );
      return( -1 );
    }
    break;
  case 5 :
    if ( Compute3DNormalizedChamfer5x5x5( (void*)theDistance, USHORT,
                                (void*)theDistance, USHORT,
                                theDim ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to compute distance map\n", proc );
      vtfree( theDistance );
      return( -1 );
    }
    break;
  }



  if ( _verbose_ >= 1 ) 
    fprintf( stderr, "d\n" );




  exectime[tmeas++] = (float)clock();
  if ( _time_ )
    _PrintTimeFrom2Clocks( "total   passed time", exectime[tmeas-2], exectime[tmeas-1] );
 
  ret = valueBasedThinning( theBuf, (void*)theDistance, USHORT, theDim, p );

  vtfree( theDistance );

  return( ret );
}









/* theBuf is assumed to be pre-processed/threholded
 * it contains 3 values :
 * - _BACKGROUND_VALUE_      0
 * - _CANBECHANGED_VALUE_  150
 * - _ANCHOR_VALUE_        255
 * only point with _CANBECHANGED_VALUE_ will be processed
 * theValues will be used to sort points and to process
 * them according to the values, typically it can be
 * a distance map fro thinning (and the points are
 * processed by increasing order)
 */
int valueBasedThinning( unsigned char *theBuf,
                        void *theValues,
                        bufferType valuesType,
                        int *theDim,
                        typeThinningParameters *par )
{
  char *proc = "valueBasedThinning";

  int minPossibleValue, maxPossibleValue;
  int firstPossibleIndex, lastPossibleIndex;
  int incIndexList;
  int index;

  topologicalPointListList pointListList;
  int nbPoints;

  typeNeighborhood offsetForSimplicity;
  typeCheckThickness _checkThickness = (typeCheckThickness)NULL;
  typeIsPointSimple _isPointSimple = (typeIsPointSimple)NULL;
  typeEndConditionSimplePoint _endConditionSimplePoint = &_defaultEndConditionCurveSimplePoint;
  typeEndConditionNonSimplePoint _endConditionNonSimplePoint = (typeEndConditionNonSimplePoint)NULL;

  int ndirection = 0;

  int successfuliteration, iteration;

  int cycle;
  int totDelPts = 0;
  int totDelPtsCycle, totEndPtsCycle;
  int nbDelPtsCycle, nbEndPtsCycle;
  int direction;
  int nbDelPtsDirection, nbEndPtsDirection;

  int nbMarkedPts, nbPtsValue;
  int p, i;
  topologicalPointList *valueList;
  /*
  topologicalPoint *arrayPoint;
  */
  topologicalPoint *ptrPoint, tmpPoint;

  int neighb[27];
  int tback, tfore;

  int v = theDim[0]*theDim[1]*theDim[2];

  float exectime[4];
  int tmeas = 0;





  /*--------------------------------------------------
   *
   * start
   *
   --------------------------------------------------*/


  exectime[tmeas++] = (float)clock();

  
  
  

  /*--------------------------------------------------
   *
   * build point list
   *
   --------------------------------------------------*/

  if ( _debug_ || _verbose_ >= 2 )
    fprintf( stderr, " ... list construction" );

  initTopologicalPointListList( &pointListList );
  if ( par->binLength > 1 )
    pointListList.binlength = par->binLength;

  if ( _getMaxValue( theValues, valuesType, theDim,
                     &minPossibleValue, &maxPossibleValue ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute min and max values\n", proc );
      return( -1 );
  }

  if ( _debug_ || _verbose_ >= 2 ) {
      fprintf( stderr, ", values in [%d %d]", minPossibleValue, maxPossibleValue );
  }

  if ( allocTopologicalPointListList( &pointListList, 0, maxPossibleValue ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate point list array\n", proc );
      return( -1 );
  }

  /* put all the '_CANBECHANGED_VALUE_' points
   * in the list, according to their value
   */
  nbPoints = _buildPointList( &pointListList, theBuf,
                              theValues, valuesType, theDim );

  if ( nbPoints < 0 ) {
      freeTopologicalPointListList( &pointListList );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when building list of points\n", proc );
      return( -1 );
  }
  else if ( nbPoints == 0 ) {
      freeTopologicalPointListList( &pointListList );
      if ( _debug_ || _verbose_ >= 2 )
        fprintf( stderr, "\n" );
      if ( _verbose_ )
        fprintf( stderr, "%s: no points to be thinned\n", proc );
      return( 0 );
  }

  /* recall that
  * - pointListList.firstindex = 0
  * - pointListList.lastindex = maxPossibleValue
  */
  switch( par->typeOrdering ) {
  default :
      freeTopologicalPointListList( &pointListList );
      if ( _verbose_ )
        fprintf( stderr, "%s: such ordering type not defined yet\n", proc );
      return( -1 );
  case _INCREASING_ORDER_ :
      firstPossibleIndex = 0;
      lastPossibleIndex = pointListList.n_data-1;
      incIndexList = 1;
      break;
  case _DECREASING_ORDER_ :
      firstPossibleIndex = pointListList.n_data-1;
      lastPossibleIndex = 0;
      incIndexList = -1;
      break;
  }

  if ( _debug_ || _verbose_ >= 2 )
    fprintf( stderr, "\n" );

  exectime[tmeas++] = (float)clock();
  if ( _time_ ) {
    _PrintTimeFrom2Clocks( "partial passed time", exectime[tmeas-2], exectime[tmeas-1] );
    _PrintTimeFrom2Clocks( "total   passed time", exectime[0], exectime[tmeas-1] );
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... there are %d points, values in [%d %d] indexes in [%d %d]\n",
             nbPoints, minPossibleValue, maxPossibleValue,
             firstPossibleIndex, lastPossibleIndex );
  }


  
  
  

  /*--------------------------------------------------
   *
   * list has been built,
   * set tests, offsets, etc depending on the image dimension
   * pre-compute offsets for point access speed-up
   *
   --------------------------------------------------*/



  if ( theDim[2] == 1 ) {
    defineNeighborsForSimplicity( &offsetForSimplicity, theDim, 8 );
    _checkThickness = &_check2DThickness;
    _isPointSimple = &_is2DPointSimple;
    _endConditionNonSimplePoint = &_endCondition2DNonSimplePoint;
    ndirection = 4;
    switch ( par->typeEndPoint ) {
    default :
      break;
    case _SURFACE_ :
    case _CURVE_ :
      _endConditionSimplePoint = &_endConditionCurve2DSimplePoint;
      break;
    }
  }
  else {
    defineNeighborsForSimplicity( &offsetForSimplicity, theDim, 26 );
    _checkThickness = &_check3DThickness;
    _isPointSimple = &_is3DPointSimple;
    _endConditionNonSimplePoint = &_endCondition3DNonSimplePoint;
    ndirection = 6;
    switch ( par->typeEndPoint ) {
    default :
      break;
    case _CURVE_ :
      _endConditionSimplePoint = &_endConditionCurve3DSimplePoint;
      break;
    }
  }

  
  


  /*--------------------------------------------------
   *
   * write some information
   *
   --------------------------------------------------*/
 
  if ( _verbose_ >= 2 ) {
    _fprintfTypeThinningParameters( stderr, par );
  }





  /*--------------------------------------------------
   *
   * thinning
   *
   --------------------------------------------------*/


  /* general loop from first possible distance to last possible distance 
  */
  index = firstPossibleIndex;
  successfuliteration = iteration = 0;
  totDelPts = 0;
  do {

    /************************************************************
     * selection of a non-empty list
     ************************************************************/

    /* advance to first non-empty list
     * while updating firstPossibleIndex
     * case index == firstPossibleIndex, occurs
     * - at the beginning
     * - when points have been changed
     */
    if ( index == firstPossibleIndex ) {
      while ( pointListList.data[firstPossibleIndex].n_data == 0 ) {
        index += incIndexList;
        firstPossibleIndex += incIndexList;
        if ( (incIndexList > 0 && firstPossibleIndex > lastPossibleIndex)
             || (incIndexList < 0 && firstPossibleIndex < lastPossibleIndex) )
        break;
      }
    }
    if ( (incIndexList > 0 && firstPossibleIndex > lastPossibleIndex)
         || (incIndexList < 0 && firstPossibleIndex < lastPossibleIndex) )
    break;
    
    /* advance to next non-empty list
     * case index <> firstPossibleIndex, occurs
     * when no points have been changed at previous 'index'
     */
    while ( pointListList.data[index].n_data == 0 ) {
      index += incIndexList;
      if ( (incIndexList > 0 && index > lastPossibleIndex)
           || (incIndexList < 0 && index < lastPossibleIndex) )
      break;
    }
    if ( (incIndexList > 0 && index > lastPossibleIndex)
         || (incIndexList < 0 && index < lastPossibleIndex) )
    break;

    valueList = &(pointListList.data[index]);

    /************************************************************
     * end of selection of a non-empty list
     ************************************************************/



    /* loop on cycle for one given index
       cycle = a complete loop on all the directions
       iteration =  cycle
    */
    
    cycle = 0;
    totDelPtsCycle = totEndPtsCycle = 0;

    do { 

      /* loop on directions
       */
      nbDelPtsCycle = nbEndPtsCycle = 0;

      for ( direction = 0; direction < ndirection; direction ++ ) {

        /* loop on points for one direction
         */
        nbDelPtsDirection = nbEndPtsDirection = 0;
        nbMarkedPts = 0;
        nbPtsValue = valueList->n_data;

        for ( p = 0; p < valueList->n_data; p ++ ) {

           ptrPoint = &(valueList->data[p]);

          /* get the neighborhood
           */
          extractNeighborhood( neighb, ptrPoint,
                               theBuf, theDim,
                               &offsetForSimplicity );

          /* check the thickness for the direction
           */
          if ( (*_checkThickness)( neighb, par->typeThickness, par->typeChanging, direction ) == 0 )
            continue;

          /* test the simplicity
             if simple and not endpoint, it is marked for deletion
           */
          if ( (*_isPointSimple)( neighb, &tback, &tfore, par->typeChanging ) == 1 ) {
            /* if there is an end condition for simple points,
               it is to be tested here
            */
            if ( (*_endConditionSimplePoint)( neighb ) == 0 ) {
              ptrPoint->type = theBuf[ ptrPoint->i ] = _TOBECHANGED_;
              nbMarkedPts ++;
            }
            continue;
          }

           /* shrinking: no further test
           */
          if ( par->typeEndPoint == _NO_END_POINT_ )
            continue;

          /* end point condition for non-simple points
             - check whether it is deep enough
               1. the distance/index is large enough
               2. the number of cycles is large enough
                  recall that a cycle is the spanning of all directions
                  for a given 'index', thus a plateau of constant 'index'
                  has been eroded by this number of cycle before
                  end points may appears
              - if yes, it can be considered as a end point candidate
             else it will be deleted
           */
          if ( par->valueBeforeEnding > 0
               && ( (incIndexList > 0 && index < par->valueBeforeEnding)
                    || (incIndexList < 0 && index > par->valueBeforeEnding) ) )
            continue;
          if ( par->cyclesBeforeEnding > 0 && cycle < par->cyclesBeforeEnding )
            continue;

          /* check whether the non-simple point satisfies an end point condition
           */
          if ( (*_endConditionNonSimplePoint)( neighb, tback, tfore, par->typeEndPoint ) == 1 ) {
            ptrPoint->type = theBuf[ ptrPoint->i ] = _ENDPOINT_;
            nbMarkedPts ++;
          }
        }
        /* end of loop on points for one direction
         */

        /* update points and image for one direction
         */
        if ( nbMarkedPts > 0 ) {
          for ( p = 0; p < valueList->n_data; p ++ ) {
            ptrPoint = &(valueList->data[p]);
            switch( ptrPoint->type ) {
            default :
            case _BACKGROUND_ :
            case _ANCHOR_ :
              fprintf( stderr, "%s: WARNING this case should not occur\n", proc );
              break;
            case _CANBECHANGED_ :
              break;
            case _TOBECHANGED_ :
              switch( par->typeChanging ) {
              default :
              case _FOREGROUND_TO_BACKGROUND_ :
                  ptrPoint->type = theBuf[ ptrPoint->i ] = _BACKGROUND_VALUE_;
                  break;
              case _BACKGROUND_TO_FOREGROUND_ :
                  ptrPoint->type = theBuf[ ptrPoint->i ] = _FOREGROUND_VALUE_;
                  break;
              }
              nbDelPtsDirection ++;
            case _ENDPOINT_ :
              if ( ptrPoint->type == _ENDPOINT_ ) {
                  switch( par->typeChanging ) {
                  default :
                      break;
                  case _FOREGROUND_TO_BACKGROUND_ :
                      ptrPoint->type = theBuf[ ptrPoint->i ] = _FOREGROUND_VALUE_;
                      break;
                  }
                  nbEndPtsDirection ++;
              }
              tmpPoint = valueList->data[ valueList->n_data-1 ];
              valueList->data[ valueList->n_data-1 ] = *ptrPoint;
              *ptrPoint = tmpPoint;
              valueList->n_data --;
              p --;
              break;
            }
          }
          nbDelPtsCycle += nbDelPtsDirection;
          nbEndPtsCycle += nbEndPtsDirection;
          totDelPts += nbDelPtsDirection;
        }
        /* end of update points and image for one direction
         */

        if ( _verbose_ >= 2 ) {
          fprintf( stderr, " #%8d", iteration );
          fprintf( stderr, " Dir=%d/%d Cyc.=%2d[%d] Dist.=%3d[%d]/[%d %d]",
                   direction, ndirection,
                   cycle, par->cyclesBeforeEnding,
                   index, par->valueBeforeEnding,
                   firstPossibleIndex, lastPossibleIndex );
          fprintf( stderr, " Pts=%8d/%8d Del=%6d End=%5d",
                   valueList->n_data, nbPtsValue,
                   nbDelPtsDirection, nbEndPtsDirection );
          fprintf( stderr, " Chg=%8d", totDelPts );
          if ( _verbose_ >= 3 ) fprintf( stderr, "\n" );
          else                  fprintf( stderr, "\r" );
        }

      }
      /* end of loop on direction 
       * = end of cycle
       */

      /* go to next cycle (for the same given index) if something happens
       */
      if ( nbDelPtsCycle > 0 || nbEndPtsCycle > 0 ) {
        totDelPtsCycle += nbDelPtsCycle;
        totEndPtsCycle += nbEndPtsCycle;
        cycle ++;
        successfuliteration ++;
      }
      iteration++;
      
    } while ( (nbDelPtsCycle > 0 || nbEndPtsCycle > 0)
              && (par->maxIteration <= 0 || par->maxIteration > iteration) );

    /* processing has been ended for a given index
     * there have been no deleted simple points and no
     * detected endpoints in the last cycle
     * if some points have been deleted, go back to minDistance since points with
       smallest distances can be simple now, else increment distance
    */
    if ( totDelPtsCycle > 0 &&
         ( (incIndexList > 0 && index > firstPossibleIndex)
           || (incIndexList < 0 && index < firstPossibleIndex) ) ) {
      index = firstPossibleIndex;
    }
    else {
      index += incIndexList;
    }

  } while ( ( (incIndexList > 0 && index <= lastPossibleIndex)
              || (incIndexList < 0 && index >= lastPossibleIndex) )
            && (par->maxIteration <= 0 || par->maxIteration > iteration) );
  /* fin de la boucle generale
   */




  if ( _verbose_ == 2 ) fprintf( stderr, "\n" );


  exectime[tmeas++] = (float)clock();
  if ( _time_ ) {
    _PrintTimeFrom2Clocks( "partial passed time", exectime[tmeas-2], exectime[tmeas-1] );
    _PrintTimeFrom2Clocks( "total   passed time", exectime[0], exectime[tmeas-1] );
  }





  /* release memory
  */
  freeTopologicalPointListList( &pointListList );



  /* binarise result image
   */
  switch( par->typeChanging ) {
  default :
  case _FOREGROUND_TO_BACKGROUND_ :
      for ( i=0; i<v; i++ ) {
        if ( theBuf[i] > 0 ) theBuf[i] = 255;
      }
      break;
  case _BACKGROUND_TO_FOREGROUND_ :
      for ( i=0; i<v; i++ ) {
        if ( theBuf[i] < 255 ) theBuf[i] = 0;
      }
      break;
  }

  return( 1 );
}










/**************************************************
 *
 * time related procedure
 *
 **************************************************/



static void _ComputeTimeFrom2Clocks( float c1, float c2,
                                     int *hours, int *minutes, float *seconds )
{
  double d = ( (double)c2 / (double)CLOCKS_PER_SEC ) - 
    ( (double)c1 / (double)CLOCKS_PER_SEC );
  *hours = *minutes = 0;
  *seconds = 0.0;

  if ( d > 3600 ) {
    *hours = (int)(d / 3600);
    d -= *hours * 3600.0;
  }
  if ( d > 60 ) {
    *minutes = (int)(d / 60);
    d -= *minutes * 60.0;
  }
  *seconds = d;
}



static void _PrintTimeFrom2Clocks( char *str, float c1, float c2 )
{
  int h, m;
  float s;
  int l, b=10;
  char format[10] = "%s";
  
  _ComputeTimeFrom2Clocks( c1, c2, &h, &m, &s );

  if ( str != (char*)NULL ) {
    l = strlen(str);
    if ( l % b == 0 ) sprintf(format, "%%%ds", l / b );
    else              sprintf(format, "%%%ds", 1+(l / b) );
    fprintf(stderr, format, str);
  } else {
    fprintf(stderr, "          " );
  }
  if ( h > 0 ) {
    fprintf(stderr, " %2d h", h);
  } else {
    fprintf(stderr, "     ");
  }
  if ( m > 0 ) {
    fprintf(stderr, " %2d mn", m);
  } else {
    fprintf(stderr, "      ");
  }
  fprintf(stderr, " %9.6f s\n", s);
}
