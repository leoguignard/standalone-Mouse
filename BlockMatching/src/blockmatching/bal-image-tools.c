/*************************************************************************
 * bal-image-tools.c -
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

#include <bal-image-tools.h>

static int _verbose_ = 1;



void BAL_SetVerboseInBalImageTools( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalImageTools(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalImageTools(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}




/*--------------------------------------------------
 *
 * IMAGE GRID
 *
 --------------------------------------------------*/

int BAL_DrawGrid( bal_image *theIm,
                  bal_integerPoint *offset,
                  bal_integerPoint *spacing,
                  float value )
{
  char *proc = "BAL_DrawGrid";
  size_t i, j, k;
  int ivalue;

  ivalue = ( value > 0 ) ? (int)(value + 0.5) : (int)(value - 0.5);

  switch( theIm->type ) {

  default : 
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
    
  case UCHAR :
    {
      unsigned char ***theBuf = (unsigned char ***)(theIm->array);
      if ( ivalue < 0 )
        ivalue = 0;
      else if ( ivalue > 255 )
        ivalue = 255;
      if ( 0 ) {
        for ( k=0; k<theIm->nplanes; k++ )
        for ( j=0; j<theIm->nrows; j++ )
        for ( i=0; i<theIm->ncols; i++ )
          theBuf[k][j][i] = 0;
      }
      if ( offset->x >= 0 && spacing->x >= 1 && theIm->ncols > 1 ) {
        for ( k=0; k<theIm->nplanes; k++ )
        for ( j=0; j<theIm->nrows; j++ )
        for ( i=offset->x; i<theIm->ncols; i+=spacing->x)
          theBuf[k][j][i] = ivalue;
      }

      if ( offset->y >= 0 && spacing->y >= 1 && theIm->nrows > 1 ) {
        for ( k=0; k<theIm->nplanes; k++ )
        for ( j=offset->y; j<theIm->nrows; j+=spacing->y )
        for ( i=0; i<theIm->ncols; i++ )
          theBuf[k][j][i] = ivalue;
      }
 
      if ( offset->z >= 0 && spacing->z >= 1 && theIm->nplanes > 1 ) {
        for ( k=offset->z; k<theIm->nplanes; k+=spacing->z )
        for ( j=0; j<theIm->nrows; j++ )
        for ( i=0; i<theIm->ncols; i++ )
          theBuf[k][j][i] = ivalue;
      }

    }
    break;

  }

  return( 1 );

}





/*--------------------------------------------------
 *
 * IMAGE Mosaic
 *
 --------------------------------------------------*/

int BAL_DrawMosaic( bal_image *theIm,
                  bal_integerPoint *spacing )
{
  char *proc = "BAL_DrawMosaic";
  size_t i, j, k;
  int ix, iy, iz, ixy;
  int dx=1, dy=1 ,dz=1;
  int ivalue;

  if ( spacing->x >= 1 ) {
    dx = theIm->ncols / spacing->x;
    if ( theIm->ncols % spacing->x > 1 ) dx ++;
  }
  if ( spacing->y >= 1 ) {
    dy = theIm->nrows / spacing->y;
    if ( theIm->nrows % spacing->y > 1 ) dy ++;
  }
  if ( spacing->z >= 1 ) {
    dz = theIm->nplanes / spacing->z;
    if ( theIm->nplanes % spacing->z > 1 ) dz ++;
  }

  if ( 0 ) {
    fprintf( stderr, "%s: block dimensions = %d x %d x %d\n", proc, dx, dy, dz );
  }

  switch( theIm->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
    {
      unsigned char ***theBuf = (unsigned char ***)(theIm->array);
      if ( theIm->nplanes > 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: 3D case not handled yet\n", proc );
        return( -1 );
      }
      else {
        for ( k=0; k<theIm->nplanes; k++ )
        for ( j=0; j<theIm->nrows; j++ )
        for ( i=0; i<theIm->ncols; i++) {
          ix = ( spacing->x >= 1 ) ? i / spacing->x : 0;
          iy = ( spacing->y >= 1 ) ? j / spacing->y : 0;
          iz = ( spacing->z >= 1 ) ? k / spacing->z : 0;
          ixy = ix + iy;
          /* diagonale x+k=cte,
           * les labels utilises sont ceux necessaires pour les lignes x+y=k, avec k<cte
           * - pour une diagonale cte=k, il y a k+1 points
           * - les labels deja utilises sont sum(i=0...cte-1) (i+1)
           *     = sum(i=1...cte) (i) = cte*(cte+1)/2
           * - on ajoute 1 pour commencer la nouvelle diagonale (ix = 0)
           * - on ajoute iy
           *
           * labels a enlever selon x
           * si cte >= dx (N = cte-dx >= 0) il faut enlever
           *   sum(k=0...N) (k+1) = sum(k=1...N+1) (k) = (N+1)*(N+2)/2
           * attention: on compte la diagonale en cours
           *
           * labels a enlever selon y
           * meme calcul, meme il ne faut consider que les diagonales avant la
           * diagonale en cours
           */

          ivalue = ixy*(ixy+1)/2 + 1 + iy;
          if ( ixy >= dx ) ivalue -= (ixy-dx+1)*(ixy-dx+2)/2;
          if ( ixy >= dy ) ivalue -= (ixy-dy)*(ixy-dy+1)/2;

          while ( ivalue >= 256 ) ivalue -= 256;
          theBuf[k][j][i] = ivalue;
        }
      }

    }
    break;

  }

  return( 1 );

}
