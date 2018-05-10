/*************************************************************************
 * cspline.c - Cubic splines
 *
 * $Id: cspline.c,v 1.3 2006/04/14 08:38:55 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Alexis Roche
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Oct 11 23:10:01 MEST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vtmalloc.h>
#include <chunks.h>
#include <convert.h>

#include <cspline4x4.h>
#include <cspline.h>

static int _verbose_ = 1;

void setVerboseInCspline( int v )
{
  _verbose_ = v;
}

void incrementVerboseInCspline( )
{
  _verbose_ ++;
}

void decrementVerboseInCspline( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}








/**************************************************
 *
 * Fonctions pour l'interpolation spline cubique.
 * A. Roche
 *
 ***************************************************/


/* Cette fonction renvoie la valeur de la fonction B-spline
   cubique elementaire en x */

static double CubicBspline ( double x )
{
  double aux;
  
  if ( x >= 0.0 ) {
    if ( x >=  2.0 ) return( 0.0 );
    if ( x < 1 ) {
      aux = x*x;
      return( (double)2.0/(double)3.0 - aux + (double)0.5*aux*x );
    }
    aux = 2.0 - x;
    return( aux*aux*aux/(double)6.0 );
  } 
  
  /* ici x < 0
   */
  if ( x <= -2.0 ) return( 0.0 );
  if ( x > -1 ) {
    aux = x*x;
    return( (double)2.0/(double)3.0 - aux - (double)0.5*aux*x );
  }
  aux = 2.0 + x;
  return( aux*aux*aux/(double)6.0 );
}





static double CubicBspline_FirstDeriv ( double x ) 
{
  if ( x == 0.0 ) return( 0.0 );

  if ( x > 0.0 ) {
    if ( x >=  2.0 ) return( 0.0 );
    if ( x < 1 ) {
      return( (-2.0 + 1.5*x)*x );
    } 
    return( -0.5 * (2.0-x) * (2.0-x) );
  }

  /* ici x < 0
   */
  if ( x <= -2.0 ) return( 0.0 );
  if ( x > -1 ) {  
    return( (-2.0 - 1.5*x)*x );
  }
  return( 0.5 * (2.0+x) * (2.0+x) );
}




static double CubicBspline_SecndDeriv ( double x ) 
{
  if ( x > 0.0 ) {
    if ( x >=  2.0 ) return( 0.0 );
    if ( x < 1 ) {
      return( -2.0 + 3.0*x );
    }
    return( 2.0 - x );
  }

  /* ici x < 0
   */
  if ( x <= -2.0 ) return( 0.0 );
  if ( x > -1 ) {  
    return( -2.0 - 3.0*x );
  }
  return( 2.0 + x );
}























/**************************************************
 *
 * Fonctions pour l'interpolation spline cubique.
 * A. Roche
 *
 ***************************************************/







int Reech3DCSpline4x4WithCoefficients( typeCSplineCoefficients *theCoeff,
                       r32* resBuf,  /* result buffer */
                       int *resDim,  /* dimensions of this buffer */
                       double *mat,
                       int slice, int *derivative )
{
  char *proc = "Reech3DCSpline4x4";
  register int i, j, k=slice, ix, iy, iz;
  register double x, y, z;
  register double dx, dy, dz;
  register double res;
  
  int rdimx=resDim[0];
  int rdimy=resDim[1];

  register r32* rbuf = resBuf;
  register float *tbuf;
  float *bufCoeff = theCoeff->theCoeff;

  int tdimx  = theCoeff->theDim[0];
  int tdimxy = theCoeff->theDim[0] * theCoeff->theDim[1];

  int tdimx3 = theCoeff->theDim[0]-3;
  int tdimy3 = theCoeff->theDim[1]-3;
  int tdimz3 = theCoeff->theDim[2]-3;
 
  register double ddimx = (double)theCoeff->theDim[0]-2.0-0.5;
  register double ddimy = (double)theCoeff->theDim[1]-2.0-0.5;
  register double ddimz = (double)theCoeff->theDim[2]-2.0-0.5;

  double cx[4], cy[4], cz[4], cyz;
  int offset[64];

  int l, m, n, o;

  double (*bsplinx)(double);
  double (*bspliny)(double);
  double (*bsplinz)(double);


  if ( derivative == (int*) NULL ) {
    bsplinx = & CubicBspline;
    bspliny = & CubicBspline;
    bsplinz = & CubicBspline;
  } else {
    switch ( derivative[0] ) {
    default :
      if ( _verbose_ ) 
        fprintf( stderr, "%s: invalid derivation order in x\n", proc );
      return( -1 );
    case 0 : bsplinx = & CubicBspline;            break;
    case 1 : bsplinx = & CubicBspline_FirstDeriv; break;
    case 2 : bsplinx = & CubicBspline_SecndDeriv; break;
    }
    
    switch ( derivative[1] ) {
    default :
      if ( _verbose_ ) 
        fprintf( stderr, "%s: invalid derivation order in y\n", proc );
      return( -1 );
    case 0 : bspliny = & CubicBspline;            break;
    case 1 : bspliny = & CubicBspline_FirstDeriv; break;
    case 2 : bspliny = & CubicBspline_SecndDeriv; break;
    }
    
    switch ( derivative[2] ) {
    default :
      if ( _verbose_ ) 
        fprintf( stderr, "%s: invalid derivation order in z\n", proc );
      return( -1 );
    case 0 : bsplinz = & CubicBspline;            break;
    case 1 : bsplinz = & CubicBspline_FirstDeriv; break;
    case 2 : bsplinz = & CubicBspline_SecndDeriv; break;
    }
  }
  
  

  for ( o = 0, l = 0; l < 4; l ++ )
  for ( m = 0; m < 4; m ++ )
  for ( n = 0; n < 4; n ++, o ++ ) {
    offset[o] = l * tdimxy + m * tdimx + n;
  }
  
  

  for ( j = 0; j < rdimy; j ++ )
  for ( i = 0; i < rdimx; i ++, rbuf ++ ) {
    /* computation of the corresponding point after deformation 
       
       (x,y,x) -> (ix=(int)x,iy=(int)x,iz=(int)x)
       les coefficients de la spline a prendre en compte
       sont [1+ix+{-1,0,1,2}, 1+iy+{-1,0,1,2}, 1+iz+{-1,0,1,2}]
       
       On ajoute 1 aux coordonnees, car on a ajoute
       une bordure aux coefficients de la spline.
     */
    
    x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
    /*
    if ((x < 0) || ( x > tdimx3)) { *rbuf = 0; continue; }
    */
    if ((x < -0.5) || ( x > ddimx)) { *rbuf = 0; continue; }

    y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
    /*
    if ((y < 0) || ( y > tdimy3)) { *rbuf = 0; continue; }
    */
    if ((y < -0.5) || ( y > ddimy)) { *rbuf = 0; continue; }

    z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
    /*
    if ((z < 0) || ( z > tdimz3)) { *rbuf = 0; continue; }
    */
    if ((z < -0.5) || ( z > ddimz)) { *rbuf = 0; continue; }
    
    
    /* ici on a 0 <   x   < theCoeff->theDim[0]-3
                0 <= ix   < theCoeff->theDim[0]-3
                1 <= ix+1 < theCoeff->theDim[0]-2
                0 <= ix+1-1 && ix+1+2 < theCoeff->theDim[0]

       maintenant on a -0.5 <         x   <        theCoeff->theDim[0]-2.5
                          0 <= ix   <= theCoeff->theDim[0]-3
                          1 <= ix+1 <= theCoeff->theDim[0]-2

       les extrema pour les coefficients verifient
                          0 <= ix+1-1 <= theCoeff->theDim[0]-3
                          3 <= ix+1+2 <= theCoeff->theDim[0]
       dans ce dernier cas, on suppose que ca va etre 0 ...
    */
    ix = (int)x;    dx = x-ix;
    iy = (int)y;    dy = y-iy;
    iz = (int)z;    dz = z-iz;



    if ( x < 0.0 ) { x = 0.0; ix = 0; dx = 0.0; }
    if ( y < 0.0 ) { y = 0.0; iy = 0; dy = 0.0; }
    if ( z < 0.0 ) { z = 0.0; iz = 0; dz = 0.0; }

    /* Pre-calcul des valeurs de bsplines pour chaque dimension:
       
       (*bsplinx)(dx-1), (*bsplinx)(dx), (*bsplinx)(dx+1), (*bsplinx)(dx+2)  
       (*bspliny)(dy-1), (*bspliny)(dy), (*bspliny)(dy+1), (*bspliny)(dy+2)  
       (*bsplinz)(dz-1), (*bsplinz)(dz), (*bsplinz)(dz+1), (*bsplinz)(dz+2) 
       
       A. Roche
       
       Dans l'image des coefficients,
       (*bsplinx)(dx-1) correspond a ix, 
       (*bsplinx)(dy-1) correspond a iy, 
       (*bsplinx)(dz-1) correspond a iz, 
       
    */
    for ( l = 0; l < 4; l ++ ) {
      cx[l] = (*bsplinx)( dx - (double)(l-1) );
      cy[l] = (*bspliny)( dy - (double)(l-1) );
      cz[l] = (*bsplinz)( dz - (double)(l-1) );
    }

    tbuf = bufCoeff;
    tbuf += iz * tdimxy + iy * tdimx + ix;



    /* 
     * are we on the border or not ? 
     */
    if ( (ix < tdimx3) && (iy < tdimy3) && (iz < tdimz3) ) {
      

      for ( res = 0.0, o = 0, l = 0; l < 4; l ++ )
      for ( m = 0; m < 4; m ++ ) {
        cyz = cz[l]*cy[m];
        for ( n = 0; n < 4; n ++, o ++ ) {
          res += cx[n]*cyz*tbuf[ offset[ o ] ];
        }
      }
      *rbuf = (float)res;
      continue;
    } 

    /* 
     * here, we are sure we are on some border
     * if ix == tdimx3  => ix = theCoeff->theDim[0]-3
     */
    
    *rbuf = 0.0;

    if ( ix == tdimx3 ) {
      if ( iy == tdimy3 ) {
        if ( iz == tdimz3 ) {

          for ( res = 0.0, o = 0, l = 0; l < 3; l++, o+=4 )
          for ( m = 0; m < 3; m++, o++ ) {
            cyz = cz[l]*cy[m];
            for ( n = 0; n < 3; n ++, o ++ )
              res += cx[n]*cyz*tbuf[ offset[ o ] ];
          }
          *rbuf = (float)res;
          continue;
        }
        for ( res = 0.0, o = 0, l = 0; l < 4; l++, o+=4 )
        for ( m = 0; m < 3; m++, o++ ) {
          cyz = cz[l]*cy[m];
          for ( n = 0; n < 3; n ++, o ++ )
            res += cx[n]*cyz*tbuf[ offset[ o ] ];
        }
        *rbuf = (float)res;
        continue;
      }

      if ( iz == tdimz3 ) {
        for ( res = 0.0, o = 0, l = 0; l < 3; l++ )
        for ( m = 0; m < 4; m++, o++ ) {
          cyz = cz[l]*cy[m];
          for ( n = 0; n < 3; n ++, o ++ )
            res += cx[n]*cyz*tbuf[ offset[ o ] ];
        }
        *rbuf = (float)res;
        continue;
      }
      for ( res = 0.0, o = 0, l = 0; l < 4; l++ )
      for ( m = 0; m < 4; m++, o++ ) {
        cyz = cz[l]*cy[m];
        for ( n = 0; n < 3; n ++, o ++ )
          res += cx[n]*cyz*tbuf[ offset[ o ] ];
      }
      *rbuf = (float)res;
      continue;
      
    }

    if ( iy == tdimy3 ) {

      if ( iz == tdimz3 ) {
        for ( res = 0.0, o = 0, l = 0; l < 3; l++, o+=4 )
        for ( m = 0; m < 3; m++ ) {
          cyz = cz[l]*cy[m];
          for ( n = 0; n < 4; n ++, o ++ )
            res += cx[n]*cyz*tbuf[ offset[ o ] ];
        }
        *rbuf = (float)res;
        continue;
      }
      for ( res = 0.0, o = 0, l = 0; l < 4; l++, o+=4 )
      for ( m = 0; m < 3; m++ ) {
        cyz = cz[l]*cy[m];
        for ( n = 0; n < 4; n ++, o ++ )
          res += cx[n]*cyz*tbuf[ offset[ o ] ];
      }
      *rbuf = (float)res;
      continue;
      
    }

    if ( iz == tdimz3 ) {
      for ( res = 0.0, o = 0, l = 0; l < 3; l++ )
      for ( m = 0; m < 4; m++ ) {
        cyz = cz[l]*cy[m];
        for ( n = 0; n < 4; n ++, o ++ )
          res += cx[n]*cyz*tbuf[ offset[ o ] ];
      }
      *rbuf = (float)res;
      continue;
    }
    for ( res = 0.0, o = 0, l = 0; l < 4; l++ )
    for ( m = 0; m < 4; m++ ) {
      cyz = cz[l]*cy[m];
      for ( n = 0; n < 4; n ++, o ++ )
        res += cx[n]*cyz*tbuf[ offset[ o ] ];
    }
    *rbuf = (float)res;
    continue;
    
  }
  
  return( 1 );
}





















int Reech2DCSpline4x4WithCoefficients( typeCSplineCoefficients *theCoeff,
                       r32* resBuf,  /* result buffer */
                       int *resDim,  /* dimensions of this buffer */
                       double *mat,
                       int *derivative )
{
  char *proc = "Reech2DCSpline4x4";
  register int i, j, ix, iy;
  register double x, y;
  register double dx, dy;
  register double res;
  
  int rdimx=resDim[0];
  int rdimy=resDim[1];

  register r32* rbuf = resBuf;
  register float *tbuf;
  float *bufCoeff = theCoeff->theCoeff;

  int tdimx  = theCoeff->theDim[0];

  int tdimx3 = theCoeff->theDim[0]-3;
  int tdimy3 = theCoeff->theDim[1]-3;
 
  register double ddimx = (double)theCoeff->theDim[0]-2.0-0.5;
  register double ddimy = (double)theCoeff->theDim[1]-2.0-0.5;

  double cx[4], cy[4];
  int offset[16];

  int l, m, n, o;

  double (*bsplinx)(double);
  double (*bspliny)(double);


  if ( derivative == (int*) NULL ) {
    bsplinx = & CubicBspline;
    bspliny = & CubicBspline;
  } else {
    switch ( derivative[0] ) {
    default :
      if ( _verbose_ ) 
        fprintf( stderr, "%s: invalid derivation order in x\n", proc );
      return( -1 );
    case 0 : bsplinx = & CubicBspline;            break;
    case 1 : bsplinx = & CubicBspline_FirstDeriv; break;
    case 2 : bsplinx = & CubicBspline_SecndDeriv; break;
    }
    
    switch ( derivative[1] ) {
    default :
      if ( _verbose_ ) 
        fprintf( stderr, "%s: invalid derivation order in y\n", proc );
      return( -1 );
    case 0 : bspliny = & CubicBspline;            break;
    case 1 : bspliny = & CubicBspline_FirstDeriv; break;
    case 2 : bspliny = & CubicBspline_SecndDeriv; break;
    }
    
  }
  
  

  for ( o = 0, m = 0; m < 4; m ++ )
  for ( n = 0; n < 4; n ++, o ++ ) {
    offset[o] = m * tdimx + n;
  }
  
  

  for ( j = 0; j < rdimy; j ++ )
  for ( i = 0; i < rdimx; i ++, rbuf ++ ) {
    /* computation of the corresponding point after deformation 
       
       (x,y,x) -> (ix=(int)x,iy=(int)x,iz=(int)x)
       les coefficients de la spline a prendre en compte
       sont [1+ix+{-1,0,1,2}, 1+iy+{-1,0,1,2}, 1+iz+{-1,0,1,2}]
       
       On ajoute 1 aux coordonnees, car on a ajoute
       une bordure aux coefficients de la spline.
     */
    x = mat[0] * i +  mat[1] * j              + mat[3];
    /*
    if ((x < 0) || ( x > tdimx3)) { *rbuf = 0; continue; }
    */
    if ((x < -0.5) || ( x > ddimx)) { *rbuf = 0; continue; }

    y = mat[4] * i +  mat[5] * j              + mat[7];
    /*
    if ((y < 0) || ( y > tdimy3)) { *rbuf = 0; continue; }
    */
    if ((y < -0.5) || ( y > ddimy)) { *rbuf = 0; continue; }
    
    /* ici on a 0 <   x   < theCoeff->theDim[0]-3
                0 <= ix   < theCoeff->theDim[0]-3
                1 <= ix+1 < theCoeff->theDim[0]-2
                0 <= ix+1-1 && ix+1+2 < theCoeff->theDim[0]
       maintenant on a -0.5 <         x   <        theCoeff->theDim[0]-2.5
                          0 <= ix   <= theCoeff->theDim[0]-3
                          1 <= ix+1 <= theCoeff->theDim[0]-2

       les extrema pour les coefficients verifient
                          0 <= ix+1-1 <= theCoeff->theDim[0]-3
                          3 <= ix+1+2 <= theCoeff->theDim[0]
       dans ce dernier cas, on suppose que ca va etre 0 ...
     */
    ix = (int)x;    dx = x-ix;
    iy = (int)y;    dy = y-iy;



    if ( x < 0.0 ) { x = 0.0; ix = 0; dx = 0.0; }
    if ( y < 0.0 ) { y = 0.0; iy = 0; dy = 0.0; }

    
    /* Pre-calcul des valeurs de bsplines pour chaque dimension:
       
       (*bsplinx)(dx-1), (*bsplinx)(dx), (*bsplinx)(dx+1), (*bsplinx)(dx+2)  
       (*bspliny)(dy-1), (*bspliny)(dy), (*bspliny)(dy+1), (*bspliny)(dy+2)  
       
       A. Roche
       
       Dans l'image des coefficients,
       (*bsplinx)(dx-1) correspond a ix, 
       (*bsplinx)(dy-1) correspond a iy, 

    */
    for ( l = 0; l < 4; l ++ ) {
      cx[l] = (*bsplinx)( dx - (double)(l-1) );
      cy[l] = (*bspliny)( dy - (double)(l-1) );
    }



    tbuf = bufCoeff;
    tbuf += iy * tdimx + ix;
    
    /* 
     * are we on the border or not ? 
     */
    if ( (ix < tdimx3) && (iy < tdimy3) ) {
      

      for ( res = 0.0, o = 0, m = 0; m < 4; m ++ )
      for ( n = 0; n < 4; n ++, o ++ )
        res += cx[n]*cy[m]*tbuf[ offset[ o ] ];
      *rbuf = (float)res;
      continue;
    }

    /* 
     * here, we are sure we are on some border
     * if ix == tdimx3  => ix = theCoeff->theDim[0]-3
     */
    
    *rbuf = 0.0;
    
    if ( ix == tdimx3 ) {
      if ( iy == tdimy3 ) {

        for ( res = 0.0, o = 0, m = 0; m < 3; m ++, o++ )
        for ( n = 0; n < 3; n ++, o ++ )
          res += cx[n]*cy[m]*tbuf[ offset[ o ] ];
        *rbuf = (float)res;
        continue;

      }

      for ( res = 0.0, o = 0, m = 0; m < 4; m ++, o++ )
      for ( n = 0; n < 3; n ++, o ++ )
        res += cx[n]*cy[m]*tbuf[ offset[ o ] ];
      *rbuf = (float)res;
      continue;

    }

    if ( iy == tdimy3 ) {
      for ( res = 0.0, o = 0, m = 0; m < 3; m ++ )
      for ( n = 0; n < 4; n ++, o ++ )
        res += cx[n]*cy[m]*tbuf[ offset[ o ] ];
      *rbuf = (float)res;
      continue;
    }
    
    for ( res = 0.0, o = 0, m = 0; m < 4; m ++ )
    for ( n = 0; n < 4; n ++, o ++ )
        res += cx[n]*cy[m]*tbuf[ offset[ o ] ];
    *rbuf = (float)res;
    continue;
   

  }
  
  return( 1 );
}






















static int _ConvertSlice( r32* resSlice, int slice,
                          void* resBuf, bufferType resType, int *resDim )
{
  char *proc = "_ConvertSlice";
  size_t dimxy = resDim[0]*resDim[1];
  int z;

  switch( resType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not handled in switch\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 *tmpBuf = (u8*)resBuf;
      for ( z=0; z<slice; z++ ) tmpBuf += dimxy;
      if ( ConvertBuffer( resSlice, FLOAT, tmpBuf, resType, dimxy ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert such image type\n", proc );
        return( -1 );
      }
    }
    break;
  case SCHAR :
    {
      s8 *tmpBuf = (s8*)resBuf;
      for ( z=0; z<slice; z++ ) tmpBuf += dimxy;
      if ( ConvertBuffer( resSlice, FLOAT, tmpBuf, resType, dimxy ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert such image type\n", proc );
        return( -1 );
      }
    }
    break;
  case USHORT :
    {
      u16 *tmpBuf = (u16*)resBuf;
      for ( z=0; z<slice; z++ ) tmpBuf += dimxy;
      if ( ConvertBuffer( resSlice, FLOAT, tmpBuf, resType, dimxy ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert such image type\n", proc );
        return( -1 );
      }
    }
    break;
  case SSHORT :
    {
      s16 *tmpBuf = (s16*)resBuf;
      for ( z=0; z<slice; z++ ) tmpBuf += dimxy;
      if ( ConvertBuffer( resSlice, FLOAT, tmpBuf, resType, dimxy ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert such image type\n", proc );
        return( -1 );
      }
    }
    break;
  case FLOAT :
    {
      r32 *tmpBuf = (r32*)resBuf;
      for ( z=0; z<slice; z++ ) tmpBuf += dimxy;
      if ( ConvertBuffer( resSlice, FLOAT, tmpBuf, resType, dimxy ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert such image type\n", proc );
        return( -1 );
      }
    }
    break;
  }

  return( 1 );
}










int ReechCSpline4x4( void* theBuf, bufferType theType, int *theDim,
                     void* resBuf, bufferType resType, int *resDim,
                     double *mat,
                     int *derivative )
{
  char *proc = "Reech3DCSpline4x4";
  typeCSplineCoefficients *c;
  r32 *resSlice = NULL;



  if ( theDim[2] == 1 && resDim[2] > 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, input is 2D and output is 3D\n", proc );
    return( -1 );
  }



  c = ComputeCSplineCoefficients( theBuf, theType, theDim );
  if ( c == NULL ) {
    if( _verbose_ )
      fprintf( stderr, "%s: unable to compute cubic spline coefficients.\n", proc );
    return( -1 );
  }
  


  /* 2D case
   */
  if ( theDim[2] == 1 ) {

    if ( resType != FLOAT ) {
      resSlice = (r32*)vtmalloc( resDim[0]*resDim[1]*sizeof( r32 ), "resSlice", proc );
      if ( resSlice == NULL ) {
        FreeTypeCSplineCoefficients( &c );
        if( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary slice\n", proc );
        return( -1 );
      }
    }

    if ( Reech2DCSpline4x4WithCoefficients( c, resSlice, resDim, mat,
                                            derivative ) != 1 ) {
      FreeTypeCSplineCoefficients( &c );
      if ( resType != FLOAT ) vtfree( resSlice );
      if( _verbose_ )
        fprintf( stderr, "%s: unable to compute 2D slice\n", proc );
    }

    if ( _ConvertSlice( resSlice, 0, resBuf, resType, resDim ) != 1 ) {
      FreeTypeCSplineCoefficients( &c );
      if ( resType != FLOAT ) vtfree( resSlice );
      if( _verbose_ )
        fprintf( stderr, "%s: unable to compute 2D slice\n", proc );
    }

    if ( resType != FLOAT ) vtfree( resSlice );

  }

  /* 3D case
   */
  else {

    switch( resType ) {
    default :
      if ( _verbose_ )
        fprintf(stderr, "%s: such type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      if ( Cspline3D4x4_u8( c, resBuf, resDim, mat, derivative) != 1 ) {
        if ( _verbose_ )
          fprintf(stderr, "%s: error when resampling\n", proc );
        return( -1 );
      }
      break;
    case SCHAR :
      if ( Cspline3D4x4_s8( c, resBuf, resDim, mat, derivative) != 1 ) {
        if ( _verbose_ )
          fprintf(stderr, "%s: error when resampling\n", proc );
        return( -1 );
      }
      break;
    case USHORT :
      if ( Cspline3D4x4_u16( c, resBuf, resDim, mat, derivative) != 1 ) {
        if ( _verbose_ )
          fprintf(stderr, "%s: error when resampling\n", proc );
        return( -1 );
      }
      break;
    case SSHORT :
      if ( Cspline3D4x4_s16( c, resBuf, resDim, mat, derivative) != 1 ) {
        if ( _verbose_ )
          fprintf(stderr, "%s: error when resampling\n", proc );
        return( -1 );
      }
      break;
    case FLOAT :
      if ( Cspline3D4x4_r32( c, resBuf, resDim, mat, derivative) != 1 ) {
        if ( _verbose_ )
          fprintf(stderr, "%s: error when resampling\n", proc );
        return( -1 );
      }
      break;
    }

  }
  


  FreeTypeCSplineCoefficients( &c );

  return( 1 );
}
