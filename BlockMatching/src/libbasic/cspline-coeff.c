/*************************************************************************
 * cspline-coeff.c - Cubic splines coefficients
 *
 * Copyright (c) INRIA 2017
 *
 * AUTHOR:
 * Alexis Roche
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  4 mai 2017 11:14:19 CEST
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

#include <cspline-coeff.h>



static int _verbose_ = 1;

void setVerboseInCsplineCoeff( int v )
{
  _verbose_ = v;
}

void incrementVerboseInCsplineCoeff( )
{
  _verbose_ ++;
}

void decrementVerboseInCsplineCoeff( )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}






void InitTypeCSplineCoefficients( typeCSplineCoefficients *t )
{
  t->theDim[0] = 0;
  t->theDim[1] = 0;
  t->theDim[2] = 0;
  t->theCoeff = NULL;
}



void FreeTypeCSplineCoefficients( typeCSplineCoefficients **t )
{
  vtfree( (*t)->theCoeff );
  vtfree( *t );
  *t = NULL;
}








/* Interpolation cubique spline d'un signal S, en supposant 
   les points d'echantillonnage regulierement espaces. 

   Les 2 premiers arguments sont des pointeurs sur tableaux 
   de float SUPPOSES DEJA ALLOUES de taille N.
        - s: valeurs d'echantillon
        - c: coefficients B-spline a calucler

   L'algo est recursif, donc tres rapide par rapport a une 
   methode matricielle. Il est decrit dans:
   M. Unser, "Splines : a perfect fit for signal/image processing ", 
   IEEE Signal Processing Magazine, Nov. 1999, in press. 

   A. Roche
*/

static void CubicSpline_Transform ( double * s, double * c, int N )
{

  int k;
  double cp, cm;
  double z1_k;
  double * bufs, * bufc;
  
  /* Constants */
  const double  z1 = -0.26794919243112270648; /* -2 + sqrt(3) */
  const double cz1 =  0.28867513459481288226; /* z1/(z1^2-1) */


  /*  (void)memset( (void*) c, 0, N * sizeof( double ) ); */
  
  /* Initial value for the causal recursion.
     We use a mirror symmetric boundary condition for the discrete signal,
     yielding:
     
     cp(0) = (1/2-z1^(2N-2)) \sum_{k=0}^{2N-3} s(k) z1^k s(k),
     
     where we set: s(N)=s(N-2), s(N+1)=s(N-3), ..., s(2N-3)=s(1).
  */

  bufs = s;
  cp = * bufs;
  z1_k = 1;

  for ( k = 1; k < N; k++ ) {
    
    z1_k = z1 * z1_k;   /* == z1^k */
    
    bufs ++;            /* pointe sur s[k] */
    cp += (*bufs) * z1_k;
  }
  
  /* Quand on arrive ici, z1_k = z1^(N-1) */
  for ( k = 2; k < N; k++ ) {
    
    z1_k = z1 * z1_k;  
    
    bufs --;
    cp += (* bufs) * z1_k;
  }
  
  /* Quand on arrive ici: z1_k = z1^(2N-3) */
  z1_k = z1 * z1_k;
  cp = cp / ( 1 - z1_k );


  /* Stockage du premier coeff causal */
  bufc = c;
  (*bufc) = cp;

  /* Do the causal recursion */
  bufs = s;
  for ( k = 1; k < N; k ++ ) {
    bufs ++;
    cp = (* bufs) + z1 * cp;
    bufc ++;
    (*bufc) = cp;
  }

  /* Initial value for the anticausal recursion. */
  cm = cz1 * ( 2.0 * cp - (*bufs) );
  (*bufc) = 6.0*cm;
 

  /* Do the anti causal recursion. */
  for ( k = (N-2); k >= 0; k=k-1 ) {
    bufc --;     /* bufc pointe sur l'indice k */
    cm = z1 * ( cm - (*bufc));
    (*bufc) = 6.0*cm;
  }

}






/**********************************************************************
 *
 * parallelism procedure
 *
 **********************************************************************/



typedef struct _CsplineCoeffParam {
  typeCSplineCoefficients *theCoeff;
  void *theBuf;
  bufferType theType;
  int *theDim;
} _CsplineCoeffParam;



static void *_CsplineXCoeff( void *par )
{
  char *proc = "_CsplineXCoeff";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _CsplineCoeffParam *p = (_CsplineCoeffParam *)parameter;

  typeCSplineCoefficients *theCoeff = p->theCoeff;
  void *theBuf = p->theBuf;
  bufferType theType = p->theType;
  int *theDim = p->theDim;

  int processCoeff = 0;
  int onlyCopy = 0;

  double *theLine = NULL;
  double *resLine = NULL;

  float *resBuf = theCoeff->theCoeff;
  size_t dimx, dimy, dimz;
  int cdimx, cdimxy;
  size_t dim, x, y, z;
  size_t offset = 0;



  /* here the input is the coefficient image
   * recall that its dimension are increased by 2
   */
  if ( theBuf == (void*)NULL || theType == TYPE_UNKNOWN ) {
    processCoeff = 1;
    theBuf = theCoeff->theCoeff;
    theType = FLOAT;
    theDim = theCoeff->theDim;
    dim = theDim[0] - 2;
    dimx = theDim[0] - 2;
    dimy = theDim[1] - 2;
    dimz = theDim[2] - 2;
  }
  else {
    dim = theDim[0];
    dimx = theDim[0];
    dimy = theDim[1];
    dimz = theDim[2];
  }



  cdimx = theCoeff->theDim[0];
  cdimxy = theCoeff->theDim[0] * theCoeff->theDim[1];



  /* auxiliary line allocation
   */
  theLine = (double*)vtmalloc( 2*dim*sizeof( double ), "theLine", proc );
  if ( theLine == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    chunk->ret = -1;
    return( (void*)NULL );
  }
  resLine = theLine;
  resLine += dim;


  /* skip first slice
   */
  resBuf += cdimxy;
  if ( processCoeff ) offset += cdimxy;


  for ( z=0; z<dimz; z++ ) {

    /* skip first row in coefficients
     * and in input if required
     */
    resBuf += cdimx;
    if ( processCoeff ) offset += cdimx;


    /* rows to be processed are  [first, last]
     * skip [0 first-1] rows
     */
    for ( y=0; y<first; y++ ) {
      resBuf += cdimx;
      if ( processCoeff ) offset += cdimx;
      else offset += dimx;
    }

    /* process [first, last] rows
     */
    for ( y=first; y<=last; y++ ) {
      resBuf ++;
      if ( processCoeff ) offset += 1;

#define _GETXLINE( TYPE ) {                       \
  TYPE *tmpBuf = (TYPE*)theBuf;                   \
  tmpBuf += offset;                               \
  for ( x=0; x<dim; x++ ) theLine[x] = tmpBuf[x]; \
}

      switch( theType ) {
      default :
        vtfree( theLine );
        if ( _verbose_ )
          fprintf( stderr, "%s: such type not handled yet\n", proc );
        chunk->ret = -1;
        return( (void*)NULL );
      case UCHAR :
        _GETXLINE( u8 );
        break;
      case SCHAR :
        _GETXLINE( s8 );
        break;
      case USHORT :
        _GETXLINE( u16 );
        break;
      case SSHORT :
        _GETXLINE( s16 );
        break;
      case FLOAT :
        _GETXLINE( r32 );
        break;
      }

      if ( onlyCopy ) {
        memcpy( resLine, theLine, dim*sizeof( double ) );
      } else {
        CubicSpline_Transform( theLine, resLine, dim );
      }

      for ( x=0; x<dim; x++ ) resBuf[x] = resLine[x];

      resBuf += dim;
      offset += dim;

      resBuf ++;
      if ( processCoeff ) offset += 1;
    }

    /* rows to be processed are  [first, last]
     * skip [last+1 dimy-1] rows
     */

    for ( y=last+1; y<dimy; y++ ) {
      resBuf += cdimx;
      if ( processCoeff ) offset += cdimx;
      else offset += dimx;
    }

    /* skip last row in coefficients
     * and in input if required
     */
    resBuf += cdimx;
    if ( processCoeff ) offset += cdimx;

  }


  vtfree( theLine );

  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_CsplineYCoeff( void *par )
{
  char *proc = "_CsplineYCoeff";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _CsplineCoeffParam *p = (_CsplineCoeffParam *)parameter;

  typeCSplineCoefficients *theCoeff = p->theCoeff;
  void *theBuf = p->theBuf;
  bufferType theType = p->theType;
  int *theDim = p->theDim;

  int processCoeff = 0;
  int onlyCopy = 0;

  double *theLine = NULL;
  double *resLine = NULL;

  float *resBuf = theCoeff->theCoeff;
  size_t dimx, dimy, dimz;
  int cdimx, cdimy, cdimxy;
  size_t dim, x, y, z;
  size_t offset = 0;



  /* here the input is the coefficient image
   * recall that its dimension are increased by 2
   */
  if ( theBuf == (void*)NULL || theType == TYPE_UNKNOWN ) {
    processCoeff = 1;
    theBuf = theCoeff->theCoeff;
    theType = FLOAT;
    theDim = theCoeff->theDim;
    dim = theDim[1] - 2;
    dimx = theDim[0] - 2;
    dimy = theDim[1] - 2;
    dimz = theDim[2] - 2;
  }
  else {
    dim = theDim[1];
    dimx = theDim[0];
    dimy = theDim[1];
    dimz = theDim[2];
  }



  cdimx = theCoeff->theDim[0];
  cdimy = theCoeff->theDim[1];
  cdimxy = theCoeff->theDim[0] * theCoeff->theDim[1];



  /* auxiliary line allocation
   */
  theLine = (double*)vtmalloc( 2*dim*sizeof( double ), "theLine", proc );
  if ( theLine == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    chunk->ret = -1;
    return( (void*)NULL );
  }
  resLine = theLine;
  resLine += dim;


  /* skip first slice
   */
  resBuf += cdimxy;
  if ( processCoeff ) offset += cdimxy;


  for ( z=0; z<dimz; z++ ) {

    /* advance to the second row in coefficients
     * and in input if required
     */
    resBuf += cdimx;
    if ( processCoeff ) offset += cdimx;

    /* skip first col in coefficient
     * and in input if required
     */
    resBuf ++;
    if ( processCoeff ) offset ++;

    /* cols to be processed are  [first, last]
     * skip [0 first-1] cols
     */
    for ( x=0; x<first; x++ ) {
      resBuf ++;
      offset ++;
    }

    /* process [first, last] rows
     */
    for ( x=first; x<=last; x++ ) {

#define _GETYLINE( TYPE ) {                       \
  TYPE *tmpBuf = (TYPE*)theBuf;                   \
  tmpBuf += offset;                               \
  if ( processCoeff ) {                           \
    for ( y=0; y<dim; y++ ) theLine[y] = tmpBuf[y*cdimx]; \
  }                                               \
  else {                                          \
    for ( y=0; y<dim; y++ ) theLine[y] = tmpBuf[y*dimx]; \
  }                                               \
}

      switch( theType ) {
      default :
        vtfree( theLine );
        if ( _verbose_ )
          fprintf( stderr, "%s: such type not handled yet\n", proc );
        chunk->ret = -1;
        return( (void*)NULL );
      case UCHAR :
        _GETYLINE( u8 );
        break;
      case SCHAR :
        _GETYLINE( s8 );
        break;
      case USHORT :
        _GETYLINE( u16 );
        break;
      case SSHORT :
        _GETYLINE( s16 );
        break;
      case FLOAT :
        _GETYLINE( r32 );
        break;
      }

      if ( onlyCopy ) {
        memcpy( resLine, theLine, dim*sizeof( double ) );
      } else {
        CubicSpline_Transform( theLine, resLine, dim );
      }

      for ( y=0; y<dim; y++ ) resBuf[y*cdimx] = resLine[y];

      resBuf ++;
      offset ++;
    }

    /* rows to be processed are  [first, last]
     * skip [last+1 dimy-1] rows
     */

    for ( x=last+1; x<dimx; x++ ) {
      resBuf ++;
      offset ++;
    }

    /* skip last row in coefficients
     * and in input if required
     */
    resBuf ++;
    if ( processCoeff ) offset ++;

    /* here, we are at the beginning of the third row in the coefficients
     * and either at the beginning of the third row in the input (if
     * equal to coefficients) or at the beginning of the second row
     *
     * go to the beginning of the next slice
     */
    resBuf += (cdimy-2)*cdimx;
    if ( processCoeff ) offset += (cdimy-2)*cdimx;
    else offset += (dimy-1)*dimx;
  }


  vtfree( theLine );

  chunk->ret = 1;
  return( (void*)NULL );
}





static void *_CsplineZCoeff( void *par )
{
  char *proc = "_CsplineZCoeff";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _CsplineCoeffParam *p = (_CsplineCoeffParam *)parameter;

  typeCSplineCoefficients *theCoeff = p->theCoeff;
  void *theBuf = p->theBuf;
  bufferType theType = p->theType;
  int *theDim = p->theDim;

  int processCoeff = 0;
  int onlyCopy = 0;

  double *theLine = NULL;
  double *resLine = NULL;

  float *resBuf = theCoeff->theCoeff;
  size_t dimx, dimy, dimxy;
  int cdimx, cdimxy;
  size_t dim, x, y, z;
  size_t offset = 0;



  /* here the input is the coefficient image
   * recall that its dimension are increased by 2
   */
  if ( theBuf == (void*)NULL || theType == TYPE_UNKNOWN ) {
    processCoeff = 1;
    theBuf = theCoeff->theCoeff;
    theType = FLOAT;
    theDim = theCoeff->theDim;
    dim = theDim[2] - 2;
    dimx = theDim[0] - 2;
    dimy = theDim[1] - 2;
  }
  else {
    dim = theDim[2];
    dimx = theDim[0];
    dimy = theDim[1];
  }



  cdimx = theCoeff->theDim[0];
  cdimxy = theCoeff->theDim[0] * theCoeff->theDim[1];
  dimxy = dimx * dimy;



  /* auxiliary line allocation
   */
  theLine = (double*)vtmalloc( 2*dim*sizeof( double ), "theLine", proc );
  if ( theLine == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    chunk->ret = -1;
    return( (void*)NULL );
  }
  resLine = theLine;
  resLine += dim;


  /* skip first slice
   */
  resBuf += cdimxy;
  if ( processCoeff ) offset += cdimxy;


  /* advance to the second row in coefficients
   * and in input if required
   */
  resBuf += cdimx;
  if ( processCoeff ) offset += cdimx;



  for ( y=0; y<dimy; y++ ) {

    /* skip first col in coefficient
     * and in input if required
     */
    resBuf ++;
    if ( processCoeff ) offset ++;

    /* z to be processed are  [first, last]
     * skip [0 first-1] z
     */

    for ( x=0; x<first; x++ ) {
      resBuf ++;
      offset ++;
    }

    /* process [first, last] rows
     */
    for ( x=first; x<=last; x++ ) {

#define _GETZLINE( TYPE ) {                       \
  TYPE *tmpBuf = (TYPE*)theBuf;                   \
  tmpBuf += offset;                               \
  if ( processCoeff ) {                           \
    for ( z=0; z<dim; z++ ) theLine[z] = tmpBuf[(size_t)z*(size_t)cdimxy]; \
  }                                               \
  else {                                          \
    for ( z=0; z<dim; z++ ) theLine[z] = tmpBuf[(size_t)z*(size_t)dimxy]; \
  }                                               \
}

      switch( theType ) {
      default :
        vtfree( theLine );
        if ( _verbose_ )
          fprintf( stderr, "%s: such type not handled yet\n", proc );
        chunk->ret = -1;
        return( (void*)NULL );
      case UCHAR :
        _GETZLINE( u8 );
        break;
      case SCHAR :
        _GETZLINE( s8 );
        break;
      case USHORT :
        _GETZLINE( u16 );
        break;
      case SSHORT :
        _GETZLINE( s16 );
        break;
      case FLOAT :
        /* _GETZLINE( r32 ); */
      {
        r32 *tmpBuf = (r32*)theBuf;
        tmpBuf += offset;
        if ( processCoeff ) {
          for ( z=0; z<dim; z++ ) theLine[z] = tmpBuf[(size_t)z*(size_t)cdimxy];
        }
        else {
          for ( z=0; z<dim; z++ ) theLine[z] = tmpBuf[(size_t)z*(size_t)dimxy];
        }
      }
        break;
      }

      if ( onlyCopy ) {
        memcpy( resLine, theLine, dim*sizeof( double ) );
      } else {
        CubicSpline_Transform( theLine, resLine, dim );
      }

      for ( z=0; z<dim; z++ ) resBuf[(size_t)z*(size_t)cdimxy] = resLine[z];

      resBuf ++;
      offset ++;
    }

    /* rows to be processed are  [first, last]
     * skip [last+1 dimy-1] rows
     */

    for ( x=last+1; x<dimx; x++ ) {
      resBuf ++;
      offset ++;
    }

    /* skip last row in coefficients
     * and in input if required
     */
    resBuf ++;
    if ( processCoeff ) offset ++;
  }


  vtfree( theLine );

  chunk->ret = 1;
  return( (void*)NULL );
}





static int _Compute1DCSplineCoefficients( typeCSplineCoefficients *theCoeff,
                                          void *theBuf,
                                          bufferType theType,
                                          int *theDim,
                                          int direction )
{
  char *proc = "_Compute1DCSplineCoefficients";
  int dim, i;
  typeChunks chunks;
  _chunk_callfunction ftn;
  _CsplineCoeffParam p;



  if ( direction < 0 || direction > 2 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: invalid direction\n", proc );
    return( -1 );
  }



  switch( direction ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: invalid direction\n", proc );
    return( -1 );
  case 0 :
    ftn = &_CsplineXCoeff;
    dim = theCoeff->theDim[1] - 2;
    break;
  case 1 :
    ftn = &_CsplineYCoeff;
    dim = theCoeff->theDim[0] - 2;
    break;
  case 2 :
    ftn = &_CsplineZCoeff;
    dim = theCoeff->theDim[0] - 2;
    break;
  }



  initChunks( &chunks );
  if ( buildChunks( &chunks, 0, dim-1, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }

  p.theCoeff = theCoeff;
  p.theBuf = theBuf;
  p.theType = theType;
  p.theDim = theDim;

  for ( i=0; i<chunks.n_allocated_chunks; i++ )
      chunks.data[i].parameters = (void*)(&p);

  /* processing
   */
  if ( processChunks( ftn, &chunks, proc ) != 1 ) {
    freeChunks( &chunks );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute spline coefficients along #%d\n", proc, direction );
    return( -1 );
  }

  freeChunks( &chunks );
  return( 1 );
}










/* calcul des coefficients de la spline cubique

   - cas 3D : theDim[2] > 1
     les coefficients sont un buffer 3D
     de taille (theDim[0]+2) * (theDim[1]+2) * (theDim[2]+2)

   - cas 2D : theDim[2] <= 1 (donc theDim[2]==1)
     les coefficients sont un buffer 3D
     de taille (theDim[0]+2) * (theDim[1]+2) * theDim[2]

   Aux bords, on duplique les valeurs calculees.
*/

typeCSplineCoefficients *ComputeCSplineCoefficients( void *theBuf,
                                                     bufferType theType,
                                                     int *theDim )
{
  char *proc = "ComputeCSplineCoefficients";
  typeCSplineCoefficients *theCSpline = NULL;

  int min_elements = getMinElementsInChunks();

  int dimy = theDim[1];
  int dimz = theDim[2];

  int dx2  = theDim[0]+2;
  int dxy2  = (theDim[0]+2)*(theDim[1]+2);

  int y, z;
  float *resBuf;


  setMinElementsInChunks( 1 );

  theCSpline = (typeCSplineCoefficients *)vtmalloc( sizeof(typeCSplineCoefficients),
                                                    "theCSpline", proc );
  if ( theCSpline == (typeCSplineCoefficients *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate cubic spline coefficients structure\n",
               proc );
    return( NULL );
  }


  theCSpline->theDim[0] = theDim[0] + 2;
  theCSpline->theDim[1] = theDim[1] + 2;
  if ( theDim[2] > 1 ) theCSpline->theDim[2] = theDim[2] + 2;
  else                 theCSpline->theDim[2] = 1;




  /* on alloue un buffer de taille (dimx+2)*(dimy+2)*(dimz+2)
     dans le cas 3D, ou (dimx+2)*(dimy+2) dans le cas 2D
  */

  theCSpline->theCoeff = (float*)vtmalloc( (size_t)theCSpline->theDim[0] *
                                           (size_t)theCSpline->theDim[1] *
                                           (size_t)theCSpline->theDim[2] * (size_t)sizeof( float ),
                                           "theCSpline->theCoeff", proc );
  if ( theCSpline->theCoeff == (float*)NULL ) {
    vtfree( theCSpline );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate cubic spline coefficients buffer\n",
               proc );
    return( NULL );
  }


  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: processing X direction\n", proc );
  }
  if ( _Compute1DCSplineCoefficients( theCSpline, theBuf, theType, theDim, 0 ) != 1 ) {
    FreeTypeCSplineCoefficients( &theCSpline );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing spline coefficients along X\n", proc );
    return( NULL );
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: processing Y direction\n", proc );
  }
  if ( _Compute1DCSplineCoefficients( theCSpline, (void*)NULL, TYPE_UNKNOWN, (int*)NULL, 1 ) != 1 ) {
    FreeTypeCSplineCoefficients( &theCSpline );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing spline coefficients along Y\n", proc );
    return( NULL );
  }

  if ( theDim[2] > 1 ) {

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: processing Z direction\n", proc );
    }
    if ( _Compute1DCSplineCoefficients( theCSpline, (void*)NULL, TYPE_UNKNOWN, (int*)NULL, 2 ) != 1 ) {
      FreeTypeCSplineCoefficients( &theCSpline );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when computing spline coefficients along Z\n", proc );
      return( NULL );
    }
  }


  /* on fait l'effet miroir => on duplique les bords
     [0] = [2] et [dimx-1] = [dimx-3]
     1. on fait les bords X et Y dans les plans calcules
     2. on ajoute des plans selon Z dans le cas 3D
   */
  resBuf = theCSpline->theCoeff;
  if ( theDim[2] > 1 ) resBuf += dxy2;

  for ( z=0; z<dimz; z++, resBuf+=dxy2 ) {
    for ( y=1; y<=dimy; y++ ) {
      resBuf[y*dx2]         = resBuf[y*dx2+2];
      resBuf[y*dx2 + dx2-1] = resBuf[y*dx2 + dx2-3];
    }
    memcpy( &(resBuf[0]), &(resBuf[2*dx2]), dx2*sizeof(float) );
    memcpy( &(resBuf[(dimy+1)*dx2]), &(resBuf[(dimy-1)*dx2]), dx2*sizeof(float) );
  }

  if ( theDim[2] > 1 ) {
    resBuf = theCSpline->theCoeff;
    memcpy( &(resBuf[0]), &(resBuf[2*dxy2]), dxy2*sizeof(float) );
    memcpy( &(resBuf[(size_t)(dimz+1)*(size_t)dxy2]), &(resBuf[(size_t)(dimz-1)*(size_t)dxy2]), dxy2*sizeof(float) );
  }


  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: done\n", proc );
  }

  setMinElementsInChunks( min_elements );

  return( theCSpline );
}



