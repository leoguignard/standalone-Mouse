/****************************************************
 * topological-component-image.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Ven  9 jui 2017 16:05:53 CEST
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
#include <math.h>
#include <time.h>
#include <string.h>



#include <vtmalloc.h>
#include <connexe.h>



#include <topological-component-image.h>








/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalComponentImage( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalComponentImage(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalComponentImage(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalComponentImage( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalComponentImage(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalComponentImage(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}














/**************************************************
 *
 *
 *
 **************************************************/



int countNeighborsInImage( void *theBuf, bufferType theType,
                               void *resBuf, bufferType resType, int *theDim )
{
  char *proc = "_countNeighborsInImage";
  int x, y, z;
  int i, j, k;
  int n;
  int neighbors;

  /* 2D case
   */
#define _NEIGHBORS2D_THETYPE( TYPE ) {               \
  TYPE *buf = (TYPE*)theBuf;                         \
  for ( n=0, y=0; y<theDim[1]; y++ )                 \
  for ( x=0; x<theDim[0]; x++, n++ ) {               \
    if ( buf[n] == 0 ) {                             \
      res[n] = 0;                                    \
      continue;                                      \
    }                                                \
    neighbors = 0;                                   \
    for ( j=-1; j<=1; j++ ) {                        \
      if ( y+j<0 || y+j>=theDim[1] ) continue;       \
      for ( i=-1; i<=1; i++ ) {                      \
        if ( x+i<0 || x+i>=theDim[0] ) continue;     \
        if ( j == 0 && i == 0 ) continue;            \
        if ( buf[n+j*theDim[0]+i] > 0 ) neighbors++; \
      }                                              \
    }                                                \
    res[n] = neighbors;                              \
  }                                                  \
}


#define _NEIGHBORS3D_THETYPE( TYPE ) {               \
  TYPE *buf = (TYPE*)theBuf;                         \
  for ( n=0, z=0; z<theDim[2]; z++ )                 \
  for ( y=0; y<theDim[1]; y++ )                      \
  for ( x=0; x<theDim[0]; x++, n++ ) {               \
    if ( buf[n] == 0 ) {                             \
      res[n] = 0;                                    \
      continue;                                      \
    }                                                \
    neighbors = 0;                                   \
    for ( k=-1; k<=1; k++ ) {                        \
      if ( z+k<0 || z+k>=theDim[2] ) continue;       \
      for ( j=-1; j<=1; j++ ) {                      \
        if ( y+j<0 || y+j>=theDim[1] ) continue;     \
        for ( i=-1; i<=1; i++ ) {                    \
          if ( x+i<0 || x+i>=theDim[0] ) continue;   \
          if ( k == 0 && j == 0 && i == 0 ) continue; \
          if ( buf[n+(k*theDim[1]+j)*theDim[0]+i] > 0 ) neighbors++; \
        }                                            \
      }                                              \
    }                                                \
    res[n] = neighbors;                              \
  }                                                  \
}


#define _NEIGHBORS2D_RESTYPE( TYPE ) { \
  TYPE *res = (TYPE*)resBuf;           \
  switch( theType ) {                  \
  default :                            \
    if ( _verbose_ )                   \
      fprintf( stderr, "%s: such input image type not handled yet\n", proc ); \
    return( -1 );                      \
  case UCHAR :                         \
    _NEIGHBORS2D_THETYPE( u8 );        \
    break;                             \
  case SSHORT :                        \
    _NEIGHBORS2D_THETYPE( s16 );       \
    break;                             \
  case USHORT :                        \
    _NEIGHBORS2D_THETYPE( u16 );       \
    break;                             \
  }                                    \
}

#define _NEIGHBORS3D_RESTYPE( TYPE ) { \
  TYPE *res = (TYPE*)resBuf;           \
  switch( theType ) {                  \
  default :                            \
    if ( _verbose_ )                   \
      fprintf( stderr, "%s: such input image type not handled yet\n", proc ); \
    return( -1 );                      \
  case UCHAR :                         \
    _NEIGHBORS3D_THETYPE( u8 );        \
    break;                             \
  case SSHORT :                        \
    _NEIGHBORS3D_THETYPE( s16 );       \
    break;                             \
  case USHORT :                        \
    _NEIGHBORS3D_THETYPE( u16 );       \
    break;                             \
  }                                    \
}



  /* 2D case
   */
  if ( theDim[2] == 1 ) {
    switch ( resType ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such output image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _NEIGHBORS2D_RESTYPE( u8 );
      break;
    }
  }
  /* 3D case
   */
  else {
    switch ( resType ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such output image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _NEIGHBORS3D_RESTYPE( u8 );
      break;
    }
  }

  return( 1 );
}








/**************************************************
 *
 * Component image
 * this is an image where components are labeled
 * first the edges, then the junctions
 *
 **************************************************/



static void initComponent( typeComponent *c )
{
  c->label = -1;
  c->size = 0;
  c->type = _UNKNOWN_COMPONENT_;
}



void initComponentImage( typeComponentImage *c )
{
  c->component = (typeComponent*)NULL;
  c->n_component = 0;
  c->n_allocated_component = 0;

  c->componentBuf = (void*)NULL;
  c->componentBufferType = TYPE_UNKNOWN;
  c->theDim[0] = 0;
  c->theDim[1] = 0;
  c->theDim[2] = 0;
  c->voxelSize[0] = 1.0;
  c->voxelSize[1] = 1.0;
  c->voxelSize[2] = 1.0;
}



int allocComponentImage( typeComponentImage *c, bufferType type, int *theDim )
{
  char *proc = "allocComponentImage";
  size_t i, v;

  v = (size_t)theDim[0] * (size_t)theDim[1] * (size_t)theDim[2];
  switch( type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 *theBuf;
      c->componentBuf = (void*)vtmalloc( v*sizeof(u8),
                                         "c->componentBuf", proc );
      if ( c->componentBuf == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
        return( -1 );
      }
      theBuf = (u8*)c->componentBuf;
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = 0;
    }
    break;
  case SSHORT :
    {
      s16 *theBuf;
      c->componentBuf = (void*)vtmalloc( v*sizeof(s16),
                                         "c->componentBuf", proc );
      if ( c->componentBuf == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
        return( -1 );
      }
      theBuf = (s16*)c->componentBuf;
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = 0;
    }
    break;
  case USHORT :
    {
      u16 *theBuf;
      c->componentBuf = (void*)vtmalloc( v*sizeof(u16),
                                         "c->componentBuf", proc );
      if ( c->componentBuf == (void*)NULL ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: allocation failed\n", proc );
        return( -1 );
      }
      theBuf = (u16*)c->componentBuf;
      for ( i=0; i<v; i++, theBuf++ )
        *theBuf = 0;
    }
    break;
  }

  c->componentBufferType = type;
  c->theDim[0] = theDim[0];
  c->theDim[1] = theDim[1];
  c->theDim[2] = theDim[2];

  return( 1 );
}





static int _size_to_be_allocated_ = 10;

static int _addComponentToComponentImage( typeComponentImage *theIm, int label, enumComponentType type )
{
  char * proc = "_addComponentToComponentImage";
  int i, s =  theIm->n_allocated_component;
  typeComponent *data;

  if ( label >= theIm->n_allocated_component ) {
    while ( label >= s ) s+= _size_to_be_allocated_;
    data = (typeComponent*)vtmalloc( s * sizeof(typeComponent), "data", proc );
    if ( data == (typeComponent*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( theIm->n_allocated_component > 0 ) {
      (void)memcpy( data, theIm->component, theIm->n_allocated_component*sizeof(typeComponent) );
      vtfree( theIm->component );
    }
    for ( i = theIm->n_allocated_component; i<s; i++ ) {
      initComponent( &(data[i]) );
    }
    theIm->n_allocated_component = s;
    theIm->component = data;
  }

  if ( theIm->component[label].type != _UNKNOWN_COMPONENT_ ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, component already typed\n", proc );
  }
  theIm->component[label].type = type;
  if ( label >= theIm->n_component )
    theIm->n_component  = label + 1;

  return( 1 );
}





void freeComponentImage( typeComponentImage *c )
{
  if ( c->component != (typeComponent*)NULL )
      vtfree( c->component );
  if ( c->componentBuf != (void*)NULL )
      vtfree( c->componentBuf );
  initComponentImage( c );
}





/**************************************************
 *
 * *
 *
 **************************************************/



static int fprintfComponentImage( FILE *f, typeComponentImage *c )
{
  char *proc = "fprintfComponentImage";
  int i, v;

  v = c->theDim[0] * c->theDim[1] * c->theDim[2];
  for ( i=0; i<c->n_component; i++ )
    c->component[i].size = 0;

#define FPRINTFComponentImage( TYPE ) {     \
  TYPE *theBuf = (TYPE*)c->componentBuf;    \
  for ( i=0; i<v; i++ ) {                   \
    if ( theBuf[i] == 0 ) continue;         \
    c->component[ theBuf[i] ].size ++;      \
  }                                         \
}

  switch( c->componentBufferType ) {
  default :
    if ( _verbose_ )
        fprintf( stderr, "%s: such component type not handled yet\n", proc );
      return( -1 );
  case UCHAR :
    FPRINTFComponentImage( u8 );
    break;
  case SSHORT :
    FPRINTFComponentImage( s16 );
    break;
  case USHORT :
    FPRINTFComponentImage( u16 );
    break;
  }

  for ( i=1; i<c->n_component; i++ ) {
    fprintf( f, "#%3d ", i );
    switch( c->component[i].type ) {
    default :                   fprintf( f, "- default  - "); break;
    case _UNKNOWN_COMPONENT_ :  fprintf( f, "- unknown  - "); break;
    case _EDGE_COMPONENT_ :     fprintf( f, "- edge     - "); break;
    case _JUNCTION_COMPONENT_ : fprintf( f, "- junction - "); break;
    }
    fprintf( f, "has size %5d\n", c->component[i].size );
  }

  return( 1 );
}



static void fprintfSummaryComponentImage( FILE *f, typeComponentImage *c, char *desc  )
{
  char *proc = "fprintfSummaryComponentImage";
  int i;
  int nedge, fedge, ledge;
  int nunknown, funknown, lunknown;
  int njunction, fjunction, ljunction;

  fedge = ledge = -1;
  funknown = lunknown = -1;
  fjunction = ljunction = -1;

  nedge = nunknown = njunction = 0;

  for ( i=1; i<c->n_component; i++ ) {
    switch( c->component[i].type ) {
    default :
    case _UNKNOWN_COMPONENT_ :
      if ( nunknown == 0 ) {
        funknown = lunknown = i;
      }
      else {
        if ( funknown > i ) funknown = i;
        else if ( lunknown < i ) lunknown = i;
      }
      nunknown ++;
      break;
    case _EDGE_COMPONENT_ :
      if ( nedge == 0 ) {
        fedge = ledge = i;
      }
      else {
        if ( fedge > i ) fedge = i;
        else if ( ledge < i ) ledge = i;
      }
      nedge ++;
      break;
    case _JUNCTION_COMPONENT_ :
      if ( njunction == 0 ) {
        fjunction = ljunction = i;
      }
      else {
        if ( fjunction > i ) fjunction = i;
        else if ( ljunction < i ) ljunction = i;
      }
      njunction ++;
      break;
    }
  }

  if ( desc != (char*)NULL )
    fprintf( f, "%s: ", desc );
  else
    fprintf( f, "%s: ", proc );

  fprintf( f, "%d edge labels in [%d %d]\n",
           nedge, fedge, ledge );
  fprintf( f, "\t %d junctions labels in [%d %d]\n",
           njunction, fjunction, ljunction );
  if ( nunknown > 0 ) {
    fprintf( f, "\t %d unknown labels in [%d %d]\n",
             nunknown, funknown, lunknown );
  }

}


/**************************************************
 *
 *
 *
 **************************************************/


/* the typeComponentImage *tree is supposed to be either
 * 0, or only having values at the same points than theBuf
 * (theBuf can be the componentBuf of tree)
 */
int imageToComponentImage( void *theBuf, bufferType theType,
                           typeComponentImage *tree )
{
  char *proc = "imageToComponentImage";
  int i, v;
  int theDim[3];
  unsigned char *tmpBuf;
  unsigned short int *tmpComponents;
  int lastedge, lastjunction;
  int njunctionpoints;

  /* allocations
   */
  for (i=0; i<3; i++) theDim[i] = tree->theDim[i];
  v = theDim[2] * theDim[1] * theDim[0];

  tmpBuf = (unsigned char*)vtmalloc( v * sizeof(unsigned char), "tmpBuf", proc );
  if ( tmpBuf == (unsigned char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary input buffer\n", proc );
    return( -1 );
  }

  tmpComponents = (unsigned short int*)vtmalloc( v * sizeof(unsigned short int),
                                                 "tmpComponents", proc );
  if ( tmpComponents == (unsigned short int*)NULL ) {
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary output buffer\n", proc );
    return( -1 );
  }



  /* tree branches
   */
  countNeighborsInImage( theBuf, theType, tmpBuf, UCHAR, theDim );
  for ( i=0; i<v; i++ ) {
      tmpComponents[i] = ( tmpBuf[i] == 1 || tmpBuf[i] == 2 ) ? 255 : 0;
  }
  lastedge = CountConnectedComponents( tmpComponents, USHORT, tmpComponents, USHORT, theDim );
  if ( lastedge == -1 ) {
    vtfree( tmpComponents );
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when counting tree branch components\n", proc );
    return( -1 );
  }
  if ( lastedge == 0 ) {
    vtfree( tmpComponents );
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: no tree branch components\n", proc );
    return( -1 );
  }
  if ( _debug_ ) {
      fprintf( stderr, "%s: found %d tree branches\n", proc, lastedge );
  }

#define _EDGE_COMPONENT_IMAGE_( TYPE, MAX ) { \
  TYPE *theCC = (TYPE*)tree->componentBuf;    \
  if ( lastedge > MAX ) {                     \
    vtfree( tmpComponents );                  \
    vtfree( tmpBuf );                         \
    if ( _verbose_ )                          \
      fprintf( stderr, "%s: too many components for the component type\n", proc ); \
    return( -1 );                             \
  }                                           \
  for ( i=0; i<v; i++ ) {                     \
      theCC[i] = tmpComponents[i];            \
  }                                           \
}

  switch( tree->componentBufferType ) {
  default :
    vtfree( tmpComponents );
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: such component type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _EDGE_COMPONENT_IMAGE_( u8, 255 );
    break;
  case SSHORT :
    _EDGE_COMPONENT_IMAGE_( s16, 32767 );
    break;
  case USHORT :
    _EDGE_COMPONENT_IMAGE_( u16, 65535 );
    break;
  }

  for ( i=1; i<=lastedge; i++ ) {
    if ( _addComponentToComponentImage( tree, i, _EDGE_COMPONENT_ ) != 1 ) {
      vtfree( tmpComponents );
      vtfree( tmpBuf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add component #%d (edge) to component image\n", proc, i );
      return( -1 );
    }
  }



  /* tree junctions
   */
  for ( i=0, njunctionpoints=0; i<v; i++ ) {
    if ( tmpBuf[i] >= 3 ) {
      tmpComponents[i] = 255;
      njunctionpoints ++;
    }
    else {
      tmpComponents[i] = 0;
    }
  }

  vtfree( tmpBuf );

  /* no junctions
   */
  if ( njunctionpoints == 0 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: warning, no tree branch junctions\n", proc );
    }
    vtfree( tmpComponents );
    return( 1 );
  }

  lastjunction = CountConnectedComponents( tmpComponents, USHORT, tmpComponents, USHORT, theDim );
  if ( lastjunction == -1 ) {
    vtfree( tmpComponents );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when counting tree branch junctions\n", proc );
    return( -1 );
  }
  if ( lastjunction == 0 ) {
    vtfree( tmpComponents );
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, this should not happen (no tree branch junctions)\n", proc );
    return( -1 );
  }
  if ( _debug_ ) {
      fprintf( stderr, "%s: found %d tree junctions\n", proc, lastjunction );
  }

#define _JUNCTION_COMPONENT_IMAGE_( TYPE, MAX ) { \
  TYPE *theCC = (TYPE*)tree->componentBuf;        \
  if ( lastjunction+lastedge > 255 ) {            \
    vtfree( tmpComponents );                      \
    if ( _verbose_ )                              \
      fprintf( stderr, "%s: too many components for the component type\n", proc ); \
    return( -1 );                                 \
  }                                               \
  for ( i=0; i<v; i++ ) {                         \
    if ( tmpComponents[i] > 0 )                   \
      theCC[i] = tmpComponents[i]+lastedge;       \
  }                                               \
}

  switch( tree->componentBufferType ) {
  default :
    vtfree( tmpComponents );
    if ( _verbose_ )
      fprintf( stderr, "%s: such component type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _JUNCTION_COMPONENT_IMAGE_( u8, 255 );
    break;
  case SSHORT :
    _JUNCTION_COMPONENT_IMAGE_( s16, 32767 );
    break;
  case USHORT :
    _JUNCTION_COMPONENT_IMAGE_( u16, 65535 );
    break;
  }

  for ( i=1; i<=lastjunction; i++ ) {
    if ( _addComponentToComponentImage( tree, i+lastedge, _JUNCTION_COMPONENT_ ) != 1 ) {
      vtfree( tmpComponents );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add component #%d (edge) to component image\n", proc, i );
      return( -1 );
    }
  }



  vtfree( tmpComponents );

  if ( _verbose_ >= 2 ) {
    fprintfSummaryComponentImage( stderr, tree, proc );
    if ( 0 && _debug_ )
      (void)fprintfComponentImage( stderr, tree );
  }

  return( 1 );
}











