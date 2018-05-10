

#include <string.h>
#include <strings.h>

#include <klb_Cwrapper.h>
#include <klb.h>


static int _verbose_ = 1;
static int _debug_ = 0;





/************************************************************
 *
 *
 *
 ************************************************************/



int testKlbHeader( char *magic __attribute__ ((unused)), const char *name __attribute__ ((unused)) )
{
    return -1;
}





/************************************************************
 *
 *
 *
 ************************************************************/

int readKlbImage( const char *name, _image *im )
{
  char *proc = "readKlbImage";
  uint32_t	xyzctR[KLB_DATA_DIMS];
  uint32_t	blockSizeR[KLB_DATA_DIMS];
  char metadataR[KLB_METADATA_SIZE];
  enum KLB_DATA_TYPE dataTypeR;
  enum KLB_COMPRESSION_TYPE compressionTypeR;
  float32_t pixelSizeR[KLB_METADATA_SIZE];
  void *imRead = (void*)NULL;

  ImageIO_close(im);
  im->fd = NULL;

  /* -1 is numThreads
   * numThreads <= 0 means as many threads as possible
   */
  imRead = (void*)readKLBstack( name, xyzctR, &dataTypeR, -1, pixelSizeR, blockSizeR, &compressionTypeR, metadataR );
  if ( imRead == (void*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading file '%s'\n", proc, name );
    return( -1 );
  }

  if ( _debug_ ) {
    fprintf( stderr, "dimensions = %d x %d x %d x %d x %d\n", xyzctR[0], xyzctR[1], xyzctR[2], xyzctR[3], xyzctR[4] );
    fprintf( stderr, "block size = %d x %d x %d x %d x %d\n", blockSizeR[0], blockSizeR[1], blockSizeR[2], blockSizeR[3], blockSizeR[4] );
    fprintf( stderr, "compression type = " );
    switch( compressionTypeR ) {
    default :
      fprintf( stderr, "unknown\n" ); break;
    case NONE :
      fprintf( stderr, "NONE\n" ); break;
    case BZIP2 :
      fprintf( stderr, "BZIP2\n" ); break;
    /* case ZLIB : */
    case 2 :
      fprintf( stderr, "ZLIB\n" ); break;
    }
  }

  /* data type
   */
  switch( dataTypeR ) {
  default :
      free( imRead );
      if ( _verbose_ )
        fprintf( stderr, "%s: such datatype (%d) not handled yet\n",
                 proc, dataTypeR );
      return( -1 );
  case INT8_TYPE :
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case UINT8_TYPE :
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case INT16_TYPE :
      im->wdim = 2;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case UINT16_TYPE :
      im->wdim = 2;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case INT32_TYPE :
      im->wdim = 4;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case UINT32_TYPE :
      im->wdim = 4;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case INT64_TYPE :
      im->wdim = 8;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case UINT64_TYPE :
      im->wdim = 8;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case FLOAT32_TYPE :
      im->wdim = 4;
      im->wordKind = WK_FLOAT;
      im->sign = SGN_SIGNED;
      break;
  case FLOAT64_TYPE :
      im->wdim = 8;
      im->wordKind = WK_FLOAT;
      im->sign = SGN_SIGNED;
      break;
  }

  /* dimensions
   */
  im->xdim = xyzctR[0];
  im->ydim = xyzctR[1];
  im->zdim = xyzctR[2];
  if ( xyzctR[3] > 1 && xyzctR[4] > 1 ) {
      free( imRead );
      if ( _verbose_ )
        fprintf( stderr, "%s: can not handle both c>1 and t>1\n",
                 proc );
      return( -1 );
  }
  if ( xyzctR[3] > 1 )
    im->vdim = xyzctR[3];
  if ( xyzctR[4] > 1 )
    im->vdim = xyzctR[4];

  /* voxel size
   */
  im->vx = pixelSizeR[0];
  im->vy = pixelSizeR[1];
  im->vz = pixelSizeR[2];

  im->data = (unsigned char *) imRead;

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/

int writeKlbImage( char *name, _image *im )
{
  char *proc = "writeKlbImage";
  uint32_t	xyzctR[KLB_DATA_DIMS];
  enum KLB_DATA_TYPE dataTypeR;
  enum KLB_COMPRESSION_TYPE compressionTypeR = BZIP2;
  float32_t pixelSizeR[KLB_METADATA_SIZE];
  int error;

  switch ( im->wordKind ) {
  case WK_FIXED :
    switch( im->sign ) {
    case SGN_UNSIGNED :
      if ( im->wdim  == sizeof( unsigned char ) ) {
          dataTypeR = UINT8_TYPE;
      }
      else if ( im->wdim  == sizeof( unsigned short int ) ) {
          dataTypeR = UINT16_TYPE;
      }
      else if ( im->wdim  == sizeof( unsigned int ) ) {
          dataTypeR = UINT32_TYPE;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown WK_FIXED UNSIGNED word dim (%d)\n ", proc, im->wdim );
        return( -1 );
      }
      break;
    case SGN_SIGNED :
      if ( im->wdim  == sizeof( char ) ) {
          dataTypeR = INT8_TYPE;
      }
      else if ( im->wdim  == sizeof( short int ) ) {
          dataTypeR = INT16_TYPE;
      }
      else if ( im->wdim  == sizeof( int ) ) {
          dataTypeR = INT32_TYPE;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown WK_FIXED SIGNED word dim (%d)\n ", proc, im->wdim );
        return( -1 );
      }
      break;
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown wordSign\n ", proc );
      return( -1 );
    }
    break;
  case WK_FLOAT :
    if ( im->wdim  == sizeof( float ) ) {
        dataTypeR = FLOAT32_TYPE;
    }
    else if ( im->wdim  == sizeof( double ) ) {
        dataTypeR = FLOAT64_TYPE;
    }
    else {
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown WK_FLOAT word dim\n ", proc );
      return( -1 );
    }
    break;
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: unknown wordKind\n ", proc );
    return( -1 );
  }

  xyzctR[0] = im->xdim;
  xyzctR[1] = im->ydim;
  xyzctR[2] = im->zdim;
  if ( im->vdim == 3 && dataTypeR == UINT8_TYPE ) {
    xyzctR[3] = im->vdim;
    xyzctR[4] = 1;
  }
  else if ( im->vdim > 1 ) {
    xyzctR[3] = 1;
    xyzctR[4] = im->vdim;
  }
  else {
    xyzctR[3] = 1;
    xyzctR[4] = 1;
  }

  pixelSizeR[0] = im->vx;
  pixelSizeR[1] = im->vy;
  pixelSizeR[2] = im->vz;
  pixelSizeR[3] = 1;
  pixelSizeR[4] = 1;

  error = writeKLBstack( (const void*)im->data, name, xyzctR, dataTypeR, -1, pixelSizeR, NULL, compressionTypeR, NULL);
  if ( error > 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when wrting file '%s'\n ", proc, name );
      return( -1 );
  }


  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



PTRIMAGE_FORMAT createKlbFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat = &testKlbHeader;
  f->readImageHeader = &readKlbImage;
  f->readImageData = NULL;
  f->writeImageHeader = NULL;
  f->writeImage = &writeKlbImage;

  strcpy(f->fileExtension,".klb");
  strcpy(f->realName,"Klb");
  return f;
}








