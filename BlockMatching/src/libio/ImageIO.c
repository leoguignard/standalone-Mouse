
#ifdef WIN32
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <stdio.h>
#endif

#include <stdio.h>
int fileno(FILE *stream);
#include <stdlib.h>
#include <string.h>
char *strdup(const char *s);
#include <limits.h>
#include <math.h>

#include <locale.h>

#include <ImageIO.h>

/* formats actuellement lus

   format     |   extension(s)   |  lecture  | ecriture
   INRIMAGE   |  .inr[.gz]       |     X     |     X     -> + .gradient[.gz] + .gradient_direction[.gz]
   GIS        |  .dim, .ima[.gz] |     X     |     X
   ANALYZE    |  .hdr, .img[.gz] |     X     |     X
   PNM        |  .ppm, .pgm      |     X     |     X
   GIF        |  .gif            |     X     |
   BMP        |  .bmp            |     X     |
   RAW        |  .raw            |           |     X
*/

#include "inr.h"
#include "gif.h"
#include "pnm.h"
#include "bmp.h"
#include "iris.h"
#ifdef _NIFTICLIB_
#include "nifti.h"
#endif
#ifdef _KLB_
#include "klb.h"
#endif
#ifdef _LIBLIBTIFF_
#include "tif.h"
#endif
#include "analyze.h"
#include "raw.h"
#include "nachos.h"
#include <metaImage.h>

#ifdef MINC_FILES
#  include "mincio.h"
#endif





static int _ImageIO_debug_ = 0;
static int _ImageIO_verbose_ = 0;



void _SetVerboseInImageIO( int v )
{
  _ImageIO_verbose_ = v;
}

void _IncrementVerboseInImageIO(  )
{
  _ImageIO_verbose_ ++;
}

void _SetDebugInImageIO( int d )
{
  _ImageIO_debug_ = d;
}

void _IncrementDebugInImageIO(  )
{
  _ImageIO_debug_ ++;
}





static char LC_NUMERIC_LOCALE[128];
static int LC_NUMERIC_LOCALE_is_set = 0;

static void _set_LC_NUMERIC_LOCALE_to_C()
{
  char *proc = "_get_LC_NUMERIC_LOCALE";
  int i;

  if ( _ImageIO_debug_ >= 1 ) {
    fprintf( stderr, "%s (entering): LC_NUMERIC_LOCALE is", proc );
    if ( LC_NUMERIC_LOCALE_is_set == 0 ) {
      fprintf( stderr, " not set\n" );
    }
    else {
      fprintf( stderr, " set to '%s'\n", LC_NUMERIC_LOCALE );
    }
  }

  if ( LC_NUMERIC_LOCALE_is_set == 0 ) {

    for ( i=0; i<128; i++ ) LC_NUMERIC_LOCALE[i] = '\0';
    (void)strcpy( LC_NUMERIC_LOCALE, setlocale( LC_NUMERIC, (char*)NULL ) );
    LC_NUMERIC_LOCALE_is_set = 1;

    if ( _ImageIO_debug_ >= 1 ) {
      fprintf( stderr, "%s: LC_NUMERIC is '%s'\n", proc, setlocale( LC_NUMERIC, (char*)NULL ) );
    }

    if ( LC_NUMERIC_LOCALE[0] != 'C' ) {
      if ( setlocale( LC_NUMERIC, "C" ) == (char*)NULL ) {
        if ( _ImageIO_verbose_ || _ImageIO_debug_ )
          fprintf( stderr, "%s: unable to change LC_NUMERIC locale to 'C'\n", proc );
      }
    }
  }

  if ( _ImageIO_debug_ >= 1 ) {
    fprintf( stderr, "%s (exiting): LC_NUMERIC is '%s'\n", proc, setlocale( LC_NUMERIC, (char*)NULL ) );
  }

}

static void _set_LC_NUMERIC_LOCALE_back()
{
  char *proc = "_set_LC_NUMERIC_LOCALE_back";
  int i;

  if ( _ImageIO_debug_ >= 1 ) {
    fprintf( stderr, "%s (entering): LC_NUMERIC_LOCALE is", proc );
    if ( LC_NUMERIC_LOCALE_is_set == 0 ) {
      fprintf( stderr, " not set\n" );
    }
    else {
      fprintf( stderr, " set to '%s'\n", LC_NUMERIC_LOCALE );
    }
  }

  if ( LC_NUMERIC_LOCALE_is_set ) {

    if ( setlocale( LC_NUMERIC, LC_NUMERIC_LOCALE ) == (char*)NULL ) {
      if ( _ImageIO_verbose_ || _ImageIO_debug_ )
        fprintf( stderr, "%s: unable to change back LC_NUMERIC locale to its original value '%s'\n", proc, LC_NUMERIC_LOCALE );
    }
    else {
      if ( _ImageIO_debug_ >= 1 ) {
        fprintf( stderr, "%s: LC_NUMERIC is set back to '%s'\n", proc, setlocale( LC_NUMERIC, (char*)NULL ) );
      }
    }

    for ( i=0; i<128; i++ ) LC_NUMERIC_LOCALE[i] = '\0';
    LC_NUMERIC_LOCALE_is_set = 0;
  }

  if ( _ImageIO_debug_ >= 1 ) {
    fprintf( stderr, "%s (exiting): LC_NUMERIC is '%s'\n", proc, setlocale( LC_NUMERIC, (char*)NULL ) );
  }

}



static void _nifti_mat44_to_quatern( _image *im );
static void _nifti_quatern_to_mat44( _image *im );

/*--------------------------------------------------
 *
 * image reader/writer generic format management
 *
 --------------------------------------------------*/



/** the first file format is initialized to null */
static PTRIMAGE_FORMAT firstFormat=NULL;

/** the Inrimage file format (default format) is initialized to null */
static PTRIMAGE_FORMAT InrimageFormat=NULL;


void _initImageFormat( PTRIMAGE_FORMAT format )
{
  int i;
  format->testImageFormat = NULL;

  format->readImageHeader = NULL;
  format->readImageData = NULL;

  format->writeImageHeader = NULL;
  format->writeImage = NULL;

  for ( i=0; i<IMAGE_FORMAT_NAME_LENGTH; i++ )
    format->fileExtension[i] = '\0';
  for ( i=0; i<IMAGE_FORMAT_NAME_LENGTH; i++ )
    format->realName[i] = '\0';

  format->next = NULL;
}



/** adds a format at the beginning of the list of image formats.
    Test if all mandatory fields have been filled */
static int addImageFormat( PTRIMAGE_FORMAT format)
{

  if ( ! format->testImageFormat ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr,"addImageFormat: no test format procedure in %s\n",
        format->realName);
    return -1;
  }

  if ( format->writeImage ||
       ( (format->readImageHeader)
         && (strlen(format->fileExtension)>0)
         && (strlen(format->realName)>0) ) ) {

    format->next=firstFormat;
    firstFormat=format;

    return 0;

  }
  else {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr,"addImageFormat: information missing in file format %s\n",
              format->realName);
    return -1;
  }
}



/** adds a format at the end of the list of image formats.
    Test if all mandatory fields have been filled */
static int addImageFormatAtEnd( PTRIMAGE_FORMAT format)
{
  PTRIMAGE_FORMAT f;

  if ( ! format->testImageFormat ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr,"addImageFormatAtEnd: no test format procedure in %s\n",
        format->realName);
    return -1;
  }

  if ( format->writeImage ||
       ( (format->readImageHeader)
         && (strlen(format->fileExtension)>0)
         && (strlen(format->realName)>0) ) ) {

    format->next = NULL;

    if (firstFormat == NULL) {
      firstFormat=format;
    }
    else {
      for(f=firstFormat;(f->next!=NULL);f=f->next)
        ;
      f->next=format;
    }

    return 0;

  }
  else {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr,"addImageFormatAtEnd: information missing in file format %s\n",
              format->realName);
    return -1;
  }
}



/** creates supported image formats */
static void initSupportedFileFormat()
{
  PTRIMAGE_FORMAT f;
  if ( InrimageFormat == NULL ) {

#ifdef _NIFTICLIB_
    f = createNiftiFormat();
    addImageFormatAtEnd( f );
#else
    f = createAnalyzeFormat();
    addImageFormatAtEnd( f );
#endif

#ifdef _KLB_
    f = createKlbFormat();
    addImageFormatAtEnd( f );
#endif

#ifdef _LIBLIBTIFF_
    f = createTifFormat();
    addImageFormatAtEnd( f );
#endif


    f = createBMPFormat();
    addImageFormatAtEnd( f );

    f = createGifFormat();
    addImageFormatAtEnd( f );

    f = createIrisFormat();
    addImageFormatAtEnd( f );

    f = createPnmFormat();
    addImageFormatAtEnd( f );

    f = createRawFormat();
    addImageFormatAtEnd( f );

    f = createMetaImageFormat();
    addImageFormatAtEnd( f );

    f = createNachosFormat();
    addImageFormatAtEnd( f );

    InrimageFormat = createInrimageFormat();
    addImageFormat( InrimageFormat );
  }
}



#ifdef _UNUSED_

PTRIMAGE_FORMAT firstImageFormat() {
  return firstFormat;
}

#endif



/** prints supported image formats */
void printSupportedFileFormat() {
  PTRIMAGE_FORMAT f;
  int i;

  initSupportedFileFormat();

  for(i=0, f=firstFormat;(f!=NULL);i++, f=f->next) {
    if ( (f->testImageFormat) &&
         (f->readImageHeader) &&
         (strlen(f->fileExtension)>0) &&
         (strlen(f->realName)>0)) {
      fprintf( stderr, "#%2d: format name ='%s', extensions='%s'",
              i, f->realName, f->fileExtension );
      if (f->readImageHeader)
        fprintf( stderr, ", read" );
      if (f->writeImage)
        fprintf( stderr, ", write" );
      fprintf( stderr, "\n" );
    }
  }
}



#ifdef _UNUSED_

/** remove supported image formats */
void removeSupportedFileFormat() {
  PTRIMAGE_FORMAT f=firstFormat;

  while( f != NULL) {
    PTRIMAGE_FORMAT f_old = f;
    f = f->next;
    ImageIO_free( f_old);
  }
  InrimageFormat=NULL;

}

#endif


#ifdef _UNUSED_

/* check the type of image in fileName
 * with respect to magic numbers
 */
PTRIMAGE_FORMAT imageType(const char *fileName) {
  _ImageIO_file f;
  char magic[5];
  PTRIMAGE_FORMAT format;

  if(!fileName) {
#ifdef ZLIB
    f = gzdopen(fileno(stdin), "rb");
#else
    f = fdopen(fileno(stdin), "rb");
#endif
  }
  else {
#ifdef ZLIB
    f = gzopen(fileName, "rb");
#else
    f = fopen(fileName, "rb");
#endif
  }

  if(!f) return NULL;

#ifdef ZLIB
  gzread( f, (void *) magic, 4);
#else
  fread( (void *) magic, 1, 4, f );
#endif


  magic[4] = '\0';

#ifdef ZLIB
  gzclose( f );
#else
  if(fileName) fclose( f );
#endif

  if (firstFormat==NULL)
    initSupportedFileFormat();

  for(format=firstFormat;(format!=NULL);format=format->next) {
    /* test if it is the correct header based on magic and file extension */
    if (((*format->testImageFormat)(magic,fileName)) >=0) {
      return format;
    }
  }
  return 0;

}

#endif





/*--------------------------------------------------
 *
 * find file format
 *
 --------------------------------------------------*/


/* get the adequate reader from the open file
 * test magic strings
 */
static PTRIMAGE_FORMAT _getImageFormatFromMagicString( _image *im, const char *name )
{
  char *proc = "_getImageFormatFromMagicString";
  char magic[5];
  PTRIMAGE_FORMAT f;

  /* standard input
   */
  if ( im->openMode == OM_STD) {
      return( InrimageFormat );
  }

  /* standard file
   * get magic string for disk files
   */
  ImageIO_read(im, magic, 4);
  magic[4] = '\0';
  ImageIO_seek(im, 0L, SEEK_SET);
  /** test each format */
  for(f=firstFormat;(f!=NULL)&& (im->imageFormat==NULL);f=f->next) {
      /* test whether 'testImageFormat' exists
       */
      if ( ! f->testImageFormat )
          continue;
      /* test if it is the correct format based on magic and file extension
       */
      if (((*f->testImageFormat)(magic, name)) >=0) {
        if ( _ImageIO_debug_ )
          fprintf( stderr, "%s: format '%s' recognized for '%s'\n",
                   proc, f->realName, name );
        return( f );
    }
  }

  return( NULL );
}



static PTRIMAGE_FORMAT _getImageFormatFromName( const char *name )
{
  char *proc = "_getImageFormatFromName";
  int i, length, extLength;
  PTRIMAGE_FORMAT f;
  char ext[IMAGE_FORMAT_NAME_LENGTH];
  char *ptr;

  /* different conventions for the standard input/output
   */
  if ( name == (char*)NULL || name[0] == '\0'
       || (name[0] == '-' && name[1] == '\0')
       || (name[0] == '>' && name[1] == '\0')
       || (name[0] == '<' && name[1] == '\0') ) {
    return( InrimageFormat );
  }

  /* scan all formats; */
  length=strlen( name );

  for( f=firstFormat; f!=NULL; f=f->next ) {
    /* scan all extensions for that format */
    ptr=&f->fileExtension[0];

    do {
      /* get next file extension */
      for(i=0;((*ptr)!=',' && (*ptr)!='\0');i++,ptr++) {
        ext[i]=(*ptr);
      }
      if ((*ptr)==',') {
        ext[i]='\0';
        ptr++;
      }
      else {
        ext[i]='\0';
      }
      extLength=strlen(ext);

      /* test if the tail of name matches the extension */
      if ( (length > extLength) && (!strcmp( name + length - extLength, ext)) ) {
        if ( _ImageIO_debug_ )
            fprintf( stderr, "%s: format '%s' recognized for '%s'\n",
                     proc, f->realName, name );
        return( f );
      }

    } while ( (*ptr)!='\0' );
  }

  return( InrimageFormat );
  /* default is Inrimage !
   * return( NULL );
   */
}



static PTRIMAGE_FORMAT _getImageFormatForReading( _image *im, const char *name )
{
  PTRIMAGE_FORMAT f;

  f = _getImageFormatFromMagicString( im, name );
  if ( f != NULL ) {
    return( f );
  }

  f = _getImageFormatFromName( name );
  if ( f != NULL ) {
    return( f );
  }

  return( NULL );
}




/*--------------------------------------------------
 *
 * endianness stuff
 *
 --------------------------------------------------*/


ENDIANNESS _getEndianness()
{
  union {
    unsigned char uc[2];
    unsigned short us;
  } twobytes;
  twobytes.us = 255;
  /* on linux or dec
   */
  if ( twobytes.uc[1] == 0 ) return( END_LITTLE );
  /* on solaris or sgi
   */
  return( END_BIG );
}





/*--------------------------------------------------
 *
 * image definition
 *
 --------------------------------------------------*/



/* Allocates and initializes an image descriptor */
_image *_initImage() {
  _image *im;
  int i, j;

  im = (_image *) ImageIO_alloc(sizeof(_image));
  if ( im == NULL ) return( im );

  /* default image size is 1*1*1 */
  im->xdim = im->ydim = im->zdim = im->vdim = 1;
  /* default image voxel size is 1.0*1.0*1.0 */
  im->vx = im->vy = im->vz = 1.0;

  /* default image offset  is 0 0 0 */
  im->t_is_set = 0;
  im->tx = im->ty = im->tz = 0.0;


  /* default image rotation  is 0 0 0 */
  im->q_is_set = 0;
  im->qb = im->qc = im->qd = 0.0;

  /* default image qfac  is 1.0 */
  im->qfac_is_set = 0;
  im->qfac = 1.0;

  /* qform_code and sform_code
   */
  im->qform_code = 0;
  im->sform_code = 0;

  /* default qform_toreal matrix
   */
  for ( i=0; i<4; i++ )
  for ( j=0; j<4; j++ )
    im->qform_toreal[i][j] = 0.0;
  for ( i=0; i<4; i++ )
    im->qform_toreal[i][i] = 1.0;

  /* default sform_toreal matrix
   */
  for ( i=0; i<4; i++ )
  for ( j=0; j<4; j++ )
    im->sform_toreal[i][j] = 0.0;
  for ( i=0; i<4; i++ )
    im->sform_toreal[i][i] = 1.0;


  /* no data yet */
  im->data = NULL;

  /* no file associated to image */
  im->fd = NULL;
  im->openMode = OM_CLOSE;
  im->endianness = _getEndianness();

  /* unknown data kind
     default is binary
   */
  im->dataMode = DM_BINARY;

  /* no user string */
  im->user = NULL;
  im->nuser = 0;

  /* unknown word kind */
  im->wdim = 0;
  im->wordKind = WK_UNKNOWN;
  im->vectMode = VM_SCALAR;
  im->sign = SGN_UNKNOWN;
  im->imageFormat = NULL;

  /** eventually initializes the supported file formats */
  if (firstFormat==NULL)
    initSupportedFileFormat();
  /* return image descriptor */
  return im;
}



/* Free an image descriptor */
void _freeImageStructure( _image *im )
{
  unsigned int i;

  if ( im == NULL ) return;

  /* close image if opened */
  if(im->openMode != OM_CLOSE) ImageIO_close(im);

  /* free user string array if any */
  if( (im->nuser > 0) && (im->user != NULL) ) {
    for(i = 0; i < im->nuser; i++)
      if ( im->user[i] != NULL ) ImageIO_free(im->user[i]);
    ImageIO_free(im->user);
  }
  im->nuser = 0;
  im->user = NULL;

  /* free given descriptor */
  ImageIO_free( im );
}



/* Free an image descriptor */
void _freeImage( _image *im )
{
  if ( im == NULL ) return;

  /* free data if any */
  if(im->data != NULL) ImageIO_free(im->data);
  im->data = NULL;

  _freeImageStructure( im );
}



void _swapImageData( _image *im )
{
  char *proc = "_swapImageData";
  unsigned char *ptr1, *ptr2, b[8];
  unsigned short int si, *ptr3, *ptr4;
  unsigned int        i, *ptr5, *ptr6;
  unsigned int long size, length;

  if( _getEndianness() != im->endianness) {

    if ( _ImageIO_debug_ >= 2 )
      fprintf( stderr, "%s: swap image data\n", proc );

    size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;
    if ( size <= 0 ) return;

    length = size / im->wdim;
    ptr1 = ptr2 = (unsigned char *) im->data;

    /* 2 bytes swap */
    if(im->wdim == 2) {
      /*
        while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        *ptr2++ = b[1];
        *ptr2++ = b[0];
        }
      */
      ptr3 = ptr4 = (unsigned short int *) im->data;
      while( length-- ) {
        si = *ptr3++;
        *ptr4++ = ((si >> 8) & 0xff) | (si << 8);
      }
    }

    /* 4 bytes swap */
    else if(im->wdim == 4) {
      /*
        while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        b[2] = *ptr1++;
        b[3] = *ptr1++;
        *ptr2++ = b[3];
        *ptr2++ = b[2];
        *ptr2++ = b[1];
        *ptr2++ = b[0];
        }
      */
      ptr5 = ptr6 = (unsigned int *) im->data;
      while( length-- ) {
        i = *ptr5++;
        *ptr6++ =  (i << 24) | ((i & 0xff00) << 8) | ((i >> 8) & 0xff00) | ((i >> 24) & 0xff);
      }
    }
    /* 8 bytes swap */
    else if(im->wdim == 8) {
      while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        b[2] = *ptr1++;
        b[3] = *ptr1++;
        b[4] = *ptr1++;
        b[5] = *ptr1++;
        b[6] = *ptr1++;
        b[7] = *ptr1++;
        *ptr2++ = b[7];
        *ptr2++ = b[6];
        *ptr2++ = b[5];
        *ptr2++ = b[4];
        *ptr2++ = b[3];
        *ptr2++ = b[2];
        *ptr2++ = b[1];
        *ptr2++ = b[0];
      }
    }
  }
}





/*--------------------------------------------------
 *
 * mimics standard routines
 *
 --------------------------------------------------*/


/* default allocation routine */
static void *(*allocRoutine)(size_t) = 0;
/* default deallocation routine */
static void (*deleteRoutine)(void *) = 0;



/* set allocation and deallocation routines */
void setImageIOAllocationRoutines(ALLOCATION_FUNCTION alloc,
                                  DEALLOCATION_FUNCTION del) {
  if(alloc != NULL) allocRoutine = alloc;
  if(del != NULL) deleteRoutine = del;
}

void *ImageIO_alloc(size_t s) {
  if(!allocRoutine) allocRoutine = malloc;
  return ( (*allocRoutine)(s) );
}
/* call allocation routine */
void ImageIO_free(void *m) {
  if(!deleteRoutine) deleteRoutine = free;
  (*deleteRoutine)(m);
}
/* call deallocation routine */



/* mimics fwrite() function.
   According to _openWriteImage(), openMode will has one
   of the following value:
   - OM_STD (for stdout)
   - OM_GZ
   - OM_FILE
*/
size_t ImageIO_write(const _image *im, const void *buf, size_t len)
{
  char *proc = "ImageIO_write";
  size_t to_be_written = len;
  size_t stepWriteMax = INT_MAX - 1;
  size_t l;
  char *b = (char*)buf;

  if ( buf == (void*)NULL ) {
    if ( _ImageIO_verbose_ )
      fprintf( stderr, "%s: NULL buffer\n", proc );
    return( 0 );
  }

  size_t stepWrite = (stepWriteMax < to_be_written) ? stepWriteMax : to_be_written;

  if ( _ImageIO_debug_ >= 2 )
    fprintf( stderr, "%s: length to be written = %lu ; stepWrite =  %lu bytes ; \n", proc, to_be_written, stepWrite );

  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return 0;
  case OM_STD :
#ifdef ZLIB
    while ( (to_be_written > 0) && ((l = gzwrite(im->fd, (void *) b, stepWrite)) > 0) ) {
      to_be_written -= l;
      stepWrite = (stepWriteMax < to_be_written) ? stepWriteMax : to_be_written;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: gzwrite has written %lu, remains %lu, step=%lu\n", proc, l, to_be_written, stepWrite );
      b += l;
    }
#else
    while ( (to_be_written > 0) && ((l = fwrite(b, 1, stepWrite, (FILE*)im->fd)) > 0) ) {
      to_be_written -= l;
      stepWrite = (stepWriteMax < to_be_written) ? stepWriteMax : to_be_written;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: fwrite has written %lu, remains %lu, step=%lu\n", proc, l, to_be_written, stepWrite );
      b += l;
    }
#endif
    return ( len - to_be_written );
    break;
#ifdef ZLIB
  case OM_GZ :
    while ( (to_be_written > 0) && ((l = gzwrite(im->fd, (void *) b, stepWrite)) > 0) ) {
      to_be_written -= l;
      stepWrite = (stepWriteMax < to_be_written) ? stepWriteMax : to_be_written;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: gzwrite has written %lu, remains %lu, step=%lu\n", proc, l, to_be_written, stepWrite );
      b += l;
    }
    return ( len - to_be_written );
    break;
#endif
  case OM_FILE:
    while ( (to_be_written > 0) && ((l = fwrite(b, 1, stepWrite, (FILE*)im->fd)) > 0) ) {
      to_be_written -= l;
      stepWrite = (stepWriteMax < to_be_written) ? stepWriteMax : to_be_written;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: fwrite has written %lu, remains %lu, step=%lu\n", proc, l, to_be_written, stepWrite );
      b += l;
    }
    return ( len - to_be_written );
    break;
  }
  return 0;
}



/* mimics fread() function.
   According to _openReadImage(), openMode will has one
   of the following value:
   - OM_STD (for stdin)
   - OM_GZ *or* OM_FILE
*/
size_t ImageIO_read(const _image *im, void *buf, size_t len)
{
  char *proc = "ImageIO_read";
  size_t to_be_read = len;
  size_t stepReadMax = INT_MAX - 1;
#ifdef ZLIB
  int lgzread;
#else
  size_t lfread;
#endif
  char *b = (char*)buf;

  size_t stepRead = (stepReadMax < to_be_read) ? stepReadMax : to_be_read;

  if ( _ImageIO_debug_ >= 2 )
    fprintf( stderr, "%s: length to be read = %lu ; stepRead =  %lu bytes ; \n", proc, to_be_read, stepRead );

  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return 0;
  case OM_STD :
#ifdef ZLIB
    while ( (to_be_read > 0) && ((lgzread = gzread(im->fd, (void *) b, stepRead)) > 0) ) {
      to_be_read -= lgzread;
      stepRead = (stepReadMax < to_be_read) ? stepReadMax : to_be_read;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: gzread has read %d, remains %lu, step=%lu\n", proc, lgzread, to_be_read, stepRead );
      b += lgzread;
    }
#else
    while ( (to_be_read > 0) && ((lfread = fread( b, 1, stepRead, im->fd )) > 0) ) {
      to_be_read -= lfread;
      stepRead = (stepReadMax < to_be_read) ? stepReadMax : to_be_read;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: fread has read %lu, remains %lu, step=%lu\n", proc, lfread, to_be_read, stepRead );
      b += lfread;
    }
#endif
    return ( len - to_be_read );
    break;
#ifdef ZLIB
  case OM_GZ :
    while ( (to_be_read > 0) && ((lgzread = gzread(im->fd, (void *) b, stepRead)) > 0) ) {
      to_be_read -= lgzread;
      stepRead = (stepReadMax < to_be_read) ? stepReadMax : to_be_read;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: gzread has read %d, remains %lu, step=%lu\n", proc, lgzread, to_be_read, stepRead );
      b += lgzread;
    }
    return ( len - to_be_read );
    break;
#else
  case OM_FILE :
    while ( (to_be_read > 0) && ((lfread = fread( b, 1, stepRead, im->fd )) > 0) ) {
      to_be_read -= lfread;
      stepRead = (stepReadMax < to_be_read) ? stepReadMax : to_be_read;
      if ( _ImageIO_debug_ >= 2 )
        fprintf( stderr, "%s: fread has read %lu, remains %lu, step=%lu\n", proc, lfread, to_be_read, stepRead );
      b += lfread;
    }
    return ( len - to_be_read );
    break;
#endif
  }

  return 0;
}



/* mimics fgets() function.
   According to _openReadImage(), openMode will has one
   of the following value:
   - OM_STD (for stdout)
   - OM_GZ *or* OM_FILE
*/
char *ImageIO_gets( const _image *im, char *str, int size )
{
  char *ret = NULL;
  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return NULL;
  case OM_STD :
#ifdef ZLIB
    ret = (char *) gzgets(im->fd, str, size );
#else
    ret = fgets(str, size, im->fd);
#endif
    break;
#ifdef ZLIB
  case OM_GZ :
    ret = (char *) gzgets(im->fd, str, size);
    break;
#else
  case OM_FILE :
    ret = fgets(str, size, im->fd);
    break;
#endif
  }
  return ret;
}



int ImageIO_seek( const _image *im, long offset, int whence ) {
  switch(im->openMode) {
  case OM_CLOSE :
  default :
    return -1;
#ifdef ZLIB
  case OM_GZ:
    return gzseek(im->fd, offset, whence );
#endif
  case OM_FILE:
    return fseek( (FILE*)im->fd, offset, whence );
  }
}

/* return non 0 in case of error
 */
int ImageIO_error( const _image *im )
{
  static int errnum;
  switch(im->openMode) {
  case OM_CLOSE :
  default :
    return 0;
#ifdef ZLIB
  case OM_GZ :
    (void)gzerror(im->fd, &errnum);
    return( (errnum != Z_OK) || gzeof(im->fd) );
#endif
  case OM_FILE :
    return( ferror( (FILE*)im->fd ) || feof( (FILE*)im->fd ) );
  }
  return 0;
}



/* Upon successful completion 0 is returned.

   Closing the standard output with gzclose()
   is necessary as it will flush the pending output.
 */
int ImageIO_close( _image* im )
{
  int ret=0;

  switch ( im->openMode ) {
  default :
  case OM_CLOSE :
    break;
#ifdef ZLIB
  case OM_GZ :
  case OM_STD :
    ret = gzclose( im->fd );
    break;
#else
  case OM_STD :
    break;
#endif
  case OM_FILE :
    ret = fclose( (FILE*)im->fd );
  }
  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return ret;
}





/*--------------------------------------------------
 *
 * image I/O
 *
 --------------------------------------------------*/



/* given an initialized file descriptor and a file name,
   open file from stdin (if name == NULL, or name == "-", or name == "<"),
   or a standard/gzipped file otherwise (gzipped files are handled assuming
   that it is compiled and linked with zlib).
   openMode will have one of the following value:
   - OM_STD (for stdin)
   - OM_GZ *or* OM_FILE
*/
void _openReadImage(_image* im, const char *name) {
  if(im->openMode == OM_CLOSE) {

    /* open from stdin */
    if( name == NULL || name[0] == '\0'
        || (name[0] == '-' && name[1] == '\0')
        || (name[0] == '<' && name[1] == '\0') ) {
#ifdef ZLIB
      im->fd = gzdopen(fileno(stdin), "rb");
#else
      im->fd = fdopen(fileno(stdin), "rb");
#endif
      im->openMode = OM_STD;
    }

    else {
#ifdef ZLIB
      im->fd = gzopen(name, "rb");
      if(im->fd) im->openMode = OM_GZ;
#else
      im->fd = fopen(name, "rb");
      if(im->fd) im->openMode = OM_FILE;
#endif

    }

  }
}





/* given an initialized file descriptor and a file name,
   open file from stdout (if name == NULL, or name == "-", or name == ">"),
   a gzipped pipe (if name got the extension ".gz")
   or a standard file otherwise.
   openMode will have one of the following value:
   - OM_STD (for stdout)
   - OM_GZ
   - OM_FILE
*/
void _openWriteImage(_image* im, const char *name)
{
  im->openMode = OM_CLOSE;

  if( name == NULL || name[0] == '\0'
      || (name[0] == '-' && name[1] == '\0')
      || (name[0] == '>' && name[1] == '\0') ) {

#ifdef ZLIB
#if (defined _LINUX_) || (defined _SOLARIS_) || (defined _SGI_)
    im->fd = gzdopen(1, "wb");
#else
    im->fd = gzdopen(fileno(stdout), "wb");
#endif
#else
    im->fd = (_ImageIO_file) stdout;
#endif
    im->openMode = OM_STD;
  }

  else{
#ifdef ZLIB

    /* from gzopen() doc:
       ... The mode parameter is as in fopen ("rb" or "wb") but can
       also include a compression level ("wb9") or a strategy: 'f' for
       filtered data as in "wb6f", 'h' for Huffman only compression as
       in "wb1h" ...
       However, a small .gz header will be written ... thus gz(d)open can not
       be used for regular files.
    */

    if( !strncmp(name+strlen(name)-3, ".gz", 3) )
      {
#ifdef _MSC_VER
        int ffd=_open(name,_O_RDWR | _O_CREAT| _O_TRUNC | _O_BINARY, _S_IREAD|_S_IWRITE);
        im->fd = gzdopen( ffd, "wb" );
#else
        im->fd = gzopen( name, "wb" );
#endif
        im->openMode = OM_GZ;
      }
    else
#endif
    {
      im->fd = (_ImageIO_file) fopen(name, "wb");
      im->openMode = OM_FILE;
    }
  }
}






/*--------------------------------------------------
 *
 * image reading procedures
 *
 *--------------------------------------------------*/



static _image *_readImageHeaderAndGetError( const char *name_to_be_read, int *error )
{
  char *proc = "_readImageHeaderAndGetError";
  _image *im;
  char *name = NULL;
  int res;

  *error = ImageIO_NO_ERROR;

  _set_LC_NUMERIC_LOCALE_to_C();

  /* open image file */
  im = _initImage();
  if ( name_to_be_read == NULL || name_to_be_read[0] == '\0'
       || (name_to_be_read[0] == '-' && name_to_be_read[1] == '\0')
       || (name_to_be_read[0] == '<' && name_to_be_read[1] == '\0') ) {
    name = NULL;
  }
  else {
    name = strdup( name_to_be_read );
  }


  _openReadImage(im, name);

  if(!im->fd) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: error: unable to open file \'%s\'\n", proc, name );
    _freeImage(im);
    *error = ImageIO_OPENING;
    if ( name != NULL ) free( name );
    _set_LC_NUMERIC_LOCALE_back();
    return NULL;
  }

  /* get the image format reader
   */
  initSupportedFileFormat();
  im->imageFormat = _getImageFormatForReading( im, name );


  if ( im->imageFormat == NULL ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: does not find image format for \'%s\'\n", proc, name );
    ImageIO_close( im );
    _freeImage(im);
    *error = ImageIO_UNKNOWN_TYPE;
    if ( name != NULL ) free( name );
    _set_LC_NUMERIC_LOCALE_back();
    return NULL;
  }

  /* now tests if the header can be read correctly
   */

  if ( *(im->imageFormat)->readImageHeader == NULL ) {
      if ( _ImageIO_verbose_ || _ImageIO_debug_ )
        fprintf( stderr, "%s: no reader for image format '%s'\n",
                 proc, (im->imageFormat)->realName );
      ImageIO_close( im );
      _freeImage(im);
      *error = ImageIO_UNKNOWN_TYPE;
      if ( name != NULL ) free( name );
      _set_LC_NUMERIC_LOCALE_back();
      return NULL;
  }


  res=(*(im->imageFormat)->readImageHeader)(name,im);
  /* could read header only */
  if (res == 0) {
    if ( name != NULL ) free( name );
    _set_LC_NUMERIC_LOCALE_back();
    return( im );
  }
  /* could read header and data */
  else if ( res > 0 ) {
    ImageIO_close(im);
    if ( name != NULL ) free( name );
    _set_LC_NUMERIC_LOCALE_back();
    return im;
  }

  /* could not read error : throw error */
  if ( _ImageIO_verbose_ || _ImageIO_debug_ )
    fprintf(stderr, "%s: an error occurs when reading image\n", proc );
  if ( name == NULL || im->openMode == OM_STD) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr, "\t from \'standard input\'" );
  }
  else {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr, "\t from file \'%s\'", name );
  }
  if ( _ImageIO_verbose_ || _ImageIO_debug_ )
    fprintf(stderr, " using format \'%s\'\n", (im->imageFormat)->realName );
  ImageIO_close( im );
  _freeImage(im);
  *error = ImageIO_READING_HEADER;
  if ( name != NULL ) free( name );
  _set_LC_NUMERIC_LOCALE_back();
  return NULL;
}





/* read header from an image file

   if standard input, it's an inrimage
   if not, get a magic string
   and try to find the good format

   if data are in a separate file,
   the header reading procedure will open
   the data file.

   error:
   0  success
   -1 unknown image type
   -2 error while opening
   -3 error while reading header
   -4 error while reading header or data
 */
_image *_readImageHeader( const char *name ) {
  int error = 0;
  return( _readImageHeaderAndGetError( name, &error ) );
}





/* Read data of an inrimage.
   If im->data is not NULL, assume that the buffer was previously allocated
   Swap bytes depending on the endianness and the current architecture  */
int _readImageData(_image *im) {
  char *proc = "_readImageData";
  size_t size, nread;

  if(im->openMode != OM_CLOSE) {
    if ( _ImageIO_debug_ ) {
      fprintf( stderr, "%s: %lu * %lu * %lu * %d *%d\n", proc, im->xdim, im->ydim, im->zdim, im->vdim, im->wdim );
      }
    size = (size_t)im->xdim * (size_t)im->ydim * (size_t)im->zdim * (size_t)im->vdim * (size_t)im->wdim;

    if ( size <= 0 ) return -3;

    /* image buffer may have been allocated
     */
    if( !im->data ) {
      im->data = (unsigned char *) ImageIO_alloc(size);
      if ( _ImageIO_debug_ >= 2 ) {
        fprintf( stderr, "%s: allocate %lu bytes at %p\n", proc, size, im->data );
      }

      if( !im->data ) {
          if ( _ImageIO_verbose_ )
            fprintf( stderr, "%s: image buffer allocation failed\n", proc );
          return( -2 );
      }
    }

    nread = ImageIO_read(im, im->data, size);
    if ( nread != size ) {
        if ( _ImageIO_verbose_ )
          fprintf( stderr, "%s: image buffer reading failed\n", proc );
        return( -1 );
    }



    /* architecture is big endian and data little endian
       length = nb of points
     */
    _swapImageData( im );

  }

  return 1;
}





/* Reads an image from a file and returns an image descriptor or NULL if
   reading failed.
   Reads from stdin if image name is NULL. */
_image* _readImage(const char *name) {
  _image *im;

  _set_LC_NUMERIC_LOCALE_to_C();

  /* read header
   */
  im = _readImageHeader( name );

  /* set geometry
   * if qform_code is set,
   *   compute quaternion, translation and orientation
   * else if either quaternion, translation or orientation are set
   *   compute matrix
   *
   */
  if ( im != NULL ) {
    if ( im->qform_code ) {
      if ( ! im->t_is_set || ! im->q_is_set || ! im->qfac_is_set ) {
        _nifti_mat44_to_quatern( im );
        im->t_is_set = im->q_is_set = im->qfac_is_set = 1;
      }
    }
    else {
      if ( im->t_is_set || im->q_is_set || im->qfac_is_set ) {
        _nifti_quatern_to_mat44( im );
        im->qform_code = 1;
      }
    }
  }


  /* read data
   */
  if (im != NULL && im->openMode != OM_CLOSE) {
    /* read body */
    /* switch from generic reader to specific ones
       if(_readImageData(im) < 0) {
       GM Lun  6 mai 2013 12:19:28 CEST */

    if ( *(im->imageFormat)->readImageData
         && (*(im->imageFormat)->readImageData)(im) < 0 ) {

      if ( _ImageIO_verbose_ || _ImageIO_debug_ )
        fprintf(stderr, "_readImage: error: invalid data encountered in \'%s\'\n",
                name);

      _freeImage(im);
      _set_LC_NUMERIC_LOCALE_back();
      return NULL;
    }
    ImageIO_close(im);
  }

  _set_LC_NUMERIC_LOCALE_back();

  return im;
}



/*--------------------------------------------------
 *
 * image writing procedures
 *
 *--------------------------------------------------*/




int _writeImageHeader( _image *im, char *name_to_be_written )
{
  char *proc = "_writeImageHeader";
  int r = ImageIO_NO_ERROR;

  if ( im == NULL ) return -1;

  /* get the image format writer
   */
  initSupportedFileFormat();
  im->imageFormat = _getImageFormatFromName( name_to_be_written );

  if ( im->imageFormat == NULL ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: does not find image format for \'%s\'\n", proc, name_to_be_written );
    return( ImageIO_WRITING_HEADER );
  }

  if ( *(im->imageFormat)->writeImageHeader == NULL ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: no header writer for image format '%s'\n",
               proc, (im->imageFormat)->realName );
    return( ImageIO_WRITING_HEADER );
  }

  if ( (*im->imageFormat->writeImageHeader)(name_to_be_written, im) < 0 ) {
      if ( _ImageIO_verbose_ || _ImageIO_debug_ )
        fprintf(stderr, "%s: error: unable to write header \'%s\'\n",
                proc, name_to_be_written );
      return( ImageIO_WRITING_HEADER );
  }


  /* close file descriptor */
  ImageIO_close( im );

  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return r;
}






/* Write inrimage given in inr in file name. If file name's suffix is
   .gz, the image is gziped. If file name's suffix is .hdr, the image
   is written in ANALYZE format. If file name is NULL, image is written
   on stdout */
int _writeImage(_image *im, char *name_to_be_written )
{
  char *proc = "_writeImage";
  int r = ImageIO_NO_ERROR;

  if ( im == NULL ) return -1;

  _set_LC_NUMERIC_LOCALE_to_C();

  /* set geometry
   */
  if ( im->qform_code ) {
    if ( ! im->t_is_set || ! im->q_is_set || ! im->qfac_is_set ) {
      _nifti_mat44_to_quatern( im );
      im->t_is_set = im->q_is_set = im->qfac_is_set = 1;
    }
  }
  else {
    if ( im->t_is_set && (im->q_is_set || im->qfac_is_set) ) {
      _nifti_quatern_to_mat44( im );
      im->qform_code = 1;
    }
  }

  /* get the image format writer
   */
  initSupportedFileFormat();
  im->imageFormat = _getImageFormatFromName( name_to_be_written );

  if ( im->imageFormat == NULL ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: does not find image format for \'%s\'\n", proc, name_to_be_written );
    _set_LC_NUMERIC_LOCALE_back();
    return( ImageIO_WRITING_HEADER );
  }

  if ( *(im->imageFormat)->writeImage == NULL ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf( stderr, "%s: no image writer for image format '%s'\n",
               proc, (im->imageFormat)->realName );
    _set_LC_NUMERIC_LOCALE_back();
    return( ImageIO_WRITING_HEADER );
  }

  if ( (*im->imageFormat->writeImage)(name_to_be_written, im) < 0 ) {
    if ( _ImageIO_verbose_ || _ImageIO_debug_ )
      fprintf(stderr, "%s: error: unable to write \'%s\'\n",
              proc, name_to_be_written );
    _set_LC_NUMERIC_LOCALE_back();
    return( ImageIO_WRITING_HEADER );
  }


  /* close file descriptor */
  ImageIO_close( im );

  im->fd = NULL;
  im->openMode = OM_CLOSE;

  _set_LC_NUMERIC_LOCALE_back();

  return r;
}



/****************************************************************
 *
 * these procedures came from the nifti c lib
 *
 ****************************************************************/



typedef struct {                   /** 3x3 matrix struct **/
  float m[3][3] ;
} _mat33 ;



/*----------------------------------------------------------------------*/
/*! compute the determinant of a 3x3 matrix
*//*--------------------------------------------------------------------*/
static float _nifti_mat33_determ( _mat33 R )   /* determinant of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   return r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;
}



/*----------------------------------------------------------------------*/
/*! compute the max row norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
static float _nifti_mat33_rownorm( _mat33 A )  /* max row norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = fabs(A.m[0][0])+fabs(A.m[0][1])+fabs(A.m[0][2]) ;
   r2 = fabs(A.m[1][0])+fabs(A.m[1][1])+fabs(A.m[1][2]) ;
   r3 = fabs(A.m[2][0])+fabs(A.m[2][1])+fabs(A.m[2][2]) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}



/*----------------------------------------------------------------------*/
/*! compute the max column norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
static float _nifti_mat33_colnorm( _mat33 A )  /* max column norm of 3x3 matrix */
{
   float r1,r2,r3 ;

   r1 = fabs(A.m[0][0])+fabs(A.m[1][0])+fabs(A.m[2][0]) ;
   r2 = fabs(A.m[0][1])+fabs(A.m[1][1])+fabs(A.m[2][1]) ;
   r3 = fabs(A.m[0][2])+fabs(A.m[1][2])+fabs(A.m[2][2]) ;
   if( r1 < r2 ) r1 = r2 ;
   if( r1 < r3 ) r1 = r3 ;
   return r1 ;
}



/*----------------------------------------------------------------------*/
/*! compute the inverse of a 3x3 matrix
*//*--------------------------------------------------------------------*/
static _mat33 _nifti_mat33_inverse( _mat33 R )   /* inverse of 3x3 matrix */
{
   double r11,r12,r13,r21,r22,r23,r31,r32,r33 , deti ;
   _mat33 Q ;
                                                       /*  INPUT MATRIX:  */
   r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
   r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
   r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

   deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
         +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

   if( deti != 0.0l ) deti = 1.0l / deti ;

   Q.m[0][0] = deti*( r22*r33-r32*r23) ;
   Q.m[0][1] = deti*(-r12*r33+r32*r13) ;
   Q.m[0][2] = deti*( r12*r23-r22*r13) ;

   Q.m[1][0] = deti*(-r21*r33+r31*r23) ;
   Q.m[1][1] = deti*( r11*r33-r31*r13) ;
   Q.m[1][2] = deti*(-r11*r23+r21*r13) ;

   Q.m[2][0] = deti*( r21*r32-r31*r22) ;
   Q.m[2][1] = deti*(-r11*r32+r31*r12) ;
   Q.m[2][2] = deti*( r11*r22-r21*r12) ;

   return Q ;
}



/*---------------------------------------------------------------------------*/
/*! polar decomposition of a 3x3 matrix

   This finds the closest orthogonal matrix to input A
   (in both the Frobenius and L2 norms).

   Algorithm is that from NJ Higham, SIAM J Sci Stat Comput, 7:1160-1174.
*//*-------------------------------------------------------------------------*/
static _mat33 _nifti_mat33_polar( _mat33 A )
{
   _mat33 X , Y , Z ;
   float alp,bet,gam,gmi , dif=1.0 ;
   int k=0 ;

   X = A ;

   /* force matrix to be nonsingular */

   gam = _nifti_mat33_determ(X) ;
   while( gam == 0.0 ){        /* perturb matrix */
     gam = 0.00001 * ( 0.001 + _nifti_mat33_rownorm(X) ) ;
     X.m[0][0] += gam ; X.m[1][1] += gam ; X.m[2][2] += gam ;
     gam = _nifti_mat33_determ(X) ;
   }

   while(1){
     Y = _nifti_mat33_inverse(X) ;
     if( dif > 0.3 ){     /* far from convergence */
       alp = sqrt( _nifti_mat33_rownorm(X) * _nifti_mat33_colnorm(X) ) ;
       bet = sqrt( _nifti_mat33_rownorm(Y) * _nifti_mat33_colnorm(Y) ) ;
       gam = sqrt( bet / alp ) ;
       gmi = 1.0 / gam ;
     } else {
       gam = gmi = 1.0 ;  /* close to convergence */
     }
     Z.m[0][0] = 0.5 * ( gam*X.m[0][0] + gmi*Y.m[0][0] ) ;
     Z.m[0][1] = 0.5 * ( gam*X.m[0][1] + gmi*Y.m[1][0] ) ;
     Z.m[0][2] = 0.5 * ( gam*X.m[0][2] + gmi*Y.m[2][0] ) ;
     Z.m[1][0] = 0.5 * ( gam*X.m[1][0] + gmi*Y.m[0][1] ) ;
     Z.m[1][1] = 0.5 * ( gam*X.m[1][1] + gmi*Y.m[1][1] ) ;
     Z.m[1][2] = 0.5 * ( gam*X.m[1][2] + gmi*Y.m[2][1] ) ;
     Z.m[2][0] = 0.5 * ( gam*X.m[2][0] + gmi*Y.m[0][2] ) ;
     Z.m[2][1] = 0.5 * ( gam*X.m[2][1] + gmi*Y.m[1][2] ) ;
     Z.m[2][2] = 0.5 * ( gam*X.m[2][2] + gmi*Y.m[2][2] ) ;

     dif = fabs(Z.m[0][0]-X.m[0][0])+fabs(Z.m[0][1]-X.m[0][1])
          +fabs(Z.m[0][2]-X.m[0][2])+fabs(Z.m[1][0]-X.m[1][0])
          +fabs(Z.m[1][1]-X.m[1][1])+fabs(Z.m[1][2]-X.m[1][2])
          +fabs(Z.m[2][0]-X.m[2][0])+fabs(Z.m[2][1]-X.m[2][1])
          +fabs(Z.m[2][2]-X.m[2][2])                          ;

     k = k+1 ;
     if( k > 100 || dif < 3.e-6 ) break ;  /* convergence or exhaustion */
     X = Z ;
   }

   return Z ;
}





/*---------------------------------------------------------------------------*/
/*! Given the 3x4 upper corner of the matrix R, compute the quaternion
   parameters that fit it.

     - Any NULL pointer on input won't get assigned (e.g., if you don't want
       dx,dy,dz, just pass NULL in for those pointers).
     - If the 3 input matrix columns are NOT orthogonal, they will be
       orthogonalized prior to calculating the parameters, using
       the polar decomposition to find the orthogonal matrix closest
       to the column-normalized input matrix.
     - However, if the 3 input matrix columns are NOT orthogonal, then
       the matrix produced by nifti_quatern_to_mat44 WILL have orthogonal
       columns, so it won't be the same as the matrix input here.
       This "feature" is because the NIFTI 'qform' transform is
       deliberately not fully general -- it is intended to model a volume
       with perpendicular axes.
     - If the 3 input matrix columns are not even linearly independent,
       you'll just have to take your luck, won't you?

   \see "QUATERNION REPRESENTATION OF ROTATION MATRIX" in nifti1.h

   \see nifti_quatern_to_mat44, nifti_make_orthog_mat44,
       nifti_mat44_to_orientation
*//*-------------------------------------------------------------------------*/
static void _nifti_mat44_to_quatern( _image *im )
{
   double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
   double xd,yd,zd, a,b,c,d ;
   _mat33 P,Q ;


   /* offset outputs are read write out of input matrix  */

   im->tx = im->qform_toreal[0][3] ; im->ty = im->qform_toreal[1][3] ; im->tz = im->qform_toreal[2][3] ;

   /* load 3x3 matrix into local variables */

   r11 =  im->qform_toreal[0][0] ; r12 =  im->qform_toreal[0][1] ; r13 =  im->qform_toreal[0][2] ;
   r21 =  im->qform_toreal[1][0] ; r22 =  im->qform_toreal[1][1] ; r23 =  im->qform_toreal[1][2] ;
   r31 =  im->qform_toreal[2][0] ; r32 =  im->qform_toreal[2][1] ; r33 =  im->qform_toreal[2][2] ;

   /* compute lengths of each column; these determine grid spacings  */

   xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
   yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
   zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;

   /* if a column length is zero, patch the trouble */

   if( xd == 0.0l ){ r11 = 1.0l ; r21 = r31 = 0.0l ; xd = 1.0l ; }
   if( yd == 0.0l ){ r22 = 1.0l ; r12 = r32 = 0.0l ; yd = 1.0l ; }
   if( zd == 0.0l ){ r33 = 1.0l ; r13 = r23 = 0.0l ; zd = 1.0l ; }

   /* xd, yd, and zd are new estimates for voxel size
    * we assume we already have the good ones
    */

   /* normalize the columns */

   r11 /= xd ; r21 /= xd ; r31 /= xd ;
   r12 /= yd ; r22 /= yd ; r32 /= yd ;
   r13 /= zd ; r23 /= zd ; r33 /= zd ;

   /* At this point, the matrix has normal columns, but we have to allow
      for the fact that the hideous user may not have given us a matrix
      with orthogonal columns.

      So, now find the orthogonal matrix closest to the current matrix.

      One reason for using the polar decomposition to get this
      orthogonal matrix, rather than just directly orthogonalizing
      the columns, is so that inputting the inverse matrix to R
      will result in the inverse orthogonal matrix at this point.
      If we just orthogonalized the columns, this wouldn't necessarily hold. */

   Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
   Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
   Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;

   P = _nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

   r11 = P.m[0][0] ; r12 = P.m[0][1] ; r13 = P.m[0][2] ; /* unload */
   r21 = P.m[1][0] ; r22 = P.m[1][1] ; r23 = P.m[1][2] ;
   r31 = P.m[2][0] ; r32 = P.m[2][1] ; r33 = P.m[2][2] ;

   /*                            [ r11 r12 r13 ]               */
   /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
   /*                            [ r31 r32 r33 ]               */

   /* compute the determinant to determine if it is proper */

   zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
       +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */

   if( zd > 0 ){             /* proper */
     im->qfac = 1.0;
   } else {                  /* improper ==> flip 3rd column */
     im->qfac = -1.0;
     r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
   }

   /* now, compute quaternion parameters */

   a = r11 + r22 + r33 + 1.0l ;

   if( a > 0.5l ){                /* simplest case */
     a = 0.5l * sqrt(a) ;
     b = 0.25l * (r32-r23) / a ;
     c = 0.25l * (r13-r31) / a ;
     d = 0.25l * (r21-r12) / a ;
   } else {                       /* trickier case */
     xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
     yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
     zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
     if( xd > 1.0 ){
       b = 0.5l * sqrt(xd) ;
       c = 0.25l* (r12+r21) / b ;
       d = 0.25l* (r13+r31) / b ;
       a = 0.25l* (r32-r23) / b ;
     } else if( yd > 1.0 ){
       c = 0.5l * sqrt(yd) ;
       b = 0.25l* (r12+r21) / c ;
       d = 0.25l* (r23+r32) / c ;
       a = 0.25l* (r13-r31) / c ;
     } else {
       d = 0.5l * sqrt(zd) ;
       b = 0.25l* (r13+r31) / d ;
       c = 0.25l* (r23+r32) / d ;
       a = 0.25l* (r21-r12) / d ;
     }
     if( a < 0.0l ){ b=-b ; c=-c ; d=-d; a=-a; }
   }

   im->qb = b;
   im->qc = c;
   im->qd = d;
}




static void _nifti_quatern_to_mat44( _image *im )
{
   double a,b=im->qb,c=im->qc,d=im->qd;

   /* last row is always [ 0 0 0 1 ] */

   im->qform_toreal[3][0] = im->qform_toreal[3][1] = im->qform_toreal[3][2] = 0.0 ;
   im->qform_toreal[3][3]= 1.0 ;

   /* compute a parameter from b,c,d */

   a = 1.0l - (b*b + c*c + d*d) ;
   if( a < 1.e-7l ){                   /* special case */
     a = 1.0l / sqrt(b*b+c*c+d*d) ;
     b *= a ; c *= a ; d *= a ;        /* normalize (b,c,d) vector */
     a = 0.0l ;                        /* a = 0 ==> 180 degree rotation */
   } else{
     a = sqrt(a) ;                     /* angle = 2*arccos(a) */
   }

   /* load rotation matrix, including scaling factors for voxel sizes */

   im->qform_toreal[0][0] =        (a*a+b*b-c*c-d*d) * im->vx ;
   im->qform_toreal[0][1] = 2.0l * (b*c-a*d        ) * im->vy ;
   im->qform_toreal[0][2] = 2.0l * (b*d+a*c        ) * im->vz ;
   im->qform_toreal[1][0] = 2.0l * (b*c+a*d        ) * im->vx ;
   im->qform_toreal[1][1] =        (a*a+c*c-b*b-d*d) * im->vy ;
   im->qform_toreal[1][2] = 2.0l * (c*d-a*b        ) * im->vz ;
   im->qform_toreal[2][0] = 2.0l * (b*d-a*c        ) * im->vx ;
   im->qform_toreal[2][1] = 2.0l * (c*d+a*b        ) * im->vy ;
   im->qform_toreal[2][2] =        (a*a+d*d-c*c-b*b) * im->vz ;
   if ( im->qfac < 0 ) {
     im->qform_toreal[0][2] = -im->qform_toreal[0][2];
     im->qform_toreal[1][2] = -im->qform_toreal[1][2];
     im->qform_toreal[2][2] = -im->qform_toreal[2][2];
   }

   /* load offsets */

   im->qform_toreal[0][3] = im->tx ; im->qform_toreal[1][3] = im->ty ; im->qform_toreal[2][3] = im->tz ;

}
