#ifndef IMAGEIO_H
#define IMAGEIO_H

#include <stdlib.h>
#include <stdio.h>



#ifdef NOZLIB
#undef ZLIB
#else
#define ZLIB
#endif

#ifdef ZLIB
#include <zlib.h>
/* see http://www.gzip.org/zlib/
   for details and documentation
*/
#endif



#ifdef __cplusplus
extern "C" {
#endif


#ifndef LONGINT

#if (defined _ALPHA_ || (defined _SGI_ && (defined _64_ || defined _64_M4_ || defined _64_M3_)))
/* the 64-bits type on 64-bits platform (long int) */
#define LONGINT long  int
#else
#ifdef __GNUC__
/* the 64-bits type on 32-bits platform (long long int) */
#define LONGINT long long int
#else
/*#define LONGINT __int64 */
#define LONGINT long int
#endif
#endif

#endif




/** set verbose and debug indicators */
extern void _SetVerboseInImageIO( int v );
extern void _IncrementVerboseInImageIO(  );
extern void _SetDebugInImageIO( int v );
extern void _IncrementDebugInImageIO(  );





/*--------------------------------------------------
 *
 * image reader/writer generic format
 *
 --------------------------------------------------*/

struct point_image;
#define IMAGE_FORMAT_NAME_LENGTH 100



/** defines the type of function called to test if an image is of a given
    format. The first parameter is an array of char of size 5 (ends with
    character 0) that describes the first file character (magic string). The
    second parameter is the filename. The output value is >=0 if the image is
    of that given format and <0 otherwise */
typedef int (*TEST_IMAGE_FORMAT)(char *,const char *);

/** defines the type of function called to read an image or an image header
    from a file corresponding to a given format. The first parameter is the
    file name whereas the second parameter is an _image structure. Note that
    the file has been already opened (the file descriptor fd is valid).
    The output value is >0  if  the whole image has been read, it is 0 if
    only the header has been read and it is  <0 otherwise */
typedef int (*READ_IMAGE_HEADER)(const char *, struct point_image *);
typedef int (*READ_IMAGE_DATA)(struct point_image *);

/** defines the type of function called to write an image to a file
    corresponding to a given format.
    The first parameter is the full file name whereas the second parameter
    is an _image structure.
    Note that the file has to be opened and closed in the function.
    The output value is >=0 if the whole image has been written
    correctly and it is <0 otherwise */
typedef int (*WRITE_IMAGE_HEADER)(char *,struct point_image *);
typedef int (*WRITE_IMAGE)(char *,struct point_image *);



/** Image Format descriptor */
typedef struct imformat {

  /** a pointer on a function that tests if an image is of a given format */
  TEST_IMAGE_FORMAT testImageFormat;

  /** a pointer on a function that reads the header of an image file */
  READ_IMAGE_HEADER readImageHeader;
  READ_IMAGE_DATA readImageData;

  /** a pointer on a function that writes  image of a given
      format */
  WRITE_IMAGE_HEADER writeImageHeader;
  WRITE_IMAGE writeImage;

  /* the file extension of format (including a dot ".": if several
     extensions may be used, they should be separed with a
     comma ".inr,.inr.gz" */
  char fileExtension[IMAGE_FORMAT_NAME_LENGTH];

  /** the usual name given to a format : for instance "inrimage", "gif" */
  char realName[IMAGE_FORMAT_NAME_LENGTH];
  /* pointer towards the next image format*/
  struct imformat *next;
} IMAGE_FORMAT, *PTRIMAGE_FORMAT;



extern void _initImageFormat( PTRIMAGE_FORMAT format );

extern void printSupportedFileFormat();





/*--------------------------------------------------
 *
 * endianness stuff
 *
 --------------------------------------------------*/



/** endianness */
typedef enum ENDIANNESS {
  /** Little endian processor */
  END_LITTLE,
  /** Big endian processor */
  END_BIG,
  /** Unknown endianness (unopenned file) */
  END_UNKNOWN
} ENDIANNESS;

/** returns the endianness of the hardware architecture */
ENDIANNESS  _getEndianness();



/*--------------------------------------------------
 *
 * image definition
 *
 --------------------------------------------------*/



/** file open mode */
typedef enum OPEN_MODE {
  /** no file open */
  OM_CLOSE,
  /** file is stdin or stdout */
  OM_STD,
  /** file is gzipped */
#ifdef ZLIB
  OM_GZ,
#endif
  /** normal file */
  OM_FILE
} OPEN_MODE;


/** data mode */
typedef enum DATA_MODE {
  /** data are binary */
  DM_BINARY,
  /** data are ascii */
  DM_ASCII
} DATA_MODE;


/** kind of image word */
typedef enum WORD_KIND {
  /** fixed type */
  WK_FIXED,
  /** floating point */
  WK_FLOAT,
  /** unknown (uninitialized) */
  WK_UNKNOWN
} WORD_KIND;


/** image word sign */
typedef enum SIGN {
  /** signed */
  SGN_SIGNED,
  /** unsigned */
  SGN_UNSIGNED,
  /** unknown (uninitialized or floating point words) */
  SGN_UNKNOWN
} SIGN;




/** inrimage vectorial storage mode */
typedef enum VECTORIAL_MODE {
  /** interlaced vectors (i.e. x1, y1, z1, x2, y2, z2, x3, y3, z3, ...) */
  VM_INTERLACED,
  /** non interlaced vectors (i.e. x1, x2, x3, ..., y1, y2, y3, ..., z1, z2, z3...) */
  VM_NON_INTERLACED,
  /** scalar inrimage */
  VM_SCALAR
} VECTORIAL_MODE;



#ifdef ZLIB
typedef gzFile _ImageIO_file;
#else
typedef FILE*  _ImageIO_file;
#endif



/** Image descriptor */
typedef struct point_image {
  /** Image x dimension (number of columns) */
  size_t xdim;
  /** Image y dimension (number of rows) */
  size_t ydim;
  /** Image z dimension (number of planes) */
  size_t zdim;
  /** Image vectorial dimension */
  unsigned int vdim;

  /** Image voxel size in x dimension */
  double vx;
  /** Image voxel size in y dimension */
  double vy;
  /** Image voxel size in z dimension */
  double vz;

  /* geometry stuff
   */
  int t_is_set;
  /** Image offset in x dimension */
  float tx;
  /** Image offset in y dimension */
  float ty;
  /** Image offset in z dimension */
  float tz;

  /** nifti-like quaternion definition
   */

  int q_is_set;
  float qb;
  float qc;
  float qd;

  int qfac_is_set;
  float qfac;

  /* matrix to transform pixel coordinates into
   * real coordinates
   * (default is multiplication by the voxel sizes)
   */
  int qform_code;
  int sform_code;
  float qform_toreal[4][4];
  float sform_toreal[4][4];


  /** Image data buffer */
  void *data;

  /** Image word size (in bytes) */
  unsigned int wdim;
  /** Image format to use for I/0. Should not be set by user */
  PTRIMAGE_FORMAT imageFormat;
  /** Data buffer vectors are interlaced or non interlaced */
  VECTORIAL_MODE vectMode;
  /** Image words kind */
  WORD_KIND wordKind;
  /** Image words sign */
  SIGN sign;

  /** User defined strings array. The user can use any internal purpose string.
      Each string is written at then end of header after a '#' character. */
  char **user;
  /** Number of user defined strings */
  unsigned int nuser;

  /** Image file descriptor */
  _ImageIO_file fd;


  /** Kind of image file descriptor */
  OPEN_MODE openMode;
  /** Written words endianness */
  ENDIANNESS endianness;
  /** Kind of image data encoding */
  DATA_MODE dataMode;

} _image;



/** Allocates and initializes an image descriptor */
_image *_initImage();


/* desallocates the _image structure
 * except for the data
 */
extern void _freeImageStructure( _image *im );

/** Free an image descriptor
    @param im image descriptor */
void _freeImage(_image *im);

void _swapImageData( _image *im );


/*--------------------------------------------------
 *
 * mimics standard routines
 *
 --------------------------------------------------*/

/** function prototype to allocate memory */
typedef void *(*ALLOCATION_FUNCTION)(size_t);

/** function prototype to free memory */
typedef void (*DEALLOCATION_FUNCTION)(void *);


/** set allocation and deallocation routines
    @param alloc new allocation routine
    @param del new deallocation routine */
void setImageIOAllocationRoutines(ALLOCATION_FUNCTION alloc,
          DEALLOCATION_FUNCTION del);



/** call allocation routine */
void *ImageIO_alloc(size_t);
/** call deallocation routine */
void ImageIO_free(void *);

/** replaces fwrite function
    @param im image to write
    @param buf data buffer to write
    @param len buffer length */
size_t ImageIO_write(const _image *im, const void *buf, size_t len);


/** replaces fread function
    @param im image to read
    @param buf data buffer to read
    @param len buffer length */
size_t ImageIO_read(const _image *im, void *buf, size_t len);

/** replaces fgets function
 */
char *ImageIO_gets( const _image *im, char *str, int size );

/** replaces fseek function
 */
int ImageIO_seek( const _image *im, long offset, int whence );

/** replaces ferror function
 */
int ImageIO_error( const _image *im );

int ImageIO_close( _image* im );


/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/

/** given an initialized file descriptor and a file name, open file
   from stdout (if name == NULL), a gziped pipe (if file is gziped)
   or a standard file otherwise.
   @param im initialized image descriptor
   @param name image file name */
void _openWriteImage(_image* im, const char *name) ;

/** open an image file from stdin (if name == NULL), from a pipe
   (piped with gzip if image was compressed) or from a standard file
   @param im initialized image descriptor
   @param name image file name */
void _openReadImage(_image *im, const char *name);





/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/



/** Error codes */
#define ImageIO_NO_ERROR 0
#define ImageIO_UNKNOWN_TYPE -1
#define ImageIO_OPENING -2
#define ImageIO_READING_HEADER -3
#define ImageIO_READING_IMAGE -4
#define ImageIO_WRITING_HEADER -3
#define ImageIO_WRITING_IMAGE -4
#define ImageIO_WRITING_DATA  -5


/*--------------------------------------------------
 *
 * image reading procedure
 *
 *--------------------------------------------------*/

/** Reads header from an image file<br>
    If file is an inrimage, only header is read. Otherwise, whole image<br>
    is read and image file descriptor is closed.<br>
    If name is NULL, header is read from STDIN
    @param name image file name or NULL */
_image* _readImageHeader(const char *name);

/** Reads body from an inrmage whose header has been read by
    _readImageHeader
    @param im image to read */
int _readImageData(_image *im);

/** Reads an image from a file and returns an image descriptor or NULL if<br>
    reading failed.<br>
    Reads from stdin if image name is NULL.
    The image data field points to a xdim * ydim * zdim * vdim buffer
    containing voxels in order:
    (Z1, Y1, X1, V1) (Z1, Y1, X1, V2), ... , (Z1, Y1, X1, Vt),
    (Z1, Y1, X2, V1) ...         ...       , (Z1, Y1, X2, Vt),
    ...
    (Z1, Y1, Xn, V1) ...         ...       , (Z1, Y1, Xn, Vt),
    (Z1, Y2, X1, V1) ...         ...       , (Z1, Y2, X1, Vt),
    ...
    (Z2, Y1, X1, V1) ...         ...       , (Z2, Y1, X1, Vt),
    ...
                     ...         ...       , (Zl, Ym, Xn, Vt)

    Read the following format:
    Inrimage,
    GIF,
    IRIS,
    ANALYSE,
    PGM,
    PPM,
    BMP,
    GIS (CEA, IRISA, ENST 3D image format).
    
    See also:
    http://www.dcs.ed.ac.uk/home/mxr/gfx/2d-hi.html and
    http://www.gzip.org/zlib/
    

   @param name image file name or NULL for stdin */
_image* _readImage(const char *name);



/*--------------------------------------------------
 *
 * image writing procedures
 *
 *--------------------------------------------------*/


int _writeImageHeader( _image *im, char *name_to_be_written );


/** Writes given image in file 'name'.<br>
    If name ends with '.gz', file is gzipped.<br>
    If name is NULL, image is sent to stdout.
    @param im image descriptor 
    @param name file name to store image or NULL */
int _writeImage(_image *im, char *name);


#ifdef __cplusplus
}
#endif

#endif
