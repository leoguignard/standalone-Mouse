#include <string.h>

#include <tiffio.h>

#include <tif.h>

static int _verbose_ = 1;
static int _debug_ = 0;

/** Magic header for TIFF files written in little endian format "II"
 */
#define TIFF_LE_MAGIC "\111\111\052\000"
/** Magic header for TIFF files written in big endian format "MM"
 */
#define TIFF_BE_MAGIC "\115\115\000\052"





/************************************************************
 *
 * inspired from
 * InsightToolkit-4.9.1/Modules/IO/TIFF/src/itkTIFFReaderInternal.[h,cxx]
 *
 ************************************************************/

typedef enum _imageFormat {
  FORMAT_RGB,
  FORMAT_GRAYSCALE,
  FORMAT_PALETTE_RGB,
  FORMAT_PALETTE_GRAYSCALE,
  FORMAT_OTHER
} _imageFormat;

typedef struct _internalTiffImage {
  TIFF *         m_Image;
  int           m_IsOpen;
  unsigned int       m_Width;
  unsigned int       m_Height;
  unsigned short int        m_NumberOfPages;
  unsigned short int        m_CurrentPage;
  unsigned short int        m_SamplesPerPixel;
  unsigned short int        m_Compression;
  unsigned short int        m_BitsPerSample;
  unsigned short int        m_Photometrics;
  int           m_HasValidPhotometricInterpretation;
  unsigned short int        m_PlanarConfig;
  unsigned short int        m_Orientation;
  unsigned int        m_TileRows;
  unsigned int        m_TileColumns;
  unsigned int        m_TileWidth;
  unsigned int        m_TileHeight;
  unsigned int        m_NumberOfTiles;
  unsigned int        m_SubFiles;
  unsigned int        m_IgnoredSubFiles;
  unsigned short int        m_ResolutionUnit;
  float          m_XResolution;
  float          m_YResolution;
  unsigned short int        m_SampleFormat;

  void *m_Buffer;
  _imageFormat m_ImageFormat;
  unsigned short int *m_ColorRed;
  unsigned short int *m_ColorGreen;
  unsigned short int *m_ColorBlue;
  int m_TotalColors;
} _internalTiffImage;



static void _InitInternalTiffImage( _internalTiffImage *im )
{
    im->m_Image = (TIFF*)NULL;
    im->m_Width = 0;
    im->m_Height = 0;
    im->m_SamplesPerPixel = 0;
    im->m_Compression = 0;
    im->m_BitsPerSample = 0;
    im->m_Photometrics = 0;
    im->m_HasValidPhotometricInterpretation = 0;
    im->m_PlanarConfig = 0;
    im->m_CurrentPage = 0;
    im->m_NumberOfPages = 0;
    im->m_NumberOfTiles = 0;
    im->m_Orientation = ORIENTATION_TOPLEFT;
    im->m_TileRows = 0;
    im->m_TileColumns = 0;
    im->m_TileWidth = 0;
    im->m_TileHeight = 0;
    im->m_XResolution = 1;
    im->m_YResolution = 1;
    im->m_SubFiles = 0;
    im->m_IgnoredSubFiles = 0;
    im->m_SampleFormat = 1;
    im->m_ResolutionUnit = 1; /* none */
    im->m_IsOpen = 0;

    im->m_ImageFormat = FORMAT_OTHER;
    im->m_Buffer = (void*)NULL;
    im->m_ColorRed = (unsigned short int*)NULL;
    im->m_ColorGreen = (unsigned short int*)NULL;
    im->m_ColorBlue = (unsigned short int*)NULL;
    im->m_TotalColors = 0;
}



static void _CleanInternalTiffImage( _internalTiffImage *im )
{
    if ( im->m_Image ) {
      TIFFClose(im->m_Image);
    }
    _InitInternalTiffImage( im );
}



/* inspired from TIFFImageIO::GetFormat()
 * from Modules/IO/TIFF/src/itkTIFFReaderInternal.cxx
 */
static int _InitializeInternalTiffImage( _internalTiffImage *im )
{
  char *proc = "_InitializeInternalTiffImage";
  unsigned int page;
  int subfiletype, cc;

  if ( im->m_Image ) {
    if ( !TIFFGetField(im->m_Image, TIFFTAG_IMAGEWIDTH, &im->m_Width)
         || !TIFFGetField(im->m_Image, TIFFTAG_IMAGELENGTH, &im->m_Height) ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to get image length or width\n", proc );
      return( -1 );
    }

    /* Get the resolution in each direction
     */
    TIFFGetField( im->m_Image,
                  TIFFTAG_XRESOLUTION, &im->m_XResolution);
    TIFFGetField( im->m_Image,
                  TIFFTAG_YRESOLUTION, &im->m_YResolution);
    TIFFGetField( im->m_Image,
                  TIFFTAG_RESOLUTIONUNIT, &im->m_ResolutionUnit);

    /* Check the number of pages. First by looking at the number of directories
     */
    im->m_NumberOfPages = TIFFNumberOfDirectories(im->m_Image);

    if ( im->m_NumberOfPages == 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: no directories found in TIFF file\n", proc );
      return( -1 );
    }

    if ( TIFFIsTiled(im->m_Image) ) {
      im->m_NumberOfTiles = TIFFNumberOfTiles(im->m_Image);

      if ( !TIFFGetField(im->m_Image, TIFFTAG_TILEWIDTH, &im->m_TileWidth)
           || !TIFFGetField(im->m_Image, TIFFTAG_TILELENGTH, &im->m_TileHeight) ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get tile length or width\n", proc );
        return( -1 );
      }
      else {
        im->m_TileRows = im->m_Height / im->m_TileHeight;
        im->m_TileColumns = im->m_Width / im->m_TileWidth;
      }
    }

    /* Checking if the TIFF contains subfiles
     */
    if ( im->m_NumberOfPages > 1 ) {
      im->m_SubFiles = 0;
      im->m_IgnoredSubFiles = 0;

      for ( page = 0; page < im->m_NumberOfPages; page++ ) {
        subfiletype = 6;
        if ( TIFFGetField(im->m_Image, TIFFTAG_SUBFILETYPE, &subfiletype) ) {
          if ( subfiletype == 0 ) {
            im->m_SubFiles += 1;
          }
          /* ignored flags
           */
          else if ( subfiletype & FILETYPE_REDUCEDIMAGE
                      || subfiletype & FILETYPE_MASK ) {
            ++im->m_IgnoredSubFiles;
          }

        }
        TIFFReadDirectory(im->m_Image);
      }

      /* Set the directory to the first image, and reads it
       */
      TIFFSetDirectory(im->m_Image, 0);
    }

    TIFFGetFieldDefaulted( im->m_Image, TIFFTAG_ORIENTATION,
                           &im->m_Orientation);
    TIFFGetFieldDefaulted( im->m_Image, TIFFTAG_SAMPLESPERPIXEL,
                           &im->m_SamplesPerPixel);
    TIFFGetFieldDefaulted( im->m_Image, TIFFTAG_COMPRESSION, &im->m_Compression);
    TIFFGetFieldDefaulted( im->m_Image, TIFFTAG_BITSPERSAMPLE,
                           &im->m_BitsPerSample);
    TIFFGetFieldDefaulted(im->m_Image, TIFFTAG_PLANARCONFIG, &im->m_PlanarConfig);
    TIFFGetFieldDefaulted(im->m_Image, TIFFTAG_SAMPLEFORMAT, &im->m_SampleFormat);

    /* If TIFFGetField returns false, there's no Photometric Interpretation
     * set for this image, but that's a required field so we set a warning flag.
     * (Because the "Photometrics" field is an enum, we can't rely on setting
     * im->m_Photometrics to some signal value.)
     */
    if ( TIFFGetField(im->m_Image, TIFFTAG_PHOTOMETRIC, &im->m_Photometrics) ) {
        im->m_HasValidPhotometricInterpretation = 1;
    }
    else {
        im->m_HasValidPhotometricInterpretation = 0;
    }
  }


  /* inspired from TIFFImageIO::GetFormat()
   * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
   */
  switch ( im->m_Photometrics ) {
  default :
    im->m_ImageFormat = FORMAT_OTHER;
    break;
  case PHOTOMETRIC_RGB:
  case PHOTOMETRIC_YCBCR:
    im->m_ImageFormat = FORMAT_RGB;
    break;
  case PHOTOMETRIC_MINISWHITE:
  case PHOTOMETRIC_MINISBLACK:
    im->m_ImageFormat = FORMAT_GRAYSCALE;
    break;
  case PHOTOMETRIC_PALETTE:
    /* inspired from TIFFImageIO::InitializeColors()
     * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
     */
    if ( TIFFGetField(im->m_Image, TIFFTAG_COLORMAP,
                       &(im->m_ColorRed), &(im->m_ColorGreen), &(im->m_ColorBlue)) == 0 ) {
      im->m_ColorRed = im->m_ColorGreen = im->m_ColorBlue = NULL;
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to get color map\n", proc );
      return( -1 );
    }
    switch ( im->m_BitsPerSample ) {
    case 8:
        im->m_TotalColors = 256;
        break;
    case 16:
        im->m_TotalColors = 65536;
        break;
      break;
    default:
        im->m_ColorRed = im->m_ColorGreen = im->m_ColorBlue = NULL;
        if ( _verbose_ )
          fprintf( stderr, "%s: such image type not handled for color map\n", proc );
        return( -1 );
    }
    im->m_ImageFormat = FORMAT_PALETTE_GRAYSCALE;
    for ( cc=0; cc<im->m_TotalColors; cc++ ) {
      if ( im->m_ColorRed[cc] != im->m_ColorGreen[cc]
           || im->m_ColorRed[cc] != im->m_ColorBlue[cc])
        im->m_ImageFormat = FORMAT_PALETTE_RGB;
    }
    break;
  }

  return( 1 );
}



static int _CanReadInternalTiffImage( _internalTiffImage *im )
{
    int compressionSupported = ( TIFFIsCODECConfigured(im->m_Compression) == 1 );
    return ( im->m_Image && ( im->m_Width > 0 ) && ( im->m_Height > 0 )
             && ( im->m_SamplesPerPixel > 0 )
             && compressionSupported
             && ( im->m_NumberOfTiles == 0 )
             /* just use TIFFReadRGBAImage, an
              * native optimized version would be nice
              */
             && ( im->m_HasValidPhotometricInterpretation )
             && ( im->m_Photometrics == PHOTOMETRIC_RGB
                  || im->m_Photometrics == PHOTOMETRIC_MINISWHITE
                  || im->m_Photometrics == PHOTOMETRIC_MINISBLACK
                  || ( im->m_Photometrics == PHOTOMETRIC_PALETTE
                       && im->m_BitsPerSample != 32 )
                  )
             && ( im->m_PlanarConfig == PLANARCONFIG_CONTIG )
             && ( im->m_Orientation == ORIENTATION_TOPLEFT
                  || im->m_Orientation == ORIENTATION_BOTLEFT )
             && ( im->m_BitsPerSample == 8 || im->m_BitsPerSample == 16 || im->m_BitsPerSample == 32 ) );

}



static int _OpenInternalTiffImage( const char *filename, _internalTiffImage *im )
{
  char *proc = "_OpenInternalTiffImage";

  _CleanInternalTiffImage( im );

  im->m_Image = TIFFOpen(filename, "r");
  if ( !im->m_Image ) {
    _CleanInternalTiffImage( im );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading '%s' from libtiff\n", proc, filename );
    return( 0 );
  }
  if ( _InitializeInternalTiffImage( im ) != 1 ) {
    _CleanInternalTiffImage( im );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when extraction information from '%s'\n", proc, filename );
    return( 0 );
  }

  im->m_IsOpen = 1;
  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/




/************************************************************
 *
 *
 *
 ************************************************************/


int testTifHeader( char *magic __attribute__ ((unused)),
                   const char *name )
{
    /*
    if( !memcmp(magic,TIFF_LE_MAGIC,4) ||
        !memcmp(magic,TIFF_BE_MAGIC,4))
      return 0;
    */
    if ( (!strncmp(name+strlen(name)-4, ".tif", 4))
         || (!strncmp(name+strlen(name)-4, ".TIF", 4))
         || (!strncmp(name+strlen(name)-5, ".tiff", 5))
         || (!strncmp(name+strlen(name)-5, ".TIFF", 5)) )
      return 0;
    else
      return -1;
}


/************************************************************
 *
 *
 *
 ************************************************************/



/* inspired from TIFFImageIO::ReadGenericImage()
 * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
 */
static int _ReadGenericImage( _internalTiffImage *m_InternalImage, _image *im )
{
  char *proc = "_ReadGenericImage";
  int isize, inc;
  unsigned int row, x;
  void *buf, *image;


  isize = TIFFScanlineSize(m_InternalImage->m_Image);
  buf = (void*)malloc( isize );
  if ( buf == (void*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
    return( -1 );
  }

  if ( m_InternalImage->m_PlanarConfig != PLANARCONFIG_CONTIG ) {
    free( buf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable deal with something else than PLANARCONFIG_CONTIG\n", proc );
    return( -1 );
  }
  if ( m_InternalImage->m_Orientation != ORIENTATION_TOPLEFT
       && m_InternalImage->m_Orientation != ORIENTATION_BOTLEFT ) {
    free( buf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable deal with something else than ORIENTATION_TOPLEFT or ORIENTATION_BOTLEFT\n", proc );
    return( -1 );
  }


  switch ( m_InternalImage->m_ImageFormat ) {
  case FORMAT_GRAYSCALE :
  case FORMAT_PALETTE_GRAYSCALE :
      inc = 1;
      break;
  case FORMAT_RGB :
      inc = m_InternalImage->m_SamplesPerPixel;
      break;
  case FORMAT_PALETTE_RGB :
      inc = 3;
      break;
  default :
    free( buf );
    if ( _verbose_ )
      fprintf( stderr, "%s: such image format not handled yet\n", proc );
    return( -1 );
  }


  for ( row=0; row<im->ydim; row++ ) {

    if ( TIFFReadScanline(m_InternalImage->m_Image, buf, row, 0) <= 0 ) {
        free( buf );
        if ( _verbose_ )
          fprintf( stderr, "%s: problem reading row #%d\n", proc, row );
        return( -1 );
    }

    if ( m_InternalImage->m_Orientation == ORIENTATION_TOPLEFT ) {
      image = m_InternalImage->m_Buffer + (size_t) (row) * im->xdim * im->wdim * inc;
    }
    else {
      /* bottom left
       */
      image = m_InternalImage->m_Buffer + (size_t) (im->xdim) * im->wdim * inc * ( im->ydim - ( row + 1 ) );
    }

    switch ( m_InternalImage->m_ImageFormat ) {
    case FORMAT_GRAYSCALE :
      memcpy( image, buf, im->xdim*im->vdim*im->wdim );
      break;
    case FORMAT_RGB :
      memcpy( image, buf, im->xdim*im->vdim*im->wdim );
      break;

    case FORMAT_PALETTE_GRAYSCALE :
      if ( im->wordKind == WK_FIXED && im->sign == SGN_UNSIGNED) {
        switch ( im->wdim ) {
        default :
          free( buf );
          if ( _verbose_ )
            fprintf( stderr, "%s: such word dimension not handled\n", proc );
          return( -1 );
        case 1 :
          {
            unsigned char *fbuf = (unsigned char*)buf;
            unsigned char *tbuf = (unsigned char*)image;
            for ( x=0; x<im->xdim; x++ )
              tbuf[x] = (unsigned char)(m_InternalImage->m_ColorRed[fbuf[x]]);
          }
          break;
        case 2 :
          {
            unsigned short int *fbuf = (unsigned short int*)buf;
            unsigned short int *tbuf = (unsigned short int*)image;
            for ( x=0; x<im->xdim; x++ )
              tbuf[x] = (unsigned short int)(m_InternalImage->m_ColorRed[fbuf[x]]);
          }
          break;
        }
      }
      else {
        free( buf );
        if ( _verbose_ )
          fprintf( stderr, "%s: such image type not handled\n", proc );
        return( -1 );
      }
      break;

    case FORMAT_PALETTE_RGB :
        if ( im->wordKind == WK_FIXED && im->sign == SGN_UNSIGNED) {
          switch ( im->wdim ) {
          default :
            free( buf );
            if ( _verbose_ )
              fprintf( stderr, "%s: such word dimension not handled\n", proc );
            return( -1 );
          case 1 :
            {
              unsigned char *fbuf = (unsigned char*)buf;
              unsigned char *tbuf = (unsigned char*)image;
              for ( x=0; x<im->xdim; x++ ) {
                tbuf[3*x]   = (unsigned char)(m_InternalImage->m_ColorRed[fbuf[x]]);
                tbuf[3*x+1] = (unsigned char)(m_InternalImage->m_ColorGreen[fbuf[x]]);
                tbuf[3*x+2] = (unsigned char)(m_InternalImage->m_ColorBlue[fbuf[x]]);
              }
            }
            break;
          case 2 :
            {
              unsigned short int *fbuf = (unsigned short int*)buf;
              unsigned short int *tbuf = (unsigned short int*)image;
              for ( x=0; x<im->xdim; x++ ) {
                tbuf[3*x]   = (unsigned char)(m_InternalImage->m_ColorRed[fbuf[x]]);
                tbuf[3*x+1] = (unsigned char)(m_InternalImage->m_ColorGreen[fbuf[x]]);
                tbuf[3*x+2] = (unsigned char)(m_InternalImage->m_ColorBlue[fbuf[x]]);
              }
           }
            break;
          }
        }
        else {
          free( buf );
          if ( _verbose_ )
            fprintf( stderr, "%s: such image type not handled\n", proc );
          return( -1 );
        }
        break;

    default :
      free( buf );
      if ( _verbose_ )
        fprintf( stderr, "%s: such image format not handled yet\n", proc );
      return( -1 );
    }

  }


  free( buf );
  return( 1 );
}



/* inspired from TIFFImageIO::ReadCurrentPage()
 * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
 */
static int _readTiffPage( _internalTiffImage *tiffImage, _image *im )
{
  char *proc = "_readTiffPage";

  if ( _CanReadInternalTiffImage( tiffImage ) == 0 ) {
    if ( _verbose_ )
        fprintf( stderr, "%s: this part is not implemented yet\n", proc );
    return( -1 );
  }
  else {

    /* to handle palette based color image
     * one has to implement the equivalent of
     * TIFFImageIO::InitializeColors() here
     */
    if ( ( im->wordKind == WK_FIXED && (im->wdim == 1 || im->wdim == 2) )
         || ( im->wordKind == WK_FLOAT && im->wdim == 4 ) ) {
      if ( _ReadGenericImage( tiffImage, im ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: error when reading image\n", proc );
        return( -1 );
      }
    }
    else {
        if ( _verbose_ )
          fprintf( stderr, "%s: such image type not handled\n", proc );
        return( -1 );
    }

  }

  return( 1 );
}



/* inspired from TIFFImageIO::ReadVolume()
 * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
 */
static int _readTiffVolume( _internalTiffImage *tiffImage, _image *im )
{
  char *proc = "_readTiffVolume";

  int page;
  int32 subfiletype = 6;

  for ( page = 0; page < tiffImage->m_NumberOfPages; page++ ) {
    if ( _debug_ ) {
      fprintf( stderr, "%s: read page #%d/%d\n", proc, page, tiffImage->m_NumberOfPages );
    }
    if ( tiffImage->m_IgnoredSubFiles > 0 ) {
      if ( TIFFGetField(tiffImage->m_Image, TIFFTAG_SUBFILETYPE, &subfiletype) ) {
        if ( subfiletype & FILETYPE_REDUCEDIMAGE
             || subfiletype & FILETYPE_MASK ) {
          TIFFReadDirectory(tiffImage->m_Image);
          continue;
        }
      }
    }

    _readTiffPage( tiffImage, im );

    if ( page < tiffImage->m_NumberOfPages-1 )
        tiffImage->m_Buffer += im->xdim * im->ydim * im->vdim * im->wdim;

    TIFFReadDirectory( tiffImage->m_Image );
  }

  return( 1 );
}






/* inspired from TIFFImageIO::ReadImageInformation()
 * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
 */

static int _tiffImageToImageIOImage( _internalTiffImage *m_InternalImage, _image *im )
{
  char *proc = "_tiffImageToImageIOImage";
  size_t size;

  if ( m_InternalImage->m_IsOpen == 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: tiff image in not opened\n", proc );
    return( -1 );
  }

  /* dimension
   * how to get zdim ?
   */
  im->xdim = m_InternalImage->m_Width;
  im->ydim = m_InternalImage->m_Height;
  im->zdim = 1;
  im->vdim = m_InternalImage->m_SamplesPerPixel;

  switch ( m_InternalImage->m_ImageFormat ) {
  case FORMAT_GRAYSCALE :
  case FORMAT_PALETTE_GRAYSCALE :
      im->vdim = 1;
      break;
  case FORMAT_RGB :
  case FORMAT_PALETTE_RGB :
      im->vdim = 3;
      break;
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image format not handled yet\n", proc );
    return( -1 );
  }


  /* pixel size
   */
  if ( m_InternalImage->m_ResolutionUnit > 0
       && m_InternalImage->m_XResolution > 0
       && m_InternalImage->m_YResolution > 0 ) {
    if ( m_InternalImage->m_ResolutionUnit == 2 ) {
      /* inches
       */
      im->vx = 25.4 / m_InternalImage->m_XResolution;
      im->vy = 25.4 / m_InternalImage->m_YResolution;
    }
    else if ( m_InternalImage->m_ResolutionUnit == 3 ) {
      /* cm
       */
      im->vx = 10.0 / m_InternalImage->m_XResolution;
      im->vy = 10.0 / m_InternalImage->m_YResolution;
    }
  }

  /* type
   */
  if ( m_InternalImage->m_BitsPerSample == 8 ) {
    if ( m_InternalImage->m_SampleFormat == 2 ) {
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
    }
    else {
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
    }
  }
  else   if ( m_InternalImage->m_BitsPerSample == 16 ) {
    if ( m_InternalImage->m_SampleFormat == 2 ) {
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
    }
    else {
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
    }
  }
  else if ( m_InternalImage->m_BitsPerSample == 32 ) {
    if ( m_InternalImage->m_SampleFormat == 3 ) {
        im->wdim = 4;
        im->wordKind = WK_FLOAT;
        im->sign = SGN_SIGNED;
      }
  }
  else {
      if ( _verbose_ )
        fprintf( stderr, "%s: image type not recognized\n", proc );
      return( -1 );
  }

  /* allocation
   */
  im->zdim = m_InternalImage->m_NumberOfPages;
  size = im->xdim * im->ydim * im->zdim * im->vdim;
  size *= im->wdim;
  im->data = (void*)malloc( size );
  if ( im->data == (void*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate image buffer\n", proc );
    return( -1 );
  }
  m_InternalImage->m_Buffer = im->data;

  return( 1 );
}





/* inspired from TIFFImageIO::Read(void *buffer)
 * from Modules/IO/TIFF/src/itkTIFFImageIO.cxx
 */

int readTiffImage( const char *name, _image *im )
{
  char *proc = "readTiffImage";
  _internalTiffImage tiffImage;

  _InitInternalTiffImage( &tiffImage );

  /* open file and read information
   */
  if ( _OpenInternalTiffImage( name, &tiffImage ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: can not open '%s'\n", proc, name );
    return( -1 );
  }

  if ( _tiffImageToImageIOImage( &tiffImage, im ) != 1 ) {
    _freeImage( im );
    _CleanInternalTiffImage( &(tiffImage) );
    _InitInternalTiffImage( &tiffImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate nifti header\n", proc );
    return( -1 );
  }

  if ( tiffImage.m_NumberOfPages > 0 ) {

    if ( _readTiffVolume( &tiffImage, im ) != 1 ) {
      _freeImage( im );
      _CleanInternalTiffImage( &(tiffImage) );
      _InitInternalTiffImage( &tiffImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when reading multiple pages tiff\n", proc );
      return( -1 );
    }

  }
  else if ( _readTiffPage( &tiffImage, im ) != 1 ) {
    _freeImage( im );
    _CleanInternalTiffImage( &(tiffImage) );
    _InitInternalTiffImage( &tiffImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading single page tiff\n", proc );
    return( -1 );
  }

  _CleanInternalTiffImage( &(tiffImage) );
  _InitInternalTiffImage( &tiffImage );

  return( 1 );
}






/************************************************************
 *
 *
 *
 ************************************************************/


/************************************************************
 *
 *
 *
 ************************************************************/



PTRIMAGE_FORMAT createTifFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat = &testTifHeader;
  f->readImageHeader = &readTiffImage;
  f->readImageData = NULL;
  f->writeImageHeader = NULL;
  f->writeImage = NULL;

  strcpy(f->fileExtension,".tif,.TIF,.tiff,.TIFF");
  strcpy(f->realName,"Tiff");
  return f;
}

