
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <metaImage.h>


static int _verbose_ = 1;
static int _debug_ = 0;



static off_t getFileSize( char *name )
{
  char *proc = "getFileSize";
  struct stat stbuf;

  if ( name == (char*)NULL ) return( -1 );
  if ( stat( name, &stbuf ) != 0 ) {
    fprintf( stderr, "%s: unable to fill '%s' file info \n", proc, name );
    return( -1 );
  }
  return( stbuf.st_size );
}


/* get a string from a file and discard the ending newline character
   if any */
static char *fgetns( char *str, int n, _image *im ) {
  char *ret = NULL;
  int l;

  ret = ImageIO_gets( im, str, n );

  if(!ret) return NULL;

  l = strlen(str);
  if (l > 0 && str[l-1] == '\n') str[l-1] = '\0';
  return ret;
}



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}



static int _testTag( char *line, char *tag, char **tagarg )
{
  /* test the tag, return the position of the arg
   */
  if ( strncmp( line, tag, strlen(tag) ) == 0 ) {
    if ( tagarg != (char**)NULL ) {
      *tagarg = line;
      *tagarg += strlen(tag);
      while ( **tagarg == ' ' ) (*tagarg)++;
      if ( **tagarg == '=' ) (*tagarg)++;
      while ( **tagarg == ' ' ) (*tagarg)++;
    }
    return( 1 );
  }
  return( -1 );
}



static int _testPrintfFormat( char *line )
{
  char *s = line;
  while( *s != '\0' && *s != '\n' ) {
    if ( *s == '%' ) return( 1 );
    s++;
  }
  return( 0 );
}


/************************************************************
 *
 *
 *
 ************************************************************/



int testMetaImageHeader( char *magic __attribute__ ((unused)),
                         const char *name )
{
    if (( !strncmp(name+strlen(name)-4, ".mha", 4) ) ||
        ( !strncmp(name+strlen(name)-4, ".mha.gz", 7) ) ||
        ( !strncmp(name+strlen(name)-7, ".mhd", 4) )  )
      return 0;
   else
     return -1;
}





/************************************************************
 *
 *
 *
 ************************************************************/


int readMetaImage( const char *name, _image *im )
{
  char *proc = "readMetaImage";
  char string[1024];
  char rawname[1024];
  char *retfgetns, *tagarg, *s;

  int ndims = -1;
  int elementspacing = 0;
  int headersize = -2;
  unsigned long nread, size;

  fgetns( string, 1023, im );

  do {
    retfgetns = fgetns( string, 1023, im );
    if ( retfgetns == (char*)NULL ) continue;

    if ( 0 )
      fprintf( stderr, "- has read '%s'\n", string );

    /******************************
     * MetaObject Tags
     ******************************/
    if ( _testTag( string, "Comment", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ObjectType", &tagarg ) == 1 ) {
      if ( _testTag( tagarg, "Image", (char**)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: ObjectType is not an image ('%s')\n", proc, tagarg );
      }
    }
    else if ( _testTag( string, "ObjectSubType", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "TransformType", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "NDims", &tagarg ) == 1 ) {
      if ( sscanf( tagarg, "%d", &ndims ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get 'NDims'\n", proc );
        return( -1 );
      }
    }
    else if ( _testTag( string, "Name", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ID", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ParentID", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "BinaryData", &tagarg ) == 1 && string[10] == ' ' ) {
      if ( _testTag( tagarg, "True", (char**)NULL ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: BinaryData is not true ('%s')\n", proc, tagarg );
        return( -1 );
      }
    }
    else if ( _testTag( string, "ElementByteOrderMSB", &tagarg ) == 1
              || _testTag( string, "BinaryDataByteOrderMSB", &tagarg ) == 1 ) {
      if ( _testTag( tagarg, "True", (char**)NULL ) == 1 ) {
        im->endianness = END_BIG;
      }
      else if ( _testTag( tagarg, "False", (char**)NULL ) == 1 ) {
        im->endianness = END_LITTLE;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown ByteOrderMSB ('%s')\n", proc, tagarg );
        return( -1 );
      }
    }
    else if ( _testTag( string, "Color", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "Position", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "Orientation", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "AnatomicalOrientation", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ElementSpacing", &tagarg ) == 1 ) {
      if ( ndims == 1 ) {
        if ( sscanf( tagarg, "%lf", &im->vx ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
          return( -1 );
        }
        im->vy = im->vz = 1.0;
      }
      else if ( ndims == 2 ) {
        if ( sscanf( tagarg, "%lf %lf", &im->vx, &im->vy ) != 2 ) {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
            fprintf( stderr, "\t '%s' -> '%s'\n", string, tagarg );
          }
          return( -1 );
        }
        im->vz = 1.0;
      }
      else if ( ndims == 3 ) {
        if ( sscanf( tagarg, "%lf %lf %lf", &im->vx, &im->vy, &im->vz ) != 3 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
          return( -1 );
        }
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read voxel size (undefined NDims = %d)\n", proc, ndims );
        return( -1 );
      }
      elementspacing = 1;
    }
    /******************************
     * MetaImage Tags
     ******************************/
    else if ( _testTag( string, "DimSize", &tagarg ) == 1 ) {
        if ( ndims == 1 ) {
          if ( sscanf( tagarg, "%lu", &im->xdim ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to read image dim (NDims = %d)\n", proc, ndims );
            return( -1 );
          }
          im->ydim = im->zdim = 1;
        }
        else if ( ndims == 2 ) {
          if ( sscanf( tagarg, "%lu %lu", &(im->xdim), &(im->ydim) ) != 2 ) {
            if ( _verbose_ ) {
              fprintf( stderr, "%s: unable to read image dim (NDims = %d)\n", proc, ndims );
              fprintf( stderr, "\t '%s' -> '%s'\n", string, tagarg );
            }
            return( -1 );
          }
          im->zdim = 1;
        }
        else if ( ndims == 3 ) {
          if ( sscanf( tagarg, "%lu %lu %lu", &im->xdim, &im->ydim, &im->zdim ) != 3 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to read image dim (NDims = %d)\n", proc, ndims );
            return( -1 );
          }
        }
        else {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read image dim (undefined NDims = %d)\n", proc, ndims );
          return( -1 );
        }
    }
    else if ( _testTag( string, "HeaderSize", &tagarg ) == 1 ) {
        if ( sscanf( tagarg, "%d", &headersize ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read headersize ('%s')\n", proc, tagarg );
          return( -1 );
        }
    }
    else if ( _testTag( string, "Modality", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "SequenceID", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ElementMin", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ElementMax", &tagarg ) == 1 ) {
      ;
    }
    else if ( _testTag( string, "ElementNumberOfChannels", &tagarg ) == 1 ) {
        if ( sscanf( tagarg, "%u", &im->vdim ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to read ElementNumberOfChannels ('%s')\n", proc, tagarg );
          return( -1 );
        }
    }
    else if ( _testTag( string, "ElementSize", &tagarg ) == 1 ) {
      if ( elementspacing != 1 ) {
          if ( ndims == 1 ) {
            if ( sscanf( tagarg, "%lf", &im->vx ) != 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
              return( -1 );
            }
            im->vy = im->vz = 1.0;
          }
          else if ( ndims == 2 ) {
            if ( sscanf( tagarg, "%lf %lf", &im->vx, &im->vy ) != 2 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
              return( -1 );
            }
            im->vz = 1.0;
          }
          else if ( ndims == 3 ) {
            if ( sscanf( tagarg, "%lf %lf %lf", &im->vx, &im->vy, &im->vz ) != 3 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to read voxel size (NDims = %d)\n", proc, ndims );
              return( -1 );
            }
          }
          else {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to read voxel size (undefined NDims = %d)\n", proc, ndims );
            return( -1 );
          }
      }
    }
    else if ( _testTag( string, "ElementType", &tagarg ) == 1 ) {
      if ( (_testTag( tagarg, "MET_CHAR", (char**)NULL ) == 1 && tagarg[8] == '\0')
           || (_testTag( tagarg, "MET_CHAR_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
      }
      else if ( (_testTag( tagarg, "MET_UCHAR", (char**)NULL ) == 1 && tagarg[9] == '\0')
                || (_testTag( tagarg, "MET_UCHAR_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 1;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
      }
      else if ( (_testTag( tagarg, "MET_SHORT", (char**)NULL ) == 1 && tagarg[9] == '\0')
                || (_testTag( tagarg, "MET_SHORT_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
      }
      else if ( (_testTag( tagarg, "MET_USHORT", (char**)NULL ) == 1 && tagarg[10] == '\0')
                || (_testTag( tagarg, "MET_USHORT_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 2;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
      }
      else if ( (_testTag( tagarg, "MET_INT", (char**)NULL ) == 1 && tagarg[7] == '\0')
                || (_testTag( tagarg, "MET_INT_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 4;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
      }
      else if ( (_testTag( tagarg, "MET_UINT", (char**)NULL ) == 1 && tagarg[8] == '\0')
                || (_testTag( tagarg, "MET_UINT_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 4;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
      }
      else if ( (_testTag( tagarg, "MET_LONG", (char**)NULL ) == 1 && tagarg[8] == '\0')
                || (_testTag( tagarg, "MET_LONG_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 8;
        im->wordKind = WK_FIXED;
        im->sign = SGN_SIGNED;
      }
      else if ( (_testTag( tagarg, "MET_ULONG", (char**)NULL ) == 1 && tagarg[9] == '\0')
                || (_testTag( tagarg, "MET_ULONG_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 8;
        im->wordKind = WK_FIXED;
        im->sign = SGN_UNSIGNED;
      }
      else if ( (_testTag( tagarg, "MET_FLOAT", (char**)NULL ) == 1 && tagarg[9] == '\0')
                || (_testTag( tagarg, "MET_FLOAT_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 4;
        im->wordKind = WK_FLOAT;
        im->sign = SGN_SIGNED;
      }
      else if ( (_testTag( tagarg, "MET_DOUBLE", (char**)NULL ) == 1 && tagarg[10] == '\0')
                || (_testTag( tagarg, "MET_DOUBLE_ARRAY", (char**)NULL ) == 1) ) {
        im->wdim = 8;
        im->wordKind = WK_FLOAT;
        im->sign = SGN_SIGNED;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read ElementType ('%s')\n", proc, tagarg );
        return( -1 );
      }
    }
    else if ( _testTag( string, "ElementDataFile", &tagarg ) == 1 ) {

      /* - Name of the file to be loaded
       * - A printf-style string followed by the min, max, and step values to be used
       *   to pass an argument to the string to create list of file names to be loaded
       *   (must be (N-1)D blocks of data per file).
       * - LIST [X] – This specifies that starting on the next line is a list of files
       *   (one filename per line) in which the data is stored.
       *   Each file (by default) contains an (N-1)D block of data.
       *   If a second argument is given, its first character must be a number that
       *   specifies the dimension of the data in each file.
       *   For example ElementDataFile = LIST 2D means that there will be a 2D block
       *   of data per file.
       * - LOCAL – Indicates that the data begins at the beginning of the next line.
       */

      size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

      if ( _testTag( tagarg, "LOCAL", (char**)NULL ) == 1 && tagarg[5] == '\0' ) {
        if( !im->data ) {
          im->data = (unsigned char *) ImageIO_alloc(size);
          if ( _debug_ >= 2 ) {
              fprintf( stderr, "%s: allocate %lu bytes at %p\n", proc, size, im->data );
          }
          if ( !im->data ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: image buffer allocation failed\n", proc );
            return( -2 );
          }
        }
        nread = ImageIO_read(im, im->data, size);
        if ( nread != size ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: image buffer reading failed\n", proc );
            return( -1 );
        }
        _swapImageData( im );
        return( 1 );
      }

      else if ( _testTag( tagarg, "LIST", (char**)NULL ) == 1 &&
                (tagarg[4] == ' ' || tagarg[4] == '\0') ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: 'LIST' not implemented yet for ElementDataFile\n", proc );
        return( -1 );
      }

      else if ( _testPrintfFormat( tagarg ) == 1 ) {
         if ( _verbose_ )
           fprintf( stderr, "%s: 'FORMAT' not implemented yet for ElementDataFile\n", proc );
         return( -1 );
      }

      else {
        sprintf( rawname, "%s", name );
        s = rawname;
        s += strlen( rawname );
        while( s > rawname && *s != '/' )
          s--;
        if ( *s == '/' ) s++;
        sprintf( s, "%s", tagarg );

        ImageIO_close(im);
        _openReadImage( im, rawname );

        if( !im->data ) {
          im->data = (unsigned char *) ImageIO_alloc(size);
          if ( _debug_ >= 2 ) {
              fprintf( stderr, "%s: allocate %lu bytes at %p\n", proc, size, im->data );
          }
          if ( !im->data ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: image buffer allocation failed\n", proc );
            return( -2 );
          }
        }

        /* compute offset
         */
        if ( headersize == -1 ) {
          headersize = size - getFileSize( rawname );
        }
        /* skip offset
         */
        if ( headersize > 0 ) {
          char *tmp;
          tmp = (char*)malloc( headersize * sizeof(char) );
          if ( tmp == (char*)NULL ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to allocate auxiliary buffer\n", proc );
            return( -1 );
          }
          nread = ImageIO_read( im, tmp, headersize );
          free( tmp );
        }

        nread = ImageIO_read(im, im->data, size);
        if ( nread != size ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: image buffer reading failed\n", proc );
            return( -1 );
        }
        _swapImageData( im );
        return( 1 );

      }
    }
    /******************************
     * MetaImage Tags
     ******************************/
    else if ( _testTag( string, "CompressedData", &tagarg ) == 1 && string[14] == ' ' ) {
      if ( _testTag( tagarg, "True", (char**)NULL ) == 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: can not handle compressed data\n", proc );
        return( -1 );
      }
    }

    else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown tag '%s'\n", proc, string );
    }




  } while( retfgetns != (char*)NULL );



  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



static int _writeMetaImageHeader( char *name, const _image *im )
{
  char *proc = "_writeMetaImageHeader";
  char buf[1024];
  int length;

  if ( im->openMode == OM_CLOSE ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no open file\n ", proc );
    return( -1 );
  }

  sprintf( buf, "ObjectType = Image\n" );
  /* write 3D image as 3D image
   */
  /*
  if ( im->zdim == 1 ) {
    sprintf( buf+strlen(buf), "NDims = 2\n" );
  }
  else {
    sprintf( buf+strlen(buf), "NDims = 3\n" );
  }
  */
  sprintf( buf+strlen(buf), "NDims = 3\n" );

  if ( im->vdim > 1 )
    sprintf( buf+strlen(buf), "ElementNumberOfChannels = %d\n", im->vdim );

  sprintf( buf+strlen(buf), "BinaryData = True\n" );

  switch( _getEndianness() ) {
  default :
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown endianness\n ", proc );
      return( -1 );
  case END_LITTLE :
      sprintf( buf+strlen(buf), "BinaryDataByteOrderMSB = False\n" );
      break;
  case END_BIG :
      sprintf( buf+strlen(buf), "BinaryDataByteOrderMSB = True\n" );
      break;
  }

  sprintf( buf+strlen(buf), "CompressedData = False\n" );

  /* write 3D image as 3D image
   */
  /*
  if ( im->zdim == 1 ) {
    sprintf( buf+strlen(buf), "DimSize = %lu %lu\n", im->xdim, im->ydim );
    sprintf( buf+strlen(buf), "ElementSize = %f %f\n", im->vx, im->vy );
    sprintf( buf+strlen(buf), "ElementSpacing = %f %f\n", im->vx, im->vy );
  }
  else {
    sprintf( buf+strlen(buf), "DimSize = %lu %lu %lu\n", im->xdim, im->ydim, im->zdim );
    sprintf( buf+strlen(buf), "ElementSize = %f %f %f\n", im->vx, im->vy, im->vz );
    sprintf( buf+strlen(buf), "ElementSpacing = %f %f %f\n", im->vx, im->vy, im->vz );
  }
  */
  sprintf( buf+strlen(buf), "DimSize = %lu %lu %lu\n", im->xdim, im->ydim, im->zdim );
  sprintf( buf+strlen(buf), "ElementSize = %f %f %f\n", im->vx, im->vy, im->vz );
  sprintf( buf+strlen(buf), "ElementSpacing = %f %f %f\n", im->vx, im->vy, im->vz );


  switch ( im->wordKind ) {
  case WK_FIXED :
    switch( im->sign ) {
    case SGN_UNSIGNED :
      if ( im->wdim  == sizeof( unsigned char ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_UCHAR\n" );
      }
      else if ( im->wdim  == sizeof( unsigned short int ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_USHORT\n" );
      }
      else if ( im->wdim  == sizeof( unsigned int ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_UINT\n" );
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown WK_FIXED UNSIGNED word dim (%d)\n ", proc, im->wdim );
        return( -1 );
      }
      break;
    case SGN_SIGNED :
      if ( im->wdim  == sizeof( char ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_CHAR\n" );
      }
      else if ( im->wdim  == sizeof( short int ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_SHORT\n" );
      }
      else if ( im->wdim  == sizeof( int ) ) {
        sprintf( buf+strlen(buf), "ElementType = MET_INT\n" );
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
      sprintf( buf+strlen(buf), "ElementType = MET_FLOAT\n" );
    }
    else if ( im->wdim  == sizeof( double ) ) {
      sprintf( buf+strlen(buf), "ElementType = MET_DOUBLE\n" );
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



  length = strlen( name );
  if ( !strcmp( name + length - 4, ".mha" ) || !strcmp( name + length - 7, ".mha.gz" ) ) {
    sprintf( buf+strlen(buf), "ElementDataFile = LOCAL\n" );
  }
  else if ( !strcmp( name + length - 4, ".mhd" ) ) {
      sprintf( buf+strlen(buf), "HeaderSize = -1\n" );
      sprintf( buf+strlen(buf), "ElementDataFile = %s", _BaseName( name ) );
      sprintf( buf+strlen(buf)-3, "raw\n" );
  }
  else {
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown extension of '%s'\n ", proc, name );
      return( -1 );
  }

  if ( ImageIO_write(im, buf, strlen(buf)) == (size_t)EOF ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write header\n ", proc );
      return( -1 );
  }

  return( 1 );
}





static int _writeMetaImageData( char *name, _image *im )
{
  char *proc = "_writeMetaImageData";
  unsigned long nwrt, size;
  int length;
  char buf[512];

  size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

  length = strlen( name );

  if ( !strcmp( name + length - 4, ".mhd" ) ) {
      ImageIO_close( im );
      im->fd = NULL;
      im->openMode = OM_CLOSE;
      sprintf( buf, "%s", name );
      sprintf( buf+strlen(buf)-3, "raw" );
      _openWriteImage( im, buf );
  }

  if ( _debug_ )
    fprintf( stderr, "%s: write %lu bytes\n", proc, size );

  nwrt = ImageIO_write( im, im->data, size);
  if ( nwrt != size ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: write %lu/%lu bytes\n", proc, nwrt, size );
    return( -1 );
  }

  return( 1 );
}





int writeMetaImageHeader( char *name, _image *im )
{
  char *proc = "writeMetaImageHeader";
  int res;

  _openWriteImage( im, name );

  if( !im->fd ) {
    if ( _verbose_ )
      fprintf(stderr, "%s: error: unable to open file \'%s\'\n", proc, name );
    return( ImageIO_OPENING );
  }

  res = _writeMetaImageHeader( name, im );
  if (res < 0) {
    if ( _verbose_ )
      fprintf(stderr, "%s: error: unable to write header of \'%s\'\n",
              proc, name);
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }

  ImageIO_close( im );
  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return ( res );
}





int writeMetaImage( char *name, _image *im )
{
  char *proc = "writeMetaImage";
  int res;

  _openWriteImage( im, name );

  if( !im->fd ) {
    if ( _verbose_ )
      fprintf(stderr, "%s: error: unable to open file \'%s\'\n", proc, name );
    return( ImageIO_OPENING );
  }

  res = _writeMetaImageHeader( name, im );
  if (res < 0) {
    if ( _verbose_ )
      fprintf(stderr, "%s: error: unable to write header of \'%s\'\n",
              proc, name);
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }

  res = _writeMetaImageData( name, im );
  if (res < 0) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error: unable to write data of \'%s\'\n",
               proc, name );
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }


  ImageIO_close( im );
  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return ( res );
}





/************************************************************
 *
 *
 *
 ************************************************************/



PTRIMAGE_FORMAT createMetaImageFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat = &testMetaImageHeader;
  f->readImageHeader = &readMetaImage;
  f->readImageData = NULL;
  f->writeImageHeader = &writeMetaImageHeader;
  f->writeImage = &writeMetaImage;

  strcpy(f->fileExtension,".mha,.mha.gz,.mhd,");
  strcpy(f->realName,"MetaImage");
  return f;
}
