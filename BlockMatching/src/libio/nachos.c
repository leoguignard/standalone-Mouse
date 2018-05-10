/* ad-hoc reader/writer for Nachos raw data
 */

#include <string.h>

#include <ImageIO.h>


static int _debug_ = 1;


static char *endofline_unix="\n";
/* static char *endofline_window="\r\n"; */








typedef struct _list {
  int n_data;
  int n_allocated_data;
  float *data;
} _list;

static void _initList( _list *l )
{
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->data = (float*)NULL;
}

static void _freeList( _list *l )
{
  if ( l->data != (float*)NULL )
    free( l->data );
  _initList( l );
}

static int _values_to_be_allocated_ = 100;

static int _addValueToList( _list *l, float v )
{
  char *proc = "_addValueToList";
  int s =  l->n_allocated_data;
  float *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _values_to_be_allocated_;
    data = (float*)malloc( s * sizeof(float) );
    if ( data == (float*)NULL ) {
      if ( _debug_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(float) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = v;
  l->n_data ++;
  return( 1 );
}



typedef struct _listList {
  int n_data;
  int n_allocated_data;
  _list *data;
} _listList;

static void _initListList( _listList *l )
{
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->data = (_list*)NULL;
}

static void _freeListList( _listList *l )
{
  int i;
  if ( l->data != (_list*)NULL ) {
    for ( i=0; i<l->n_data; i++ ) {
      _freeList( &(l->data[i]) );
    }
    free( l->data );
  }
  _initListList( l );
}

static int _lists_to_be_allocated_ = 100;

static _list *_getList( _listList *l )
{
  char *proc = "_getList";
  int s =  l->n_allocated_data;
  _list *data;
  int i;

  if ( l->n_data == l->n_allocated_data ) {
    s += _lists_to_be_allocated_;
    data = (_list*)malloc( s * sizeof(_list) );
    if ( data == (_list*)NULL ) {
      if ( _debug_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( (_list *)NULL );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_list) );
      free( l->data );
    }

    for ( i=l->n_allocated_data; i<s; i++ )
      _initList( &(data[i]) );

    l->n_allocated_data = s;
    l->data = data;
  }

  l->n_data ++;
  return( &(l->data[l->n_data-1]) );
}





typedef struct _str {
  int n_data;
  int n_allocated_data;
  char *data;
} _str;

static void _initStr( _str *s )
{
  s->n_data = 0;
  s->n_allocated_data = 0;
  s->data = (char*)NULL;
}

static void _freeStr( _str *s )
{
  if ( s->data != (char*)NULL )
    free( s->data );
  _initStr( s );
}



static int _chars_to_be_allocated_ = 1000;



static int _reAllocStr( _str *l, int size )
{
  char *proc = "_reAllocStr";
  int s =  l->n_allocated_data;
  char *data;

  s += size;
  data = (char*)malloc( s * sizeof(char) );
  if ( data == (char*)NULL ) {
    if ( _debug_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }
  memset( data, 0, s*sizeof(char) );

  if ( l->data != (char*)NULL && l->n_allocated_data > 0 ) {
    (void)memcpy( data, l->data, l->n_allocated_data*sizeof(char) );
    free( l->data );
  }
  l->n_allocated_data = s;
  l->data = data;

  return( 1 );
}



static int _isLineDone( char *s, int l )
{
  int i;

  if ( s == (char*)NULL ) return( -1 );
  for ( i=0; i<l; i++ ) {
    if ( s[i] == '\n' ) return( 1 );
  }
  return( -1 );
}






static int _readStr( _image *im, _str *str )
{
  char *proc = "_readStr";
  char *ret, *buf;
  int plength, nlength, length;

  if ( str->n_allocated_data > 0 ) {
    memset( str->data, 0, str->n_allocated_data*sizeof(char) );
  }
  else {
    if ( _reAllocStr( str, _chars_to_be_allocated_ ) != 1 ) {
      if ( _debug_ )
         fprintf( stderr, "%s: allocation error\n", proc );
        return( -1 );
    }
  }
  str->n_data = 0;


  buf = str->data;
  length = str->n_allocated_data;

  while ( 1 ) {

    ret = ImageIO_gets( im, buf, length );
    if ( ret == (char*)NULL ) break;

    if ( _isLineDone( buf, length ) == 1 ) {
      break;
    }

    plength = str->n_allocated_data;
    if ( _reAllocStr( str, _chars_to_be_allocated_ ) != 1 ) {
      if ( _debug_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    nlength = str->n_allocated_data;

    buf = &(str->data[plength-1]);
    length = nlength - plength + 1;
  }

  return( 1 );
}



static int _translateLine( _list *l, _str *str )
{
  char *proc = "_translateLine";
  char *s = str->data;
  float value;
  int i;

  while( 1 ) {
    if ( *s == '\0' || *s == '\n' || (s[0] == '\r' && s[1] == '\n') )
      break;
    if ( *s == '-' || *s == '.' || (*s >= '0' && *s <= '9') ) {
      if ( sscanf( s, "%f", &value ) != 1 ) {
        if ( _debug_ )
          fprintf( stderr, "%s: error when reading value\n", proc );
        return( -1 );
      }
      if ( _addValueToList( l, value ) != 1 ) {
        if ( _debug_ )
          fprintf( stderr, "%s: error when adding value to list\n", proc );
        return( -1 );
      }
      i++;
      while ( *s == '-' || *s == '.' || (*s >= '0' && *s <= '9') )
        s++;
      if ( *s == ';' )
        s++;
      else {
        if ( _debug_ )
          fprintf( stderr, "%s: weird case\n", proc );
        break;
      }
    }
    else {
      break;
    }
  }
  return( i );
}




int readNachosImage( const char *name __attribute__ ((unused)), _image *im)
{
  char *proc = "readNachosImage";
  _listList listList;
  _list *list;
  _str str;
  int y=0;
  float *buf;
  int npoints;

  _initListList( &listList );


  _initStr( &str );

  /* reading the file
   */
  while( 1 ) {
    /* reading a line
     */
    if ( _readStr( im, &str ) != 1 ) {
      _freeListList( &listList );
      _freeStr( &str );
      if ( _debug_ )
        fprintf( stderr, "%s: error when reading line\n", proc );
      return( -1 );
    }

    list = _getList( &listList );

    npoints = _translateLine( list, &str );
    if ( npoints == -1 ) {
      _freeListList( &listList );
      _freeStr( &str );
      if ( _debug_ )
        fprintf( stderr, "%s: error when translating line\n", proc );
      return( -1 );
    }

    if ( 0 ) fprintf( stderr, "%s: line #%d has %d values\n", proc, y, list->n_data );

    if ( list->n_data == 0 ) {
      _freeList( list );
      listList.n_data --;
      break;
    }
    y++;
  }

  _freeStr( &str );

  /* verifying some stuff
   */
  if ( listList.data[0].n_data <= 0 ) {
    _freeListList( &listList );
    if ( _debug_ )
      fprintf( stderr, "%s: empty first line ?!\n", proc );
    return( -1 );
  }

  for ( y=1; y<listList.n_data; y++ ) {
    if ( listList.data[y].n_data != listList.data[0].n_data ) {
      _freeListList( &listList );
      if ( _debug_ )
        fprintf( stderr, "%s: line #%d has a different number of points than line #0\n", proc, y );
      return( -1 );
    }
  }

  if ( 0 ) fprintf( stderr, "%s: dimy=%d, dimx=%d\n", proc, listList.n_data, listList.data[0].n_data );

  im->wordKind = WK_FLOAT;
  im->wdim = sizeof(float);
  im->zdim = 1;
  im->vdim = 1;

  im->xdim = listList.data[0].n_data;
  im->ydim = listList.n_data;

  im->data = ImageIO_alloc( im->xdim * im->ydim * sizeof(float) );

  buf = (float*)im->data;
  for ( y=0; y<listList.n_data; y++ ) {
      memcpy( &(buf[y*im->xdim]), listList.data[y].data, im->xdim * sizeof(float) );
  }

  _freeListList( &listList );

  return( 1 );
}










int _convertToFloat( float *res, const _image *im )
{
  char *proc = "_convertToFloat";
  unsigned long i, size;

  size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

  switch( im->wordKind ) {
  case WK_FLOAT:
    switch( im->wdim ) {
    case 4 :
      {
        float *buf = (float*)im->data;
        for ( i=0; i<size; i++ ) *res++ = *buf++;
      }
      break;
    default :
      if ( _debug_ )
        fprintf( stderr, "%s: such word dim not handled yet\n", proc  );
      return( -1 );
    }

    break;
  case WK_FIXED:
    switch(im->sign) {
    case SGN_SIGNED:
      switch( im->wdim ) {
      case 1 :
        {
          char *buf = (char*)im->data;
          for ( i=0; i<size; i++ ) *res++ = (float)*buf++;
        }
        break;
      default :
        if ( _debug_ )
          fprintf( stderr, "%s: such word dim not handled yet\n", proc  );
        return( -1 );
      }
      break;
    case SGN_UNSIGNED :
      switch( im->wdim ) {
      case 1 :
        {
          unsigned char *buf = (unsigned char*)im->data;
          for ( i=0; i<size; i++ ) *res++ = (float)*buf++;
        }
        break;
      default :
        if ( _debug_ )
          fprintf( stderr, "%s: such word dim not handled yet\n", proc  );
        return( -1 );
      }
      break;
    default :
      if ( _debug_ )
        fprintf( stderr, "%s: weird sign\n", proc  );
      return( -1 );
    }
    break;
  default :
    if ( _debug_ )
      fprintf( stderr, "%s: such word kind not handled yet\n", proc  );
    return( -1 );
  }

  return( 1 );
}



char *_advanceToNull( char *buf, int *dl )
{
  *dl = 0;
  while( *buf != '\0' ) {
    buf ++;
    (*dl) ++;
  }
  return( buf );
}





/* Writes the given image body in an already opened file.*/
int _writeNachosData(const _image *im) {
  char *proc = "_writeNachosData";
  unsigned long size, nwrt;
  float *floatBuf = (float*)NULL;
  char *buf;
  int dl;
  unsigned long i, x, y, z;
  _str str;
  char *endofline = endofline_unix;

  _initStr( &str );

  if ( im->vdim > 1 ) {
      if ( _debug_ )
        fprintf( stderr, "%s: can not deal with vectorial data (vdim=%u)\n", proc, im->vdim  );
      return( -1 );
  }
  if ( im->zdim > 1 ) {
      if ( _debug_ )
        fprintf( stderr, "%s: can not deal with 3D data (zdim=%lu)\n", proc, im->zdim  );
      return( -1 );
  }

  size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;
  floatBuf = (float*)malloc( size * sizeof(float) );
  if ( floatBuf == (float*)NULL ) {
    if ( _debug_ )
      fprintf( stderr, "%s: allocation failed\n", proc  );
    return( -1 );
  }

  if ( _convertToFloat( floatBuf, im ) != 1 ) {
    free( floatBuf );
    if ( _debug_ )
      fprintf( stderr, "%s: conversion failed\n", proc  );
    return( -1 );
  }

  if ( _reAllocStr( &str, im->xdim * 12 ) != 1 ) {
    free( floatBuf );
    if ( _debug_ )
       fprintf( stderr, "%s: string allocation error\n", proc );
    return( -1 );
  }

  buf = str.data;
  str.n_data = 0;

  for ( i=0, z=0; z<im->zdim; z++ ) {
    for ( y=0; y<im->ydim; y++ ) {
      for (x=0; x<im->xdim*im->vdim; x++, i++ ) {

        buf = _advanceToNull( buf, &dl );
        str.n_data += dl;

        if ( str.n_data > str.n_allocated_data - 50 ) {
          nwrt = ImageIO_write( im, str.data, (unsigned long)str.n_data );
          if ( nwrt != (unsigned long)str.n_data  ) {
            _freeStr( &str );
            free( floatBuf );
            if ( _debug_ )
               fprintf( stderr, "%s: writing error\n", proc );
            return( -1 );
          }
          if ( 0 )
            fprintf( stderr, "%s: write %d bytes\n", proc, str.n_data );
          buf = str.data;
          str.n_data = 0;
        }

        sprintf( buf, "%f;", floatBuf[i] );
      }
      buf = _advanceToNull( buf, &dl );
      str.n_data += dl;
      sprintf( buf, "%s", endofline );
    }
  }

  buf = _advanceToNull( buf, &dl );
  str.n_data += dl;

  if ( str.n_data > 0 ) {
    nwrt = ImageIO_write( im, str.data, (unsigned long)str.n_data );
    if ( nwrt != (unsigned long)str.n_data  ) {
      _freeStr( &str );
      free( floatBuf );
      if ( _debug_ )
         fprintf( stderr, "%s: writing error\n", proc );
      return( -1 );
    }
    if ( 0 )
      fprintf( stderr, "%s: write %d bytes\n", proc, str.n_data );
  }

  _freeStr( &str );
  free( floatBuf );

  return( 1 );
}






int testNachosHeader( char *magic __attribute__ ((unused)),
                   const char *name __attribute__ ((unused)) )
{
    return -1;
}





int writeNachosImage(char *name,_image *im) {
  int res;


  _openWriteImage( im, name );

  if(!im->fd) {
    fprintf(stderr, "writeNachos: error: unable to open file \'%s\'\n", name );
    return ImageIO_OPENING;
  }

  
  res = _writeNachosData( im );
  if (res < 0) {
    fprintf(stderr, "writeNachos: error: unable to write data of \'%s\'\n",
            name);
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





PTRIMAGE_FORMAT createNachosFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat=&testNachosHeader;
  f->readImageHeader=&readNachosImage;
  f->readImageData=0;
  f->writeImage=&writeNachosImage;
  strcpy(f->fileExtension,".dat,.dat.gz,.nachos");
  strcpy(f->realName,"NachosData");
  return f;
}
