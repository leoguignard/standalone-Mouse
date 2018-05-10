

#include <string.h>
#include <strings.h>

#include <nifti1_io.h>

#include <nifti.h>


#define NIFTI_LE_MAGIC "\000\000\001\134"
#define NIFTI_BE_MAGIC "\134\001\000\000"

static int _verbose_ = 1;
static int _debug_ = 0;





/************************************************************
 *
 *
 *
 ************************************************************/



int testNiftiHeader( char *magic, const char *name __attribute__ ((unused)) )
{
  if( !memcmp(magic,NIFTI_LE_MAGIC,4) ||
      !memcmp(magic,NIFTI_BE_MAGIC,4))
    return 0;
  else
    return -1;
}





/************************************************************
 *
 *
 *
 ************************************************************/

static int _niftiImageToImageIOImage( nifti_image *nim, _image *im )
{
  char *proc = "_niftiImageToImageIOImage";
  int i, j;

  /* data array dimension
   */
  if ( nim->ndim > 5 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: does not handle ndim=%d\n", proc, nim->ndim );
    return( -1 );
  }
  im->xdim = nim->nx;
  im->ydim = nim->ny;
  im->zdim = nim->nz;

  if ( nim->nt > 1 && nim->nu > 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: does not handle both nt>1 and nu>1\n",
               proc );
    return( -1 );
  }
  if ( nim->nv > 1 || nim->nw > 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: does not handle (nv,nw)=(%d,%d)\n",
               proc, nim->nv, nim->nw );
    return( -1 );
  }

  /* a vectorial image can be either vectorial (nu>1)
   * or temporal (nt>1)
   */
  if ( nim->nt > 1 )
    im->vdim = nim->nt;
  if ( nim->nu > 1 )
    im->vdim = nim->nu;


  /* data type
   */
  switch( nim->datatype ) {
  default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such datatype (%d) not handled yet\n",
                 proc, nim->datatype );
      return( -1 );
  case DT_INT8 :
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case DT_UINT8 :
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case DT_INT16 :
      im->wdim = 2;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case DT_UINT16 :
      im->wdim = 2;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case DT_INT32 :
      im->wdim = 4;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case DT_UINT32 :
      im->wdim = 4;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case DT_INT64 :
      im->wdim = 8;
      im->wordKind = WK_FIXED;
      im->sign = SGN_SIGNED;
      break;
  case DT_UINT64 :
      im->wdim = 8;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  case DT_FLOAT :
      im->wdim = 4;
      im->wordKind = WK_FLOAT;
      im->sign = SGN_SIGNED;
      break;
  case DT_DOUBLE :
      im->wdim = 8;
      im->wordKind = WK_FLOAT;
      im->sign = SGN_SIGNED;
      break;
  case DT_RGB :
      im->vdim = 3;
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      break;
  }

  /* voxel size
   */
  im->vx = nim->dx;
  im->vy = nim->dy;
  im->vz = nim->dz;

  /* matrix issued from the qform
   */
  if ( nim->qform_code != NIFTI_XFORM_UNKNOWN ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      im->qform_toreal[i][j] = nim->qto_xyz.m[i][j];
    im->qform_code = nim->qform_code;
  }

  /* matrix issued from the sform
   */
  if ( nim->sform_code != NIFTI_XFORM_UNKNOWN ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      im->sform_toreal[i][j] = nim->sto_xyz.m[i][j];
    im->sform_code = nim->sform_code;
  }


  if ( _debug_ ) {
    switch( nim->byteorder ) {
    default :
      fprintf( stderr, "%s: unknown byte order\n", proc );
      break;
    case 1 :
      fprintf( stderr, "%s: byte order is Least Significant Bit (little Endian)\n", proc );
      break;
    case 2 :
      fprintf( stderr, "%s: byte order is Most Significant Bit (big Endian)\n", proc );
      break;
    }
  }

  return( 1 );
}



static int _niftDataToImageIOData( nifti_image *nim, _image *im )
{
  char *proc = "_niftDataToImageIOData";
  size_t x, y, z, v;

  if ( im->vdim == 1 ) {
    (void)memcpy( im->data, nim->data, nim->nvox * nim->nbyper );
    return( 1 );
  }

#define NIFTIDATATOIMAGEIODATA( TYPE ) { \
  TYPE *niftiBuf = (TYPE*)nim->data;    \
  TYPE *ioBuf = (TYPE*)im->data;        \
  TYPE *tmpBuf;                         \
  for ( v=0; v<im->vdim; v++ ) {        \
    tmpBuf = ioBuf;                     \
    tmpBuf += v;                        \
    for ( z=0; z<im->zdim; z++ )        \
    for ( y=0; y<im->ydim; y++ )        \
    for ( x=0; x<im->xdim; x++, niftiBuf++, tmpBuf+=im->vdim ) { \
      *tmpBuf = *niftiBuf;              \
    }                                   \
  }                                     \
}

  switch( nim->datatype ) {
  default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such datatype (%d) not handled yet\n",
                 proc, nim->datatype );
      return( -1 );
  case DT_INT8 :
      NIFTIDATATOIMAGEIODATA( char );
      break;
  case DT_RGB :
  case DT_UINT8 :
      NIFTIDATATOIMAGEIODATA( unsigned char );
      break;
  case DT_INT16 :
      NIFTIDATATOIMAGEIODATA( short int );
      break;
  case DT_UINT16 :
      NIFTIDATATOIMAGEIODATA( unsigned short int );
      break;
  case DT_INT32 :
      NIFTIDATATOIMAGEIODATA( int );
      break;
  case DT_UINT32 :
      NIFTIDATATOIMAGEIODATA( unsigned int );
      break;
  case DT_INT64 :
      NIFTIDATATOIMAGEIODATA( long int );
      break;
  case DT_UINT64 :
      NIFTIDATATOIMAGEIODATA( unsigned long int );
      break;
  case DT_FLOAT :
      NIFTIDATATOIMAGEIODATA( float );
      break;
  case DT_DOUBLE :
      NIFTIDATATOIMAGEIODATA( double );
      break;
  }

  return( 1 );
}



int readNiftiImage( const char *name, _image *im )
{
  char *proc = "readNiftiImage";
  nifti_image *nim;

  ImageIO_close(im);
  im->fd = NULL;

  /* nifti_image *nifti_image_read( const char *hname , int read_data )
   *    - The data buffer will be byteswapped if necessary.
   *    - The data buffer will not be scaled.
   *    - The data buffer is allocated with calloc().
   * \param hname filename of the nifti dataset
   * \param read_data Flag, true=read data blob, false=don't read blob.
   * \return A pointer to the nifti_image data structure.
   */
  nim = nifti_image_read( name, 1 ) ;
  if ( _niftiImageToImageIOImage( nim, im ) != 1 ) {
    nifti_image_free( nim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate nifti header\n", proc );
    return( -1 );
  }

  im->data = (unsigned char *) ImageIO_alloc( nim->nvox * nim->nbyper );
  if ( im->data == (unsigned char *)NULL ) {
    nifti_image_free( nim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate data buffer \n", proc );
    return( -1 );
  }

  if ( _niftDataToImageIOData( nim, im ) != 1 ) {
    ImageIO_free( im->data );
    nifti_image_free( nim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate nifti data\n", proc );
    return( -1 );
  }

  nifti_image_free( nim );

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/

static int _imageIOImageToNiftiImage( _image *im, nifti_image *nim )
{
  char *proc = "_imageIOImageToNiftiImage";
  int i, j;

  nim->nx = nim->dim[1] = (int)im->xdim;
  nim->ny = nim->dim[2] = (int)im->ydim;
  nim->nz = nim->dim[3] = (int)im->zdim;
  nim->nt = nim->dim[4] = 1;
  nim->nu = nim->dim[5] = (int)im->vdim;
  nim->nv = nim->dim[6] = 1;
  nim->nw = nim->dim[7] = 1;
  if ( nim->nu > 1 ) {
    nim->ndim = nim->dim[0] = 5;
    nim->intent_code = NIFTI_INTENT_VECTOR;
  }
  /* 2D are treated as 3D images
   */
  else {
    nim->ndim = nim->dim[0] = 3;
  }
  /* to deal differently with 2D and
   * 3D images
   */
  /*
  else if ( nim->nz > 1 ) {
    nim->ndim = nim->dim[0] = 3;
  }
  else {
    nim->ndim = nim->dim[0] = 2;
  }
  */

  switch ( im->wordKind ) {
  case WK_FIXED :
    switch( im->sign ) {
    case SGN_UNSIGNED :
      if ( im->wdim  == sizeof( unsigned char ) ) {
        if ( im->vdim == 3 ) {
          nim->datatype = DT_RGB;
          nim->nt = nim->dim[4] = 1;
          nim->ndim = nim->dim[0] = ( nim->nz > 1 ) ? 3 : 2;
          nim->ndim = nim->dim[0] = 3;
          nim->nbyper = 3 * sizeof( unsigned char );
        }
        else {
          nim->datatype = DT_UINT8;
          nim->nbyper = sizeof( unsigned char );
        }
      }
      else if ( im->wdim  == sizeof( unsigned short int ) ) {
          nim->datatype = DT_UINT16;
          nim->nbyper = sizeof( unsigned short int );
      }
      else if ( im->wdim  == sizeof( unsigned int ) ) {
          nim->datatype = DT_UINT32;
          nim->nbyper = sizeof( unsigned int );
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unknown WK_FIXED UNSIGNED word dim (%d)\n ", proc, im->wdim );
        return( -1 );
      }
      break;
    case SGN_SIGNED :
      if ( im->wdim  == sizeof( char ) ) {
          nim->datatype = DT_INT8;
          nim->nbyper = sizeof( char );
      }
      else if ( im->wdim  == sizeof( short int ) ) {
          nim->datatype = DT_INT16;
          nim->nbyper = sizeof( short int );
      }
      else if ( im->wdim  == sizeof( int ) ) {
          nim->datatype = DT_INT32;
          nim->nbyper = sizeof( int );
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
        nim->datatype = DT_FLOAT;
        nim->nbyper = sizeof( float );
    }
    else if ( im->wdim  == sizeof( double ) ) {
        nim->datatype = DT_DOUBLE;
        nim->nbyper = sizeof( double );
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

  nim->nvox = (size_t)nim->nx * (size_t)nim->ny * (size_t)nim->nz * (size_t)nim->nt * (size_t)nim->nu;

  nim->pixdim[0] = 1.0;
  nim->dx = nim->pixdim[1] = im->vx;
  nim->dy = nim->pixdim[2] = im->vy;
  nim->dz = nim->pixdim[3] = im->vz;
  nim->dt = nim->pixdim[4] = 1.0;
  nim->du = nim->pixdim[5] = 1.0;
  nim->dv = nim->pixdim[6] = 1.0;
  nim->dw = nim->pixdim[7] = 1.0;

  nim->xyz_units = NIFTI_UNITS_MM;

  if ( 1 ) {
    nim->scl_slope = 1.0;
    nim->scl_inter = 0.0;
  }

  /* matrix issued from the qform
   */
  if ( im->qform_code || im->t_is_set ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      nim->qto_xyz.m[i][j] = im->qform_toreal[i][j];
    nifti_mat44_to_quatern( nim->qto_xyz,
                            &(nim->quatern_b), &(nim->quatern_c), &(nim->quatern_d),
                            &(nim->qoffset_x), &(nim->qoffset_y), &(nim->qoffset_z),
                            &(nim->dx), &(nim->dy), &(nim->dz),
                            &(nim->qfac) );
    if ( 0 ) {
      nim->pixdim[1] = nim->dx;
      nim->pixdim[2] = nim->dy;
      nim->pixdim[3] = nim->dz;
    }
    if ( im->qform_code > 0 )
      nim->qform_code = im->qform_code;
    else
      nim->qform_code = 5;
  }

  if ( im->sform_code ) {
    for ( i=0; i<4; i++ )
    for ( j=0; j<4; j++ )
      nim->sto_xyz.m[i][j] = im->sform_toreal[i][j];
    nim->sform_code = im->sform_code;
  }

  return( 1 );
}



static int _imageIODataToNiftiData( _image *im, nifti_image *nim )
{
  char *proc = "_imageIODataToNiftiData";
  size_t x, y, z, v;

  if ( im->vdim == 1 ) {
    (void)memcpy( nim->data, im->data, nim->nvox * nim->nbyper );
    return( 1 );
  }

#define IMAGEIODATATONIFTIDATA( TYPE ) { \
  TYPE *niftiBuf = (TYPE*)nim->data;    \
  TYPE *ioBuf = (TYPE*)im->data;        \
  TYPE *tmpBuf;                         \
  for ( v=0; v<im->vdim; v++ ) {        \
    tmpBuf = ioBuf;                     \
    tmpBuf += v;                        \
    for ( z=0; z<im->zdim; z++ )        \
    for ( y=0; y<im->ydim; y++ )        \
    for ( x=0; x<im->xdim; x++, niftiBuf++, tmpBuf+=im->vdim ) { \
      *niftiBuf = *tmpBuf;              \
    }                                   \
  }                                     \
}

  switch( nim->datatype ) {
  default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such datatype (%d) not handled yet\n",
                 proc, nim->datatype );
      return( -1 );
  case DT_INT8 :
      IMAGEIODATATONIFTIDATA( char );
      break;
  case DT_RGB :
  case DT_UINT8 :
      IMAGEIODATATONIFTIDATA( unsigned char );
      break;
  case DT_INT16 :
      IMAGEIODATATONIFTIDATA( short int );
      break;
  case DT_UINT16 :
      IMAGEIODATATONIFTIDATA( unsigned short int );
      break;
  case DT_INT32 :
      IMAGEIODATATONIFTIDATA( int );
      break;
  case DT_UINT32 :
      IMAGEIODATATONIFTIDATA( unsigned int );
      break;
  case DT_INT64 :
      IMAGEIODATATONIFTIDATA( long int );
      break;
  case DT_UINT64 :
      IMAGEIODATATONIFTIDATA( unsigned long int );
      break;
  case DT_FLOAT :
      IMAGEIODATATONIFTIDATA( float );
      break;
  case DT_DOUBLE :
      IMAGEIODATATONIFTIDATA( double );
      break;
  }

  return( 1 );
}



static int _buildNames( char *name, nifti_image *nim,
                        char **hdrname, char **imgname )
{
    char *proc = "_buildNames";
    int l;

    l = strlen( name );


    if ( strcmp( name + l - 4, ".nia" ) == 0 ) {
        nim->fname = name;
        nim->iname = name;
        nim->nifti_type = 3;
    }
    else if ( (strcmp( name + l - 4, ".nii" ) == 0)
         || (strcmp( name + l - 7, ".nii.gz" ) == 0) ) {
      nim->fname = name;
      nim->iname = name;
      nim->nifti_type = 1;
    }
    else if ( strcmp( name + l - 4, ".hdr" ) == 0 ) {
      *imgname = (char*)malloc( (l+1)*sizeof(char) );
      sprintf( *imgname, "%s", name );
      sprintf( *imgname + l - 4, ".img" );
      (*imgname)[l] = '\0';
      nim->fname = name;
      nim->iname = *imgname;
      nim->nifti_type = 0;
    }
    else if ( strcmp( name + l - 4, ".img" ) == 0 ) {
      *hdrname = (char*)malloc( (l+1)*sizeof(char) );
      sprintf( *hdrname, "%s", name );
      sprintf( *hdrname + l - 4, ".hdr" );
      (*hdrname)[l] = '\0';
      nim->fname = *hdrname;
      nim->iname = name;
      nim->nifti_type = 0;
    }
    else if ( strcmp( name + l - 7, ".img.gz" ) == 0 ) {
      *hdrname = (char*)malloc( (l+1)*sizeof(char) );
      sprintf( *hdrname, "%s", name );
      sprintf( *hdrname + l - 7, ".hdr" );
      (*hdrname)[l-3] = '\0';
      nim->fname = *hdrname;
      nim->iname = name;
      nim->nifti_type = 0;
    }
    else {
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown extension\n", proc );
      return( -1 );
    }

    return( 1 );
}





int writeNiftiImageHeader( char *name, _image *im )
{
  char *proc = "writeNiftiImageHeader";
  nifti_image nim;
  char *hdrname = (char*)NULL;
  char *imgname = (char*)NULL;
  znzFile fp;

  bzero( &nim, sizeof(nifti_image) );
  if ( _imageIOImageToNiftiImage( im, &nim) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate into nifti image\n", proc );
    return( -1 );
  }

  if ( _buildNames( name, &nim, &hdrname, &imgname ) != 1 ) {
    if ( hdrname != (char*)NULL ) free( hdrname );
    if ( imgname != (char*)NULL ) free( imgname );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to build names from '%s'\n", proc, name );
    }
    return( -1 );
  }

  nim.data = NULL;

  fp = nifti_image_write_hdr_img( &nim, 0, "wb" );
  if ( fp ) {
     free(fp);
  }

  if ( hdrname != (char*)NULL ) free( hdrname );
  if ( imgname != (char*)NULL ) free( imgname );

  return( 1 );
}





int writeNiftiImage( char *name, _image *im )
{
  char *proc = "writeNiftiImage";
  nifti_image nim;
  char *hdrname = (char*)NULL;
  char *imgname = (char*)NULL;

  bzero( &nim, sizeof(nifti_image) );
  if ( _imageIOImageToNiftiImage( im, &nim) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate into nifti image\n", proc );
    return( -1 );
  }

  if ( _buildNames( name, &nim, &hdrname, &imgname ) != 1 ) {
    if ( hdrname != (char*)NULL ) free( hdrname );
    if ( imgname != (char*)NULL ) free( imgname );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to build names from '%s'\n", proc, name );
    }
    return( -1 );
  }

  nim.data = (unsigned char *) ImageIO_alloc( nim.nvox * nim.nbyper );
  if ( nim.data == (unsigned char *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate data buffer \n", proc );
    return( -1 );
  }

  if ( _imageIODataToNiftiData( im, &nim ) != 1 ) {
    ImageIO_free( nim.data );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to translate to nifti data\n", proc );
    return( -1 );
  }

  if ( 0 ) {
    fprintf( stderr, "%s: will write\n", proc );
    fprintf( stderr, "   header '%s'\n", nim.fname );
    fprintf( stderr, "   data   '%s'\n", nim.iname );
  }

  nifti_image_write( &nim );

  ImageIO_free( nim.data );

  if ( hdrname != (char*)NULL ) free( hdrname );
  if ( imgname != (char*)NULL ) free( imgname );

  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



PTRIMAGE_FORMAT createNiftiFormat()
{
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));
  _initImageFormat( f );

  f->testImageFormat = &testNiftiHeader;
  f->readImageHeader = &readNiftiImage;
  f->readImageData = NULL;
  f->writeImageHeader = &writeNiftiImageHeader;
  f->writeImage = &writeNiftiImage;

  strcpy(f->fileExtension,".nia,.nii,.nii.gz,.hdr,.img,.img.gz");
  strcpy(f->realName,"Nifti");
  return f;
}








