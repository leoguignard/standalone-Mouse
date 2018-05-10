/*************************************************************************
 * api-cropImage.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 13 sep 2016 14:48:06 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <extract.h>
#include <vtmalloc.h>

#include <api-cropImage.h>
#include <bal-transformation-copy.h>
#include <bal-transformation-inversion.h>



static int _verbose_ = 1;
static int _debug_ = 1;


static void _API_ParseParam_cropImage( char *str, lineCmdParamCropImage *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/

/* this one is kept for historical reasons,
 * ie Stracking compilation but should disappear
 */
int cropImage( char *theim_name,
               char *resim_name,
               char *real_transformation_name,
               char *voxel_transformation_name,
               char *template_image_name, /* template */
               int analyzeFiji,
               bal_integerPoint origin,
               bal_integerPoint dim,
               bal_integerPoint slice,
               int isDebug,
               int isVerbose )
{
  char *proc = "cropImage";
  char str[256], *s;

  _verbose_ = isVerbose;
  _debug_ = isDebug;

  if ( _verbose_ )
      fprintf( stderr, "Warning, '%s' is obsolete\n", proc );

  s = str;

  if ( slice.z >= 0 ) {
    sprintf( s, "-xy %d ", slice.z );
    s = &(str[strlen(str)]);
  }

  else if ( slice.y >= 0 ) {
    sprintf( s, "-xz %d ", slice.y );
    s = &(str[strlen(str)]);
  }

  else if ( slice.x >= 0 ) {
    sprintf( s, "-yz %d ", slice.x );
    s = &(str[strlen(str)]);
  }

  else {
    if ( origin.x >= 0 ) {
      sprintf( s, "-ix %d ", origin.x );
      s = &(str[strlen(str)]);
    }
    if ( origin.y >= 0 ) {
      sprintf( s, "-iy %d ", origin.y );
      s = &(str[strlen(str)]);
    }
    if ( origin.z >= 0 ) {
      sprintf( s, "-iz %d ", origin.z );
      s = &(str[strlen(str)]);
    }
    if ( dim.x >= 1 ) {
      sprintf( s, "-x %d ", dim.x );
      s = &(str[strlen(str)]);
    }
    if ( dim.y >= 1 ) {
      sprintf( s, "-y %d ", dim.y );
      s = &(str[strlen(str)]);
    }
    if ( dim.z >= 1 ) {
      sprintf( s, "-z %d ", dim.z );
      s = &(str[strlen(str)]);
    }
  }

  if ( analyzeFiji ) {
    sprintf( s, "-analyze-fiji " );
    s = &(str[strlen(str)]);
  }

  if ( API_INTERMEDIARY_cropImage( theim_name,
                                   resim_name,
                                   real_transformation_name,
                                   voxel_transformation_name,
                                   template_image_name,
                                   s, (char*)NULL ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: some error occurs\n", proc );
      return( -1 );
  }

  return( 0 );
}





int API_INTERMEDIARY_cropImage( char *theim_name,
                                char *resim_name,
                                char *real_transformation_name,
                                char *voxel_transformation_name,
                                char *template_image_name,
                                char *param_str_1, char *param_str_2 )

{
  char *proc = "API_INTERMEDIARY_cropImage";
  bal_image theim;
  bal_image resim;
  bal_image templateim;
  bal_transformation theTrsf;
  int theLeftCorner[3] = {0, 0, 0};
  int theRightCorner[3] = {0, 0, 0};

  lineCmdParamCropImage par;



  /* parameter initialization
   */
  API_InitParam_cropImage( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cropImage( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cropImage( param_str_2, &par );



  /***************************************************
   *
   * initialization
   *
   ***************************************************/

  BAL_InitImage( &theim, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &resim, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &templateim, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitTransformation( &theTrsf );

  if ( BAL_AllocTransformation( &theTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    fprintf( stderr, "%s: Unable to allocate result transformation \n", proc);
    return( -1 );
  }
  BAL_SetTransformationToIdentity( &theTrsf );
  theTrsf.transformation_unit = VOXEL_UNIT;



  /* reading input image
   */
  if ( BAL_ReadImage( &theim, theim_name, 0 ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      fprintf( stderr, "%s: can not read input image '%s'\n", proc, theim_name );
      return( -1 );
  }



  /* initialization of result image
   * - with slice information (extraction of slices)
   * - with template image
   * - with parameters, if any
   */



  /* XY slice extraction
   */
  if ( par.slice.z >= par.originValue && par.slice.z < (int)theim.nplanes+par.originValue ) {
      if ( BAL_AllocFullImage( &resim, (char*)NULL,
             theim.ncols, theim.nrows, 1,
             theim.vdim, 1.0, 1.0, 1.0, theim.type ) != 1 ) {
        BAL_FreeTransformation( &theTrsf );
        BAL_FreeImage( &theim );
        fprintf( stderr, "%s: unable to initialize result image (XY slice)\n", proc );
        return( -1 );
      }
      resim.vx = theim.vx;
      resim.vy = theim.vy;
      resim.vz = theim.vz;
      if ( BAL_SetImageVoxelSizes( &resim, resim.vx, resim.vy, resim.vz ) != 1 ) {
        BAL_FreeTransformation( &theTrsf );
        BAL_FreeImage( &theim );
        BAL_FreeImage( &resim );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
        return( -1 );
      }

      theTrsf.mat.m[11] = par.slice.z - par.originValue;

      theLeftCorner[2] = par.slice.z - par.originValue;
  }



  /* XZ slice extraction
   */
  else if ( par.slice.y >= par.originValue && par.slice.y < (int)theim.nrows+par.originValue ) {
    if ( BAL_AllocFullImage( &resim, (char*)NULL,
           theim.ncols, theim.nplanes, 1,
           theim.vdim, 1.0, 1.0, 1.0, theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      fprintf( stderr, "%s: unable to initialize result image (XZ slice)\n", proc );
      return( -1 );
    }
    resim.vx = theim.vx;
    resim.vy = theim.vz;
    resim.vz = theim.vy;
    if ( BAL_SetImageVoxelSizes( &resim, resim.vx, resim.vy, resim.vz ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
      return( -1 );
    }

    theTrsf.mat.m[0] = theTrsf.mat.m[5] = theTrsf.mat.m[10] = 0.0;
    theTrsf.mat.m[0] = 1.0;
    theTrsf.mat.m[6] = -1.0;
    theTrsf.mat.m[9] = 1.0;
    theTrsf.mat.m[7] = (par.analyzeFiji == 0) ? par.slice.y - par.originValue : (int)theim.nrows - par.slice.y + par.originValue - 1;

    theLeftCorner[1] = (par.analyzeFiji == 0) ? par.slice.y - par.originValue : (int)theim.nrows - par.slice.y + par.originValue - 1;
  }



  /* YZ slice extraction
   */
  else if ( par.slice.x >= par.originValue && par.slice.x < (int)theim.ncols+par.originValue ) {
    if ( BAL_AllocFullImage( &resim, (char*)NULL,
           theim.nrows, theim.nplanes, 1,
           theim.vdim, 1.0, 1.0, 1.0, theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image (XZ slice)\n", proc );
      return( -1 );
    }
    resim.vx = theim.vy;
    resim.vy = theim.vz;
    resim.vz = theim.vx;
    if ( BAL_SetImageVoxelSizes( &resim, resim.vx, resim.vy, resim.vz ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
      return( -1 );
    }

    theTrsf.mat.m[0] = theTrsf.mat.m[5] = theTrsf.mat.m[10] = 0.0;
    theTrsf.mat.m[4] = 1.0;
    theTrsf.mat.m[9] = 1.0;
    theTrsf.mat.m[2] = 1.0;
    theTrsf.mat.m[3] = par.slice.x - par.originValue;

    theLeftCorner[0] = par.slice.x - par.originValue;
  }



  /* initialization with a template image
   */
  else if ( template_image_name != NULL && template_image_name[0] != '\0' ) {

    if ( BAL_ReadImage( &templateim, template_image_name, 0 ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( _verbose_ )
        fprintf( stderr, "%s: can not read template image '%s'\n", proc, template_image_name );
      return( -1 );
    }

    theLeftCorner[0] = par.origin.x - par.originValue;
    theLeftCorner[1] = (par.analyzeFiji == 0) ? par.origin.y - par.originValue: (int)theim.nrows - (int)templateim.nrows - par.origin.y + par.originValue;
    theLeftCorner[2] = par.origin.z - par.originValue;

    theRightCorner[0] = theLeftCorner[0] + templateim.ncols;
    theRightCorner[1] = theLeftCorner[1] + templateim.nrows;
    theRightCorner[2] = theLeftCorner[2] + templateim.nplanes;

    BAL_FreeImage( &templateim );

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;
    if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;

    if ( theRightCorner[0] > (int)theim.ncols ) theRightCorner[0] = theim.ncols;
    if ( theRightCorner[1] > (int)theim.nrows ) theRightCorner[1] = theim.nrows;
    if ( theRightCorner[2] > (int)theim.nplanes ) theRightCorner[2] = theim.nplanes;

    if ( BAL_AllocFullImage( &resim, (char*)NULL,
           theRightCorner[0] - theLeftCorner[0],
           theRightCorner[1] - theLeftCorner[1],
           theRightCorner[2] - theLeftCorner[2],
           theim.vdim, 1.0, 1.0, 1.0,  theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image (from template image)\n", proc );
      return( -1 );
    }

    resim.vx = theim.vx;
    resim.vy = theim.vy;
    resim.vz = theim.vz;
    if ( BAL_SetImageVoxelSizes( &resim, resim.vx, resim.vy, resim.vz ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
      return( -1 );
    }

    theTrsf.mat.m[ 3] = theLeftCorner[0];
    theTrsf.mat.m[ 7] = theLeftCorner[1];
    theTrsf.mat.m[11] = theLeftCorner[2];
  }



  /* initialisation with dimension parameters
   */
  else if ( par.dim.x > 0 && par.dim.y > 0 ) {

    theLeftCorner[0] = par.origin.x - par.originValue;
    theLeftCorner[1] = (par.analyzeFiji == 0) ? par.origin.y - par.originValue : (int)theim.nrows - par.dim.y - par.origin.y + par.originValue;

    theRightCorner[0] = theLeftCorner[0] + par.dim.x;
    theRightCorner[1] = theLeftCorner[1] + par.dim.y;

    if ( theLeftCorner[0] < 0 ) theLeftCorner[0] = 0;
    if ( theLeftCorner[1] < 0 ) theLeftCorner[1] = 0;

    if ( theRightCorner[0] > (int)theim.ncols ) theRightCorner[0] = theim.ncols;
    if ( theRightCorner[1] > (int)theim.nrows ) theRightCorner[1] = theim.nrows;

    if ( par.dim.z > 0 ) {

      theLeftCorner[2] = par.origin.z - par.originValue;
      theRightCorner[2] = theLeftCorner[2] + par.dim.z - par.originValue;

      if ( theLeftCorner[2] < 0 ) theLeftCorner[2] = 0;
      if ( theRightCorner[2] > (int)theim.nplanes ) theRightCorner[2] = theim.nplanes;

    }
    else {

      theLeftCorner[2] = 0;
      theRightCorner[2] = 1;

    }

    if ( BAL_AllocFullImage( &resim, (char*)NULL,
           theRightCorner[0] - theLeftCorner[0],
           theRightCorner[1] - theLeftCorner[1],
           theRightCorner[2] - theLeftCorner[2],
           theim.vdim, 1.0, 1.0, 1.0, theim.type ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image (from parameters)\n", proc );
      return( -1 );
    }

    resim.vx = theim.vx;
    resim.vy = theim.vy;
    resim.vz = theim.vz;
    if ( BAL_SetImageVoxelSizes( &resim, resim.vx, resim.vy, resim.vz ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
      return( -1 );
    }

    theTrsf.mat.m[ 3] = theLeftCorner[0];
    theTrsf.mat.m[ 7] = theLeftCorner[1];
    theTrsf.mat.m[11] = theLeftCorner[2];
  }

  /* do not know what to do ?!
   */
  else {
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeImage( &theim );
    if ( _verbose_ )
      fprintf( stderr, "%s: no (valid?) geometrical information?\n", proc );
    return( -1 );
  }



  /***************************************************
   *
   * here is the API call
   *
   ***************************************************/

  if ( API_cropImage( &theim, &resim, param_str_1, param_str_2 ) != 1 ) {
    BAL_FreeImage( &resim );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeImage( &theim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to crop image\n", proc );
    return( -1 );
  }


  /***************************************************
   *
   * writing results
   *
   ***************************************************/

  if ( voxel_transformation_name != (char*)NULL && voxel_transformation_name[0] != '\0' ) {
    if ( BAL_WriteTransformation( &theTrsf, voxel_transformation_name ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write result voxel transformation '%s'\n", proc, voxel_transformation_name  );
      return( -1 );
    }
  }

  if ( real_transformation_name != (char*)NULL && real_transformation_name[0] != '\0' ) {
    if ( BAL_ChangeTransformationToRealUnit( &theim, &resim, &theTrsf, &theTrsf ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to convert 'voxel' transformation into the 'real' world\n", proc );
      return( -1 );
    }
    if ( BAL_WriteTransformation( &theTrsf, real_transformation_name ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      BAL_FreeImage( &theim );
      BAL_FreeImage( &resim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write result voxel transformation '%s'\n", proc, real_transformation_name  );
      return( -1 );
    }
  }

  BAL_FreeTransformation( &theTrsf );
  BAL_FreeImage( &theim );


  if ( BAL_WriteImage( &resim, resim_name ) != 1 ) {
    BAL_FreeImage( &resim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write result image '%s'\n", proc, resim_name );
    return( -1 );
  }

  BAL_FreeImage( &resim );

  return( 1 );
}





int API_cropImage( bal_image *image, bal_image *imres, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_cropImage";
  bal_image imtmp;
  bal_image *imptr = (bal_image*)NULL;
  int theDim[3];
  int resDim[3] = {-1, -1, -1};
  int theLeftCorner[3] = {0, 0, 0};

  bal_transformation theTrsf;
  bal_transformation invTrsf;
  _MATRIX sform_to_voxel;

  lineCmdParamCropImage par;



  /* parameter initialization
   */
  API_InitParam_cropImage( &par );

  /* tests and initialization
   */
  if ( image->vdim != imres->vdim ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different vectorial dimensions\n", proc );
    return( -1 );
  }
  if ( image->vdim != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: vectorial images not handled yet\n", proc );
    return( -1 );
  }

  BAL_InitImage( &imtmp, NULL, 0, 0, 0, 0, UCHAR );
  if ( imres->type != image->type ) {
    if ( BAL_AllocImageFromImage( &imtmp, (char*)NULL,
                                      imres, image->type ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: auxiliary image allocation failed\n", proc );
      return( -1 );
    }
    imptr = &imtmp;
  }
  else {
    imptr = imres;
  }

  BAL_InitTransformation( &theTrsf );
  if ( BAL_AllocTransformation( &theTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    if ( _verbose_ )
      fprintf( stderr, "%s: Unable to allocate voxel transformation \n", proc);
    return( -1 );
  }
  BAL_SetTransformationToIdentity( &theTrsf );
  theTrsf.transformation_unit = VOXEL_UNIT;



  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_cropImage( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_cropImage( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_cropImage( stderr, proc, &par, (char*)NULL );


  theDim[0] = image->ncols;
  theDim[1] = image->nrows;
  theDim[2] = image->nplanes;

  resDim[0] = imres->ncols;
  resDim[1] = imres->nrows;
  resDim[2] = imres->nplanes;


  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  /* XY slice extraction
   */
  if ( par.slice.z >= par.originValue && par.slice.z < (int)image->nplanes+par.originValue ) {

    if ( imres->ncols != image->ncols || imres->nrows != image->nrows || imres->nplanes != 1 ) {
      if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: no compatible dimensions with XY slice extraction\n", proc );
      return( -1 );
    }

    resDim[0] = image->ncols;
    resDim[1] = image->nrows;
    resDim[2] = 1;

    theLeftCorner[2] = par.slice.z - par.originValue;

    theTrsf.mat.m[11] = theLeftCorner[2];

  }

  /* XZ slice extraction
   */
  else if ( par.slice.y >= par.originValue && par.slice.y < (int)image->nrows+par.originValue ) {

    if ( imres->ncols != image->ncols || imres->nrows != image->nplanes || imres->nplanes != 1 ) {
      if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: no compatible dimensions with XZ slice extraction\n", proc );
      return( -1 );
    }

    resDim[0] = image->ncols;
    resDim[1] = 1;
    resDim[2] = image->nplanes;

    theLeftCorner[1] = (par.analyzeFiji == 0) ? par.slice.y - par.originValue : (int)image->nrows - par.slice.y + par.originValue - 1;

    theTrsf.mat.m[0] = theTrsf.mat.m[5] = theTrsf.mat.m[10] = 0.0;
    theTrsf.mat.m[0] = 1.0;
    theTrsf.mat.m[6] = -1.0;
    theTrsf.mat.m[9] = 1.0;
    theTrsf.mat.m[7] = theLeftCorner[1];

  }

  /* YZ slice extraction
   */
  else if ( par.slice.x >= par.originValue && par.slice.x < (int)image->ncols+par.originValue ) {

    if ( imres->ncols != image->nrows || imres->nrows != image->nplanes || imres->nplanes != 1 ) {
      if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: no compatible dimensions with YZ slice extraction\n", proc );
      return( -1 );
    }

    resDim[0] = 1;
    resDim[1] = image->nrows;
    resDim[2] = image->nplanes;

    theLeftCorner[0] = par.slice.x - par.originValue;

    theTrsf.mat.m[0] = theTrsf.mat.m[5] = theTrsf.mat.m[10] = 0.0;
    theTrsf.mat.m[4] = 1.0;
    theTrsf.mat.m[9] = 1.0;
    theTrsf.mat.m[2] = 1.0;
    theTrsf.mat.m[3] = theLeftCorner[0];

  }

  /* parameters
   */
  else {

    theLeftCorner[0] = par.origin.x - par.originValue;
    theLeftCorner[1] = (par.analyzeFiji == 0) ? par.origin.y - par.originValue: (int)image->nrows - (int)imres->nrows - par.origin.y + par.originValue;
    theLeftCorner[2] = par.origin.z - par.originValue;

    if ( theLeftCorner[0]+imres->ncols > image->ncols
         || theLeftCorner[1]+imres->nrows > image->nrows
         || theLeftCorner[2]+imres->nplanes > image->nplanes ) {
      if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: incompatible dimensions\n", proc );
      return( -1 );
    }

    theTrsf.mat.m[ 3] = theLeftCorner[0];
    theTrsf.mat.m[ 7] = theLeftCorner[1];
    theTrsf.mat.m[11] = theLeftCorner[2];

  }



  /***************************************************
   *
   *
   *
   ***************************************************/


  ExtractFromBuffer( image->data, theDim, imres->data, resDim, theLeftCorner, image->type );


  /* image geometry calculation
   */

  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation \n", proc );
    return( -1 );
  }
  theTrsf.transformation_unit = VOXEL_UNIT;

  if ( BAL_InverseTransformation( &theTrsf, &invTrsf ) != 1 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: Unable to invert transformation \n", proc );
    return( -1 );
  }

  /* re-computation of Qform matrices
   */

  _mult_mat( &(invTrsf.mat), &(image->to_voxel), &(imres->to_voxel) );

  if ( InverseMat4x4( imres->to_voxel.m, imres->to_real.m ) != 4 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert qform matrix\n", proc );
    return( -1 );
  }

  imres->geometry = _BAL_QFORM_GEOMETRY_;
  imres->qform_code = image->qform_code;
  if ( imres->qform_code == 0 )
    imres->qform_code = 1;

  /* re-computation of Sform matrices
   */

  _init_mat( &(sform_to_voxel) );
  if ( _alloc_mat( &(sform_to_voxel), 4, 4 ) != 1 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate sform 'to_voxel' matrix\n", proc );
    return( -1 );
  }

  if ( InverseMat4x4( image->sform_to_real.m, sform_to_voxel.m ) != 4 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &invTrsf );
    _free_mat( &(sform_to_voxel) );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to inverse sform real matrix\n", proc );
    return( -1 );
  }

  _mult_mat( &(invTrsf.mat), &sform_to_voxel, &sform_to_voxel );

  if ( InverseMat4x4( sform_to_voxel.m, imres->sform_to_real.m ) != 4 ) {
    if ( imres->type != image->type ) BAL_FreeImage( &imtmp );
    BAL_FreeTransformation( &theTrsf );
    BAL_FreeTransformation( &invTrsf );
    _free_mat( &(sform_to_voxel) );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to inverse sform voxel matrix\n", proc );
    return( -1 );
  }

  imres->sform_code = image->sform_code;
  if ( imres->sform_code == 0 )
    imres->sform_code = 1;

  BAL_FreeTransformation( &theTrsf );
  BAL_FreeTransformation( &invTrsf );
  _free_mat( &(sform_to_voxel) );


  /***************************************************
   *
   *
   *
   ***************************************************/


  /* copy, if required
   */
  if ( imres->type != image->type ) {
    if ( BAL_CopyImage( imptr, imres ) != 1 )  {
      BAL_FreeImage( &imtmp );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy auxiliary image\n", proc );
      return( -1 );
    }
    BAL_FreeImage( &imtmp );
  }

  return( 1 );
}






/************************************************************
 *
 * static functions
 *
 ************************************************************/



static char **_Str2Array( int *argc, char *str )
{
  char *proc = "_Str2Array";
  int n = 0;
  char *s = str;
  char **array, **a;

  if ( s == (char*)NULL || strlen( s ) == 0 ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: empty input string\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* go to the first valid character
   */
  while ( *s == ' ' || *s == '\n' || *s == '\t' )
    s++;

  if ( *s == '\0' ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: weird, input string contains only separation characters\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  /* count the number of strings
   */
  for ( n = 0; *s != '\0'; ) {
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' )
      s ++;
  }

  if ( _verbose_ >= 5 )
    fprintf( stderr, "%s: found %d strings\n", proc, n );

  /* the value of the strings will be duplicated
   * so that the input string can be freed
   */
  array = (char**)vtmalloc( n * sizeof(char*) + (strlen(str)+1) * sizeof(char),
                            "array", proc );
  if ( array == (char**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    *argc = 0;
    return( (char**)NULL );
  }

  a = array;
  a += n;
  s = (char*)a;
  (void)strncpy( s, str, strlen( str ) );
  s[ strlen( str ) ] = '\0';

  while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
    *s = '\0';
    s++;
  }

  for ( n = 0; *s != '\0'; ) {
    array[n] = s;
    n ++;
    while ( *s != ' ' && *s != '\n' && *s != '\t' && *s != '\0' )
      s ++;
    while ( *s == ' ' || *s == '\n' || *s == '\t' ) {
      *s = '\0';
      s ++;
    }
  }

  *argc = n;
  return( array );
}





/************************************************************
 *
 * help / documentation
 *
 ************************************************************/



static char *usage = "[image-in] [image-out]\n\
 [-result-transformation | -res-trsf %s]\n\
 [-result-voxel-transformation | -res-voxel-trsf %s]\n\
 [-template | -t | -dims %s]\n\
 [-origin | -o %d %d [%d] | [-ix %d] [-iy %d] [-iz %d]]\n\
 [-dim %d %d [%d] | [-x %d] [-y %d] [-z %d]]\n\
 [-xy %d | -yz %d | -xz %d]\n\
 [-0 | -1]\n\
 [-analyze-fiji | -fiji]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-nv|-noverbose] [-debug|-D] [-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-notime]\n\
 [-help|-h]";



static char *detail = "\
 if 'image-in' is equal to '-', stdin will be used\n\
 if 'image-out' is not specified or equal to '-', stdout will be used\n\
 if both are not specified, stdin and stdout will be used\n\
# transformation names\n\
-result-transformation %s # resampling transformation for 'applyTrsf'\n\
    'applyTrsf %s %s -trsf this-transformation -dim %d %d %d'\n\
    yields then the same result. \n\
    Useful to compose cropping with transformations.\n\
-result-voxel-transformation %s # idem with a transformation encoded in voxels\n\
# template\n\
-template | -t | -dims %s # template image to define dimensions of the sub-image\n\
    to be extracted\n\
# sub-image definition\n\
-origin | -o %d %d [%d] # origin of the sub-image to be extrated\n\
    depending of flag '-0' or '-1', coordinates begin at 0 (default)\n\
    or 1\n\
-ix %d # origin along X\n\
-iy %d # origin along Y\n\
-iz %d # origin along Z\n\
-dim %d %d [%d] # dimensions of the sub-image to be extrated\n\
    this is an alternative to '-template %s'\n\
-x %d  # dimension along X\n\
-y %d  # dimension along Y\n\
-z %d  # dimension along Z\n\
# slice definition\n\
-xy %d # extract the XY slice #%d\n\
       # this is equivalent to '-origin 0 0 %d -dim dimx dimy 1'\n\
-xz %d # extract the XZ slice #%d\n\
       # the output image has sizes (dimx, dimz, 1)\n\
-yz %d # extract the YZ slice #%d\n\
       # the output image has sizes (dimy, dimz, 1)\n\
# general parameters\n\
-0 | -1 # specifying the origin of coordinates (ie index of first point)\n\
-analyze-fiji | -fiji # analyze image may be mirrored wrt X axis\n\
    (eg Y axis is reverted). This converted Y coordinates of origin into\n\
    (input.image.dim.y - origin.y - 1) so that values picked in fiji can be\n\
    directly used.\n\
# parallelism parameters\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general parameters \n\
  -verbose|-v: increase verboseness\n\
    parameters being read several time, use '-nv -v -v ...'\n\
    to set the verboseness level\n\
  -noverbose|-nv: no verboseness at all\n\
  -debug|-D: increase debug level\n\
  -nodebug: no debug indication\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -h: print option list\n\
  -help: print option list + details\n\
";





char *API_Help_cropImage( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_cropImage( char *proc, char *str, int flag )
{
    if ( flag >= 0 ) {
        if ( proc != (char*)NULL )
           (void)fprintf(stderr,"Usage: %s %s\n", proc, usage);
        else
            (void)fprintf(stderr,"Command line options: %s\n", usage);
    }
    if ( flag == 1 ) {
      (void)fprintf( stderr, "--------------------------------------------------\n" );
      (void)fprintf(stderr,"%s",detail);
      (void)fprintf( stderr, "--------------------------------------------------\n" );
    }
    if ( str != (char*)NULL )
      (void)fprintf(stderr,"Error: %s\n",str);
    exit( 1 );
}





/************************************************************
 *
 * parameters management
 *
 ************************************************************/



void API_InitParam_cropImage( lineCmdParamCropImage *p )
{
    (void)strncpy( p->input_name, "\0", 1 );
    (void)strncpy( p->output_name, "\0", 1 );

    (void)strncpy( p->output_real_transformation_name, "\0", 1 );
    (void)strncpy( p->output_voxel_transformation_name, "\0", 1 );

    (void)strncpy( p->template_name, "\0", 1 );

    p->origin.x = 0;
    p->origin.y = 0;
    p->origin.z = 0;

    p->dim.x = -1;
    p->dim.y = -1;
    p->dim.z = -1;

    p->slice.x = -1;
    p->slice.y = -1;
    p->slice.z = -1;

    p->originValue = 0;
    p->analyzeFiji = 0;

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_cropImage( FILE *theFile, char *proc,
                                  lineCmdParamCropImage *p, char *str )
{
  FILE *f = theFile;
  if ( theFile == (FILE*)NULL ) f = stderr;

  fprintf( f, "==================================================\n" );
  fprintf( f, "= in line command parameters" );
  if ( proc != (char*)NULL )
    fprintf( f, " for '%s'", proc );
  if ( str != (char*)NULL )
    fprintf( f, "= %s\n", str );
  fprintf( f, "\n"  );
  fprintf( f, "==================================================\n" );


  fprintf( f, "# file names\n" );

  fprintf( f, "- input image is " );
  if ( p->input_name != (char*)NULL && p->input_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image is " );
  if ( p->output_name != (char*)NULL && p->output_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output real transformation image is " );
  if ( p->output_name != (char*)NULL && p->output_real_transformation_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_real_transformation_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output voxel transformation image is " );
  if ( p->output_name != (char*)NULL && p->output_voxel_transformation_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_voxel_transformation_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- template image is " );
  if ( p->output_name != (char*)NULL && p->template_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->template_name );
  else
    fprintf( f, "'NULL'\n" );


  fprintf( f, "# sub-image definition\n" );

  fprintf( f, "- sub-image origin is [%d %d %d]\n",
           p->origin.x, p->origin.y, p->origin.z );
  fprintf( f, "- sub-image dimension is [%d %d %d]\n", p->dim.x, p->dim.y, p->dim.z );

  fprintf( f, "- YZ slice index is %d\n", p->slice.x );
  fprintf( f, "- XZ slice index is %d\n", p->slice.y );
  fprintf( f, "- XY slice index is %d\n", p->slice.z );

  /* fprintf( f, "# general image related parameters\n" ); */

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_cropImage( char *str, lineCmdParamCropImage *p )
{
  char *proc = "_API_ParseParam_cropImage";
  char **argv;
  int i, argc;

  if ( str == (char*)NULL || strlen(str) == 0 )
      return;

  argv = _Str2Array( &argc, str );
  if ( argv == (char**)NULL || argc == 0 ) {
      if ( _debug_ ) {
          fprintf( stderr, "%s: weird, no arguments were found\n", proc );
      }
      return;
  }

  if ( _debug_ > 4 ) {
      fprintf( stderr, "%s: translation from\n", proc );
      fprintf( stderr, "   '%s'\n", str );
      fprintf( stderr, "into\n" );
      for ( i=0; i<argc; i++ )
          fprintf( stderr, "   argv[%2d] = '%s'\n", i, argv[i] );
  }

  API_ParseParam_cropImage( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_cropImage( int firstargc, int argc, char *argv[],
                                  lineCmdParamCropImage *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  float t;
  int status;
  int maxchunks;

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
            if ( inputisread == 0 ) {
              (void)strcpy( p->input_name,  "<" );  /* standart input */
              inputisread = 1;
            }
            else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name,  ">" );  /* standart output */
              outputisread = 1;
            }
            else {
              API_ErrorParse_cropImage( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
            }
          }

          /* transformation names
           */
          else if ( strcmp ( argv[i], "-result-transformation") == 0
                    || (strcmp ( argv[i], "-res-trsf") == 0 && argv[i][9] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_cropImage( (char*)NULL, "parsing -result-transformation", 0 );
            (void)strcpy( p->output_real_transformation_name, argv[i] );
          }

          else if ( strcmp ( argv[i], "-result-voxel-transformation") == 0
                    || strcmp ( argv[i], "-res-voxel-trsf") == 0 ) {
            i++;
            if ( i >= argc) API_ErrorParse_cropImage( (char*)NULL, "parsing -result-voxel-transformation", 0 );
            (void)strcpy( p->output_voxel_transformation_name, argv[i] );
          }

          /* template names
           */
          else if ( strcmp ( argv[i], "-template") == 0
                    || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_cropImage( (char*)NULL, "parsing -template", 0 );
            (void)strcpy( p->template_name, argv[i] );
          }

          /* sub-image definition
           */
          else if ( (strcmp (argv[i], "-origin" ) == 0 && argv[i][7] == '\0')
                    || (strcmp (argv[i], "-o" ) == 0 && argv[i][2] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -origin %d", 0 );
            status = sscanf( argv[i], "%f", &t );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -origin %d", 0 );
            p->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -origin %d %d", 0 );
            status = sscanf( argv[i], "%f", &t );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -origin %d %d", 0 );
            p->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
            i ++;
            if ( i >= argc) p->origin.z = 0;
            else {
              status = sscanf( argv[i], "%f", &t );
              if ( status <= 0 ) {
                i--;
                p->origin.z = 0;
              }
              else {
                p->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
              }
            }
          }

          else if ( strcmp (argv[i], "-ix" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -ix %d", 0 );
            status = sscanf( argv[i], "%f", &t );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -ix %d", 0 );
            p->origin.x = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
          }
          else if ( strcmp (argv[i], "-iy" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -iy %d", 0 );
            status = sscanf( argv[i], "%f", &t );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -iy %d", 0 );
            p->origin.y = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
          }
          else if ( strcmp (argv[i], "-iz" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -iz %d", 0 );
            status = sscanf( argv[i], "%f", &t );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -iz %d", 0 );
            p->origin.z = (t > 0.0) ? (int)(t+0.5) : (int)(t-0.5);
          }

          else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -dim %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.x) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -dim %d", 0 );
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -dim %d %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.y) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -dim %d %d", 0 );
            i ++;
            if ( i >= argc) p->dim.z = 1;
            else {
              status = sscanf( argv[i], "%d", &(p->dim.z) );
              if ( status <= 0 ) {
                i--;
                p->dim.z = 1;
              }
            }
          }

          else if ( strcmp (argv[i], "-x" ) == 0 && argv[i][2] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -x %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.x) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -x %d", 0 );
          }
          else if ( strcmp (argv[i], "-y" ) == 0 && argv[i][2] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -y %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.y) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -y %d", 0 );
          }
          else if ( strcmp (argv[i], "-z" ) == 0 && argv[i][2] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -z %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.z) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -z %d", 0 );
          }

          /* slice definition
           */
          else if ( strcmp (argv[i], "-xy" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -xy %d", 0 );
            status = sscanf( argv[i], "%d", &(p->slice.z) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -xy %d", 0 );
          }
          else if ( strcmp (argv[i], "-yz" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -yz %d", 0 );
            status = sscanf( argv[i], "%d", &(p->slice.x) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -yz %d", 0 );
          }
          else if ( strcmp (argv[i], "-xz" ) == 0 && argv[i][3] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -xz %d", 0 );
            status = sscanf( argv[i], "%d", &(p->slice.y) );
            if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -xz %d", 0 );
          }

          /* other parameters
           */
          else if ( strcmp (argv[i], "-0" ) == 0 && argv[i][2] == '\0' ) {
            p->originValue = 0;
          }
          else if ( strcmp (argv[i], "-1" ) == 0 && argv[i][2] == '\0' ) {
            p->originValue = 1;
          }

          else if ( (strcmp (argv[i], "-fiji" ) == 0 && argv[i][5] == '\0')
                    || (strcmp (argv[i], "-analyze-fiji" ) == 0) ) {
            p->analyzeFiji = 1;
          }





          /* parallelism parameters
           */
          else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
             setParallelism( _DEFAULT_PARALLELISM_ );
          }

          else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
             setParallelism( _NO_PARALLELISM_ );
          }

          else if ( strcmp ( argv[i], "-parallelism-type" ) == 0 ||
                      strcmp ( argv[i], "-parallel-type" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             if ( strcmp ( argv[i], "default" ) == 0 ) {
               setParallelism( _DEFAULT_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "none" ) == 0 ) {
               setParallelism( _NO_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "openmp" ) == 0 || strcmp ( argv[i], "omp" ) == 0 ) {
               setParallelism( _OMP_PARALLELISM_ );
             }
             else if ( strcmp ( argv[i], "pthread" ) == 0 || strcmp ( argv[i], "thread" ) == 0 ) {
               setParallelism( _PTHREAD_PARALLELISM_ );
             }
             else {
               fprintf( stderr, "unknown parallelism type: '%s'\n", argv[i] );
               API_ErrorParse_cropImage( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_cropImage( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_cropImage( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
             if ( strcmp ( argv[i], "default" ) == 0 ) {
               setOmpScheduling( _DEFAULT_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "static" ) == 0 ) {
               setOmpScheduling( _STATIC_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
               setOmpScheduling( _DYNAMIC_ONE_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
               setOmpScheduling( _DYNAMIC_OMP_SCHEDULING_ );
             }
             else if ( strcmp ( argv[i], "guided" ) == 0 ) {
               setOmpScheduling( _GUIDED_OMP_SCHEDULING_ );
             }
             else {
               fprintf( stderr, "unknown omp scheduling type: '%s'\n", argv[i] );
               API_ErrorParse_cropImage( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
             }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_cropImage( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_cropImage( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
            }
          }
          else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                    || strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _debug_ <= 0 ) _debug_ = 1;
              else                _debug_ ++;
            }
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
          }

          else if ( strcmp ( argv[i], "-print-parameters" ) == 0
                    || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
             p->print_lineCmdParam = 1;
          }

          else if ( strcmp ( argv[i], "-print-time" ) == 0
                     || (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
             p->print_time = 1;
          }
          else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
             p->print_time = 0;
          }

          /* unknown option
           */
          else {
              sprintf(text,"unknown option %s\n",argv[i]);
              API_ErrorParse_cropImage( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_cropImage( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( inputisread == 0 ) {
              (void)strcpy( p->input_name, argv[i] );
              inputisread = 1;
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_cropImage( (char*)NULL, "too many file names ...\n", 0 );
          }
      }
  }

  /* if not enough file names
   */
  if ( inputisread == 0 ) {
    (void)strcpy( p->input_name,  "<" );  /* standart input */
    inputisread = 1;
  }
  if ( outputisread == 0 ) {
    (void)strcpy( p->output_name,  ">" );  /* standart output */
    outputisread = 1;
  }

}
