/*************************************************************************
 * api-copyTrsf.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 24 jui 2015 17:33:43 CEST
 *
 * ADDITIONS, CHANGES
 *
 * to generate files:
 * sed -e "s/copyTrsf/execuTable/g" \
 *     -e "s/CopyTrsf/ExecuTable/g" \
 *     -e "s/copytrsf/executable/g" \
 *     [api-]copyTrsf.[c,h] > [api-]execTable.[c,h]
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <bal-transformation.h>
#include <bal-transformation-copy.h>

#include <api-copyTrsf.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_copyTrsf( char *str, lineCmdParamCopyTrsf *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/



/* this one is kept for historical reasons,
 * ie Stracking compilation but should disappear
 */
int copyTrsf( char* thetrsf_name,
              char* restrsf_name,
              enumUnitTransfo thetrsf_unit,
              enumUnitTransfo restrsf_unit,
              char *template_image_name,
              char *floating_image_name,
              bal_integerPoint dim,
              bal_doublePoint voxel,
              enumTypeTransfo transformation_type,
              int isDebug,
              int isVerbose )
{
  char *proc = "copyTrsf";
  char str[256], *s;

  _verbose_ = isVerbose;
  _debug_ = isDebug;

  if ( _verbose_ )
      fprintf( stderr, "Warning, '%s' is obsolete\n", proc );

  s = str;

  switch( transformation_type ) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: unknown transformation type\n", proc );
      return( -1 );
  case TRANSLATION_2D :
      sprintf( s, "-transformation-type translation2D " );
      break;
  case TRANSLATION_3D :
      sprintf( s, "-transformation-type translation3D " );
      break;
  case TRANSLATION_SCALING_2D :
      sprintf( s, "-transformation-type translation-scaling2D " );
      break;
  case TRANSLATION_SCALING_3D :
      sprintf( s, "-transformation-type translation-scaling3D " );
      break;
  case RIGID_2D :
      sprintf( s, "-transformation-type rigid2D " );
      break;
  case RIGID_3D :
      sprintf( s, "-transformation-type rigid3D " );
      break;
  case SIMILITUDE_2D :
      sprintf( s, "-transformation-type similitude2D " );
      break;
  case SIMILITUDE_3D :
      sprintf( s, "-transformation-type similitude3D " );
      break;
  case AFFINE_2D :
      sprintf( s, "-transformation-type affine2D " );
      break;
  case AFFINE_3D :
      sprintf( s, "-transformation-type affine3D " );
      break;
  case VECTORFIELD_2D :
      sprintf( s, "-transformation-type vectorfield2D " );
      break;
  case VECTORFIELD_3D :
      sprintf( s, "-transformation-type vectorfield3D " );
      break;
  }
  s = &(str[strlen(str)]);

  if ( dim.x > 0 && dim.y > 0 ) {
      sprintf( s, "-dim %d %d ", dim.x, dim.y );
      s = &(str[strlen(str)]);
      if ( dim.z > 0 ) {
          sprintf( s, "%d ", dim.z );
          s = &(str[strlen(str)]);
      }
  }

  if ( voxel.x > 0 && voxel.y > 0 ) {
      sprintf( s, "-voxel %f %f ", voxel.x, voxel.y );
      s = &(str[strlen(str)]);
      if ( voxel.z > 0 ) {
          sprintf( s, "%f ", voxel.z );
          s = &(str[strlen(str)]);
      }
  }

  switch( thetrsf_unit ) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: unknown input transformation unit\n", proc );
      return( -1 );
  case UNDEF_UNIT :
  case VOXEL_UNIT :
      sprintf( s, "-input-unit voxel " );
      s = &(str[strlen(str)]);
      break;
  case REAL_UNIT :
      sprintf( s, "-input-unit real " );
      s = &(str[strlen(str)]);
      break;
  }

  switch( restrsf_unit ) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: unknown input transformation unit\n", proc );
      return( -1 );
  case UNDEF_UNIT :
  case VOXEL_UNIT :
      sprintf( s, "-output-unit voxel " );
      s = &(str[strlen(str)]);
      break;
  case REAL_UNIT :
      sprintf( s, "-output-unit real " );
      s = &(str[strlen(str)]);
      break;
  }

  /* ... */
  if ( API_INTERMEDIARY_copyTrsf( thetrsf_name,
                                    restrsf_name,
                                    template_image_name,
                                    floating_image_name,
                                    str, (char*)NULL ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: some error occurs\n", proc );
      return( -1 );
  }

  return( 0 );
}





int API_INTERMEDIARY_copyTrsf( char *thetrsf_name,
                               char *restrsf_name,
                               char *template_image_name,
                               char *floating_image_name,
                               char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_copyTrsf";
  bal_transformation theTransformation;
  bal_transformation resTransformation;
  bal_transformation *resTrsf = (bal_transformation *)NULL;
  bal_image imRead;
  bal_image imFloating;
  bal_image imReference;
  bal_image *ptrFloating = (bal_image*)NULL;
  bal_image *ptrReference = (bal_image*)NULL;

  lineCmdParamCopyTrsf par;
  int verbose_balimage;



  /* parameter initialization
   */
  API_InitParam_copyTrsf( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_copyTrsf( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_copyTrsf( param_str_2, &par );



  /***************************************************
   *
   * initialization
   *
   ***************************************************/

  BAL_InitTransformation( &theTransformation );
  BAL_InitTransformation( &resTransformation );
  BAL_InitImage( &imFloating, NULL, 0, 0, 0, 0, UCHAR );
  BAL_InitImage( &imReference, NULL, 0, 0, 0, 0, UCHAR );





  /***************************************************
   *
   * reading input transformation
   *
   ***************************************************/

  if ( thetrsf_name != (char*)NULL && thetrsf_name[0] != '\0' ) {
      if ( BAL_ReadTransformation( &theTransformation, thetrsf_name ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to read '%s'\n", proc, thetrsf_name );
        return( -1 );
      }
      theTransformation.transformation_unit = par.thetrsf_unit;
  }
  else {
      if ( _verbose_ )
        fprintf( stderr, "%s: no input transformation\n", proc );
      return( -1 );
  }

  if ( _debug_ )
    BAL_PrintTransformation( stderr, &theTransformation, "read transformation" );





  /***************************************************
   *
   * reading template (if required)
   * - the reference/template is required to define the geometry
   *   of a vector field and/or to convert from voxel to real
   *   (or the inverse) for linear transformation
   * Rq: reference is embedded in vectorfield transformation
   * - the floating is required to convert from voxel to real
   *   (or the inverse)
   *
   ***************************************************/

  if ( BAL_IsTransformationTypeLinear( theTransformation.type )
       && ( BAL_IsTransformationTypeVectorField( par.transformation_type )
            || ( par.thetrsf_unit != par.restrsf_unit ) ) ) {

    /* initializing reference image, if any
     */
    if ( template_image_name != (char*)NULL && template_image_name[0] != '\0' ) {
      verbose_balimage = BAL_GetVerboseInBalImage();
      BAL_SetVerboseInBalImage( 0 );
      if ( BAL_ReadImage( &imRead, template_image_name, 1 ) != 1 ) {
        if ( _verbose_ )
            fprintf( stderr, "%s: unable to read template/reference '%s'\n", proc, template_image_name );
        BAL_SetVerboseInBalImage( verbose_balimage );
        BAL_FreeTransformation( &theTransformation );
        return( -1 );
      }
      BAL_SetVerboseInBalImage( verbose_balimage );
      if ( BAL_InitImageFromImage( &imReference, (char*)NULL, &imRead, UCHAR ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: unable to initialize reference template\n", proc );
          BAL_FreeImage( &imRead );
          BAL_FreeTransformation( &theTransformation );
          return( -1 );
      }
      BAL_FreeImage( &imRead );
      if ( par.voxel.x > 0.0 ) imReference.vx = par.voxel.x;
      if ( par.voxel.y > 0.0 ) imReference.vy = par.voxel.y;
      if ( par.voxel.z > 0.0 ) imReference.vz = par.voxel.z;
      if ( BAL_SetImageVoxelSizes( &imReference, imReference.vx, imReference.vy, imReference.vz ) != 1 ) {
          BAL_FreeImage( &imReference );
          BAL_FreeTransformation( &theTransformation );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to initialize template/reference image geometry\n", proc );
          return( -1 );
      }
      ptrReference = &imReference;
    }

    /* initializing result image
       - with parameters, if any
    */
    else if ( par.dim.x > 0 && par.dim.y > 0 ) {
      if ( par.dim.z > 0 ) {
        if ( BAL_InitImage( &imReference, (char*)NULL, par.dim.x, par.dim.y, par.dim.z, 1, UCHAR ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to initialize template/reference image\n", proc );
          BAL_FreeTransformation( &theTransformation );
          return( -1 );
        }
      }
      else {
        if ( BAL_InitImage( &imReference, (char*)NULL, par.dim.x, par.dim.y, 1, 1, UCHAR ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to initialize template/reference image (dimz=1) \n", proc );
          BAL_FreeTransformation( &theTransformation );
          return( -1 );
        }
      }
      if ( par.voxel.x > 0.0 ) imReference.vx = par.voxel.x;
      if ( par.voxel.y > 0.0 ) imReference.vy = par.voxel.y;
      if ( par.voxel.z > 0.0 ) imReference.vz = par.voxel.z;
      if ( BAL_SetImageVoxelSizes( &imReference, imReference.vx, imReference.vy, imReference.vz ) != 1 ) {
          BAL_FreeImage( &imReference );
          BAL_FreeTransformation( &theTransformation );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to initialize template/reference image geometry\n", proc );
          return( -1 );
      }
      ptrReference = &imReference;
    }
    else {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize template/reference image\n", proc );
      BAL_FreeTransformation( &theTransformation );
      return( -1 );
    }
  }
  else if ( BAL_IsTransformationTypeVectorField( theTransformation.type ) ) {
    ptrReference = &(theTransformation.vx);
  }





  /***************************************************
   *
   * unit conversion : real / voxel or converse
   * requires a floating template
   *
   ***************************************************/

  if ( par.thetrsf_unit != par.restrsf_unit ) {

    /* initializing floating image, if any
     */
    if ( floating_image_name != (char*)NULL && floating_image_name[0] != '\0' ) {
      verbose_balimage = BAL_GetVerboseInBalImage();
      BAL_SetVerboseInBalImage( 0 );
      if ( BAL_ReadImage( &imRead, floating_image_name, 1 ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: unable to read floating '%s'\n", proc, floating_image_name );
          BAL_SetVerboseInBalImage( verbose_balimage );
          BAL_FreeTransformation( &theTransformation );
          return( -1 );
      }
      BAL_SetVerboseInBalImage( verbose_balimage );
      if ( BAL_InitImageFromImage( &imFloating, (char*)NULL, &imRead, UCHAR ) != 1 ) {
          if ( _verbose_ )
              fprintf( stderr, "%s: unable to initialize floating template\n", proc );
          BAL_FreeImage( &imRead );
          BAL_FreeTransformation( &theTransformation );
          return( -1 );
      }
      BAL_FreeImage( &imRead );
      ptrFloating = &imFloating;
    }
    else {

      /* choice:
       * 1. if no floating image, use the reference one
       * 2. force to give a floating image
      */
      if ( _verbose_ ) {
          fprintf( stderr, "%s: WARNING, use the template/reference image as floating image for unit conversion\n", proc );
      }
      ptrFloating = ptrReference;
    }

  }





  /***************************************************
   *
   * copy transformation
   *
   ***************************************************/



  if ( par.transformation_type == UNDEF_TRANSFORMATION)
    par.transformation_type = theTransformation.type;


  switch ( theTransformation.type ) {
  default :
    BAL_FreeTransformation( &theTransformation );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unknow type'", proc );
      BAL_PrintTransformationType( stderr, theTransformation.type );
      fprintf( stderr, "' for transformation '%s'\n", thetrsf_name );
    }
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    switch( par.transformation_type ) {
    default :
      BAL_FreeTransformation( &theTransformation );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unknow type'", proc );
        BAL_PrintTransformationType( stderr, par.transformation_type );
        fprintf( stderr, "' for result transformation\n" );
      }
      return( -1 );

    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      resTrsf = &theTransformation;
      break;

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      if ( BAL_AllocTransformation( &resTransformation, par.transformation_type, ptrReference ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate result transformation\n", proc );
        BAL_FreeTransformation( &theTransformation );
        return( -1 );
      }
      resTransformation.transformation_unit = theTransformation.transformation_unit;

      if ( BAL_CopyTransformation( &theTransformation, &resTransformation ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to copy transformation\n", proc );
        BAL_FreeTransformation( &theTransformation );
        return( -1 );
      }
      if ( _debug_ ) BAL_PrintTransformation( stderr, &resTransformation, "after copy" );
      resTrsf = &resTransformation;
      break;

    }
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    switch( par.transformation_type ) {
    default :
      BAL_FreeTransformation( &theTransformation );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unknow type'", proc );
        BAL_PrintTransformationType( stderr, par.transformation_type );
        fprintf( stderr, "' for result transformation\n" );
      }
      return( -1 );

    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      BAL_FreeTransformation( &theTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to convert non-linear transformation into linear\n", proc );
      return( -1 );

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      resTrsf = &theTransformation;
      break;
    }

  }





  /***************************************************
   *
   * conversion
   *
   ***************************************************/

  if ( resTrsf->transformation_unit != par.thetrsf_unit ) {
    BAL_FreeTransformation( &resTransformation );
    BAL_FreeTransformation( &theTransformation );
    if ( _verbose_ )
      fprintf( stderr, "%s: weird unwanted conversion\n", proc );
    return( -1 );
  }


  switch( resTrsf->transformation_unit ) {
  default :
    BAL_FreeTransformation( &resTransformation );
    BAL_FreeTransformation( &theTransformation );
    if ( _verbose_ )
      fprintf( stderr, "%s: unknown unit for input transformation\n", proc );
    return( -1 );

  case VOXEL_UNIT :

    switch( par.restrsf_unit ) {
    default :
      BAL_FreeTransformation( &resTransformation );
      BAL_FreeTransformation( &theTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown unit for output transformation\n", proc );
      return( -1 );

    case VOXEL_UNIT :
      break;

    case REAL_UNIT :
      if ( _debug_ )
        fprintf( stderr, "%s: converting from VOXEL to REAL\n", proc );
      if ( BAL_ChangeTransformationToRealUnit( ptrFloating, ptrReference, resTrsf, resTrsf ) != 1 ) {
        BAL_FreeTransformation( &resTransformation );
        BAL_FreeTransformation( &theTransformation );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert from VOXEL to REAL\n", proc );
        return( -1 );
      }
      break;
    }

    break;

  case REAL_UNIT :

    switch( par.restrsf_unit ) {
    default :
      BAL_FreeTransformation( &resTransformation );
      BAL_FreeTransformation( &theTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown unit for output transformation\n", proc );
      return( -1 );

    case VOXEL_UNIT :
      if ( _debug_ )
        fprintf( stderr, "%s: converting from REAL to VOXEL\n", proc );
      if ( BAL_ChangeTransformationToVoxelUnit( ptrFloating, ptrReference, resTrsf, resTrsf ) != 1 ) {
        BAL_FreeTransformation( &resTransformation );
        BAL_FreeTransformation( &theTransformation );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to convert from REAL to VOXEL\n", proc );
        return( -1 );
      }
      break;
    case REAL_UNIT :
      break;
    }

    break;

  }

  if ( _debug_ ) BAL_PrintTransformation( stderr, resTrsf, "after conversion" );



  /***************************************************
   *
   * writing transformation
   *
   ***************************************************/

  if ( restrsf_name != NULL ) {
    if ( BAL_WriteTransformation( resTrsf, restrsf_name ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to write '%s'\n", proc, restrsf_name );
      BAL_FreeTransformation( &resTransformation );
      BAL_FreeTransformation( &theTransformation );
      return( -1 );
    }
  }

  BAL_FreeTransformation( &resTransformation );
  BAL_FreeTransformation( &theTransformation );


  return( 1 );
}





/* this function should be (re)write to be called from
 * API_INTERMEDIARY_copyTrsf()
 * To be done
 */
int API_copyTrsf( bal_image *image __attribute__ ((unused)),
                  bal_image *imres __attribute__ ((unused)),
                  char *param_str_1,
                  char *param_str_2 )
{
  char *proc = "API_copyTrsf";
  lineCmdParamCopyTrsf par;



  /* parameter initialization
   */
  API_InitParam_copyTrsf( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_copyTrsf( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_copyTrsf( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_copyTrsf( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  /* ... */

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


/*---------------- longueur maximale de ligne --------------------------------*/
/*----------------------------------------------------------------------------*/
static char *usage = "[transformation-in] [transformation-out]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-template|-t|-dims|-reference|-ref %s]\n\
 [-dim %d %d [%d]]\n\
 [-voxel | -pixel | -vs %f %f [%f]]\n\
 [-input-unit | -iu %s] [-output-unit | -ou %s]\n\
 [-floating|-flo %s]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [-verbose|-v] [-nv|-noverbose] [-debug|-D] [-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-notime]\n\
 [-help|-h]";


/*---------------- longueur maximale de ligne --------------------------------*/
/*----------------------------------------------------------------------------*/
static char *detail = "\
 if 'transformation-in' is equal to '-', stdin will be used\n\
 if 'transformation-out' is not specified or equal to '-', stdout will be used\n\
 if both are not specified, stdin and stdout will be used\n\
### transformation type ###\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
### template image ###\n\
  When copying/converting a matrix into a vectorfield, the vectorial image\n\
  defining the vectorfield will be created with the geometry of the template\n\
  image.\n\
[-template|-t|-dims|-reference|-ref %s] # template image for the dimensions\n\
                         of the output image\n\
[-dim %d %d [%d]]      # output image dimensions\n\
[-voxel|-pixel|-vs %f %f [%f]]    # output image voxel sizes\n\
### to change the unit of the transformation\n\
  requires to known voxel sizes of both the floating/input\n\
  and reference/template. Recall that the transformation goes\n\
  from reference to floating.\n\
[-input-unit  [voxel|real]] #\n\
[-output-unit [voxel|real]] #\n\
### ...\n\
[-floating|-flo %s]         # template image for conversion between real and voxel units\n\
  The conversion from real to voxel units is done by calculating\n\
  H^{-1}_floating o T_input o H_template\n\
  while the conversion from voxel to real units is done by calculating\n\
  H_floating o T_input o H^{-1}_template\n\
  H_{image} being the diagonal matrix of the image voxel sizes\n\
# ...\n\
# parallelism parameters\n\
 -parallel|-no-parallel:\n\
 -max-chunks %d:\n\
 -parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread:\n\
 -omp-scheduling|-omps default|static|dynamic-one|dynamic|guided:\n\
# general image related parameters\n\
   output-image-type: -o 1    : unsigned char\n\
                      -o 2    : unsigned short int\n\
                      -o 2 -s : short int\n\
                      -o 4 -s : int\n\
                      -r      : float\n\
  -type s8|u8|s16|u16|... \n\
   default is type of input image\n\
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





char *API_Help_copyTrsf( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_copyTrsf( char *program, char *str, int flag )
{
    if ( flag >= 0 ) {
        if ( program != (char*)NULL )
           (void)fprintf(stderr,"Usage: %s %s\n", program, usage);
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



void API_InitParam_copyTrsf( lineCmdParamCopyTrsf *p )
{
    (void)strncpy( p->thetrsf_name, "\0", 1 );
    (void)strncpy( p->restrsf_name, "\0", 1 );

    p->thetrsf_unit = REAL_UNIT;
    p->restrsf_unit = REAL_UNIT;

    (void)strncpy( p->template_image_name, "\0", 1 );
    (void)strncpy( p->floating_image_name, "\0", 1 );

    p->dim.x = 256;
    p->dim.y = 256;
    p->dim.z = 256;

    p->voxel.x = 1.0;
    p->voxel.y = 1.0;
    p->voxel.z = 1.0;

    p->transformation_type = UNDEF_TRANSFORMATION;

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_copyTrsf( FILE *theFile, char *program,
                                  lineCmdParamCopyTrsf *p, char *str )
{
  FILE *f = theFile;
  if ( theFile == (FILE*)NULL ) f = stderr;

  fprintf( f, "==================================================\n" );
  fprintf( f, "= in line command parameters" );
  if ( program != (char*)NULL )
    fprintf( f, " for '%s'", program );
  if ( str != (char*)NULL )
    fprintf( f, "= %s\n", str );
  fprintf( f, "\n"  );
  fprintf( f, "==================================================\n" );


  fprintf( f, "# transformation names\n" );

  fprintf( f, "- p->thetrsf_name = " );
  if ( p->thetrsf_name != (char*)NULL && p->thetrsf_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->thetrsf_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->restrsf_name = " );
  if ( p->restrsf_name != (char*)NULL && p->restrsf_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->restrsf_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# ...\n" );

  BAL_PrintTypeTransformation( f, p->transformation_type, "p->transformation_type = " );

  fprintf( f, "- p->template_image_name = " );
  if ( p->template_image_name != (char*)NULL && p->template_image_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->template_image_name );
  else
    fprintf( f, "'NULL'\n" );

  BAL_PrintIntegerPoint( f, &(p->dim), "- p->dim = " );
  BAL_PrintDoublePoint( f, &(p->voxel), "- p->voxel = " );

  fprintf( f, "- p->thetrsf_unit = " );
  switch( p->thetrsf_unit ) {
  case UNDEF_UNIT : fprintf( f, "UNDEF_UNIT\n" ); break;
  case VOXEL_UNIT : fprintf( f, "VOXEL_UNIT\n" ); break;
  case REAL_UNIT :  fprintf( f, "REAL_UNIT\n" ); break;
  }

  fprintf( f, "- p->restrsf_unit = " );
  switch( p->restrsf_unit ) {
  case UNDEF_UNIT : fprintf( f, "UNDEF_UNIT\n" ); break;
  case VOXEL_UNIT : fprintf( f, "VOXEL_UNIT\n" ); break;
  case REAL_UNIT :  fprintf( f, "REAL_UNIT\n" ); break;
  }

  fprintf( f, "- p->floating_image_name = " );
  if ( p->floating_image_name != (char*)NULL && p->floating_image_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->floating_image_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_copyTrsf( char *str, lineCmdParamCopyTrsf *p )
{
  char *proc = "_API_ParseParam_copyTrsf";
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

  API_ParseParam_copyTrsf( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_copyTrsf( int firstargc, int argc, char *argv[],
                                  lineCmdParamCopyTrsf *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
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
            (void)strcpy( p->thetrsf_name,  "<" );  /* standart input */
            inputisread = 1;
        }
        else if ( outputisread == 0 ) {
          (void)strcpy( p->restrsf_name,  ">" );  /* standart output */
          outputisread = 1;
        }
        else {
          API_ErrorParse_copyTrsf( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
        }
      }

      else if ( strcmp ( argv[i], "-transformation-type" ) == 0
           || strcmp ( argv[i], "-transformation" ) == 0
           || strcmp ( argv[i], "-trsf-type" ) == 0 ) {
         i ++;
         if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "-transformation-type", 0 );
         if ( strcmp ( argv[i], "translation2D" ) == 0 ) {
           p->transformation_type = TRANSLATION_2D;
         }
         else if ( strcmp ( argv[i], "translation3D" ) == 0 ) {
           p->transformation_type = TRANSLATION_3D;
         }
         else if ( strcmp ( argv[i], "translation" ) == 0 && argv[i][11] == '\0') {
           p->transformation_type = TRANSLATION_3D;
         }
         else if ( strcmp ( argv[i], "translation-scaling2D" ) == 0 ) {
           p->transformation_type = TRANSLATION_SCALING_2D;
         }
         else if ( strcmp ( argv[i], "translation-scaling3D" ) == 0 ) {
           p->transformation_type = TRANSLATION_SCALING_3D;
         }
         else if ( strcmp ( argv[i], "rigid2D" ) == 0 ) {
           p->transformation_type = RIGID_2D;
         }
         else if ( strcmp ( argv[i], "rigid3D" ) == 0 ) {
           p->transformation_type = RIGID_3D;
         }
         else if ( (strcmp ( argv[i], "rigid" ) == 0 && argv[i][5] == '\0') ) {
           p->transformation_type = RIGID_3D;
         }
         else if ( strcmp ( argv[i], "similitude2D" ) == 0 ) {
           p->transformation_type = SIMILITUDE_2D;
         }
         else if ( strcmp ( argv[i], "similitude3D" ) == 0 ) {
           p->transformation_type = SIMILITUDE_3D;
         }
         else if ( strcmp ( argv[i], "similitude" ) == 0 ) {
           p->transformation_type = SIMILITUDE_3D;
         }
         else if ( strcmp ( argv[i], "affine2D" ) == 0 ) {
           p->transformation_type = AFFINE_2D;
         }
         else if ( strcmp ( argv[i], "affine3D" ) == 0 ) {
           p->transformation_type = AFFINE_3D;
         }
         else if ( strcmp ( argv[i], "affine" ) == 0 ) {
           p->transformation_type = AFFINE_3D;
         }
         /*
           else if ( strcmp ( argv[i], "spline" ) == 0 ) {
           p->transformation_type = SPLINE;
           }
         */
         else if ( strcmp ( argv[i], "vectorfield" ) == 0
             || (strcmp ( argv[i], "vector" ) == 0 && argv[i][6] == '\0') ) {
           p->transformation_type = VECTORFIELD_3D;
         }
         else if ( strcmp ( argv[i], "vectorfield3D" ) == 0
             || strcmp ( argv[i], "vector3D" ) == 0 ) {
           p->transformation_type = VECTORFIELD_3D;
         }
         else if ( strcmp ( argv[i], "vectorfield2D" ) == 0
             || strcmp ( argv[i], "vector2D" ) == 0 ) {
           p->transformation_type = VECTORFIELD_2D;
         }
         else {
           fprintf( stderr, "unknown transformation type: '%s'\n", argv[i] );
           API_ErrorParse_copyTrsf( (char*)NULL, "-transformation-type", 0 );
         }
      }

      else if ( strcmp ( argv[i], "-template") == 0
                || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
                || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0')
                ||  strcmp ( argv[i], "-reference") == 0
                || (strcmp ( argv[i], "-ref") == 0 && argv[i][4] == '\0') ) {
        i++;
        if ( i >= argc) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -template", 0 );
        (void)strcpy( p->template_image_name, argv[i] );
      }

      else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -dim %d", 0 );
        status = sscanf( argv[i], "%d", &(p->dim.x) );
        if ( status <= 0 ) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -dim %d", 0 );
        i ++;
        if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -dim %d %d", 0 );
        status = sscanf( argv[i], "%d", &(p->dim.y) );
        if ( status <= 0 ) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -dim %d %d", 0 );
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

      else if ( strcmp (argv[i], "-voxel" ) == 0
                || strcmp ( argv[i], "-pixel" ) == 0
                || (strcmp ( argv[i], "-vs" ) == 0  && argv[i][3] == '\0') ) {
        i ++;
        if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -voxel %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->voxel.x) );
        if ( status <= 0 ) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -voxel %lf", 0 );
        i ++;
        if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -voxel %lf %lf", 0 );
        status = sscanf( argv[i], "%lf", &(p->voxel.y) );
        if ( status <= 0 ) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -voxel %lf %lf", 0 );
        i ++;
        if ( i >= argc) p->voxel.z = 1;
        else {
          status = sscanf( argv[i], "%lf", &(p->voxel.z) );
          if ( status <= 0 ) {
            i--;
            p->voxel.z = 1;
          }
        }
      }

      else if ( strcmp ( argv[i], "-input-unit" ) == 0
                || ( strcmp ( argv[i], "-iu" ) == 0 && argv[i][3] == '\0' ) ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "-input-unit", 0 );
          if ( strcmp( argv[i], "real" ) == 0 ) {
            p->thetrsf_unit = REAL_UNIT;
          }
          else if ( strcmp( argv[i], "voxel" ) == 0 ) {
            p->thetrsf_unit = VOXEL_UNIT;
          }
          else {
            fprintf( stderr, "unknown unit: '%s'\n", argv[i] );
            API_ErrorParse_copyTrsf( (char*)NULL, "-input-unit", 0 );
          }
      }

      else if ( strcmp ( argv[i], "-output-unit" ) == 0
                || ( strcmp ( argv[i], "-ou" ) == 0 && argv[i][3] == '\0' ) ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "-output-unit", 0 );
          if ( strcmp( argv[i], "real" ) == 0 ) {
            p->restrsf_unit = REAL_UNIT;
          }
          else if ( strcmp( argv[i], "voxel" ) == 0 ) {
            p->restrsf_unit = VOXEL_UNIT;
          }
          else {
            fprintf( stderr, "unknown unit: '%s'\n", argv[i] );
            API_ErrorParse_copyTrsf( (char*)NULL, "-output-unit", 0 );
          }
      }

      else if ( strcmp ( argv[i], "-floating") == 0
                || (strcmp ( argv[i], "-flo") == 0 && argv[i][4] == '\0') ) {
          i++;
          if ( i >= argc) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -floating", 0 );
          (void)strcpy( p->floating_image_name, argv[i] );
      }

      /* ...
       */

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
         if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
           API_ErrorParse_copyTrsf( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
         }
      }

      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
         i ++;
         if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -max-chunks ...\n", 0 );
         status = sscanf( argv[i], "%d", &maxchunks );
         if ( status <= 0 ) API_ErrorParse_copyTrsf( (char*)NULL, "parsing -max-chunks ...\n", 0 );
         if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
      }

      else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
               ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
         i ++;
         if ( i >= argc)    API_ErrorParse_copyTrsf( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
           API_ErrorParse_copyTrsf( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
         }
      }

      /* general image related parameters
       */

      /* general parameters
       */
      else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
         API_ErrorParse_copyTrsf( (char*)NULL, (char*)NULL, 1);
      }
      else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
         API_ErrorParse_copyTrsf( (char*)NULL, (char*)NULL, 0);
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
          API_ErrorParse_copyTrsf( (char*)NULL, text, 0);
      }
    }

    /* strings beginning with a character different from '-'
     */
    else {
      if ( strlen( argv[i] ) >= STRINGLENGTH ) {
          fprintf( stderr, "... parsing '%s'\n", argv[i] );
          API_ErrorParse_copyTrsf( (char*)NULL, "too long file name ...\n", 0 );
      }
      else if ( inputisread == 0 ) {
          (void)strcpy( p->thetrsf_name, argv[i] );
          inputisread = 1;
      }
      else if ( outputisread == 0 ) {
          (void)strcpy( p->restrsf_name, argv[i] );
          outputisread = 1;
      }
      else {
          fprintf( stderr, "... parsing '%s'\n", argv[i] );
          API_ErrorParse_copyTrsf( (char*)NULL, "too many file names ...\n", 0 );
      }
    }

  }

  /* if not enough file names
   */
  if ( inputisread == 0 ) {
    (void)strcpy( p->thetrsf_name,  "<" );  /* standart input */
    inputisread = 1;
  }
  if ( outputisread == 0 ) {
    (void)strcpy( p->restrsf_name,  ">" );  /* standart output */
    outputisread = 1;
  }

}
