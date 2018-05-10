/*************************************************************************
 * interpolateImages.c - template for executable creation
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  2 jul 2015 13:53:42 CEST
 *
 * ADDITIONS, CHANGES
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <vtmalloc.h>

#include <bal-transformation-copy.h>
#include <bal-interpolation.h>

#include <api-interpolateImages.h>







static int _verbose_ = 1;


/* static function definitions
 */

static char * _outputname( char *tmp, char *format, int i, float t );



/* static function definitions
 */

#ifdef UNUSED
static char *_Array2Str( int argc, char *argv[] );
#endif
static char *_BaseName( char *p );
static double _GetTime();
static double _GetClock();






int main( int argc, char *argv[] )
{
  lineCmdParamInterpolateImages par;

  bal_transformation theTrsf;

  bal_image *ptrTemplate = (bal_image*)NULL;
  bal_image *ptrTmp = (bal_image*)NULL;
  bal_image theTemplate;

  bal_image *ptrImage0 = (bal_image*)NULL;
  bal_image theImage0;
  bal_image *ptrImage1 = (bal_image*)NULL;
  bal_image theImage1;

  bal_image *ptrResImage0 = (bal_image*)NULL;
  bal_image *ptrResImage1 = (bal_image*)NULL;
  bal_image *ptrResImage = (bal_image*)NULL;
  bal_image resImage0, resImage1, resImage;

  typeCellCorrespondence theCellCorrespondence;

  int i, imin, imax;
  float t;

  char tmpname[STRINGLENGTH];

  double time_init = _GetTime();
  double time_inter = -1.0;
  double time_inter2 = -1.0;
  double time_exit;
  double clock_init = _GetClock();
  double clock_inter = -1.0;
  double clock_inter2 = -1.0;
  double clock_exit;


  /* parameter initialization
   */
  API_InitParam_interpolateImages( &par );



  /* parameter parsing
   */
  if ( argc <= 1 )
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), (char*)NULL, 0 );
  API_ParseParam_interpolateImages( 1, argc, argv, &par );
  
  if ( par.print_lineCmdParam )
      API_PrintParam_interpolateImages( stderr, _BaseName( argv[0] ), &par, (char*)NULL );



  /* reading transformation
   */

  BAL_InitTransformation( &theTrsf );
  if ( par.input_real_transformation[0] != '\0' ) {
    if ( BAL_ReadTransformation( &theTrsf, par.input_real_transformation ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_real_transformation );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read real transformation...\n", 0 );
    }
    theTrsf.transformation_unit = REAL_UNIT;
  }
  else if ( par.input_voxel_transformation[0] != '\0' ) {
    if ( BAL_ReadTransformation( &theTrsf, par.input_voxel_transformation ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_voxel_transformation );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read voxel transformation...\n", 0 );
    }
    theTrsf.transformation_unit = VOXEL_UNIT;
  }
  else {
    API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "no input transformation...\n", 0 );
  }



  /* reading / building template
   * 1. from template image
   * 2. from parameters
   */
  if ( BAL_InitImage( &theTemplate, (char*)NULL, 1, 1, 1, 1, TYPE_UNKNOWN ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
  }

  if ( par.template_name[0] != '\0' ) {
    if ( BAL_ReadImage( &theTemplate, par.template_name, 0 ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
          fprintf( stderr, "... can not read template image '%s'\n", par.template_name );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read template image...\n", 0 );
    }
    vtfree( theTemplate.array ); theTemplate.array  = NULL;
    vtfree( theTemplate.data );  theTemplate.data  = NULL;
    vtfree( theTemplate.name );  theTemplate.name  = NULL;
    ptrTemplate = &theTemplate;
  }
  else if ( par.template_dim.x > 0 && par.template_dim.y > 0 ) {
    if ( par.template_dim.z > 0 ) {
      if ( BAL_InitImage( &theTemplate, (char*)NULL, par.template_dim.x, par.template_dim.y,
                                 par.template_dim.z, 1, TYPE_UNKNOWN ) != 1 ) {
          BAL_FreeTransformation( &theTrsf );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
        }
    }
    else {
      if ( BAL_InitImage( &theTemplate, (char*)NULL, par.template_dim.x, par.template_dim.y,
                               1, 1, TYPE_UNKNOWN ) != 1 ) {
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
      }
    }
    if ( par.template_voxel.x > 0.0 ) theTemplate.vx = par.template_voxel.x;
    if ( par.template_voxel.y > 0.0 ) theTemplate.vy = par.template_voxel.y;
    if ( par.template_voxel.z > 0.0 ) theTemplate.vz = par.template_voxel.z;
    if ( BAL_SetImageVoxelSizes( &theTemplate, theTemplate.vx, theTemplate.vy, theTemplate.vz ) != 1 ) {
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to initialize template image voxel sizes ...\n", 0 );
    }
    ptrTemplate = &theTemplate;
  }
  else if ( BAL_IsTransformationVectorField( &theTrsf) == 1 ) {
    if ( BAL_InitScalarImageFromImage( &theTemplate, (char*)NULL, &(theTrsf.vx), TYPE_UNKNOWN ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to initialize template image ...\n", 0 );
    }
    ptrTemplate = &theTemplate;
  }



  /************************************************************
   *
   * grey level images
   *
   ************************************************************/

  if ( par.input_image_0[0] != '\0' && par.input_image_1[0] != '\0' ) {

    if ( par.output_image_0[0] == '\0' && par.output_image_1[0] == '\0'
         && par.output_image[0] == '\0' ) {
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        if ( _verbose_ )
          fprintf( stderr, "... there are input grey-level images but no output images\n" );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "weird case...\n", 0 );
    }

    /* reading input images
     */
    if ( BAL_ReadImage( &theImage0, par.input_image_0, 0 ) != 1 ) {
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_image_0 );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read input image #0...\n", 0 );
    }
    ptrImage0 = &theImage0;

    if ( BAL_ReadImage( &theImage1, par.input_image_1, 0 ) != 1 ) {
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_image_1 );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read input image #1...\n", 0 );
    }
    ptrImage1 = &theImage1;


    /* change transformation unit if necessary
     */

    if ( theTrsf.transformation_unit == VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToRealUnit( ptrImage0, ptrImage1, &theTrsf, &theTrsf ) != 1 ) {
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to change transformation unit...\n", 0 );
      }
    }

    if ( ptrTemplate != (bal_image*)NULL )
      ptrTmp = ptrTemplate;
    else
      ptrTmp = &theImage1;


    /* allocating output images
     */
    if ( par.output_image_0[0] != '\0' ) {
      if ( BAL_AllocImageFromImage( &resImage0, (char*)NULL,
                                        ptrTmp, theImage0.type ) != 1 ) {
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to allocate output image #0...\n", 0 );
      }
      ptrResImage0 = &resImage0;
    }

    if ( par.output_image_1[0] != '\0' ) {
      if ( BAL_AllocImageFromImage( &resImage1, (char*)NULL,
                                        ptrTmp, theImage1.type ) != 1 ) {
        if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to allocate output image #1...\n", 0 );
      }
      ptrResImage1 = &resImage1;
    }

    if ( par.output_image[0] != '\0' ) {
      if ( BAL_AllocImageFromImage( &resImage, (char*)NULL,
                                        ptrTmp, theImage1.type ) != 1 ) {
        if ( ptrResImage1 != (bal_image*)NULL ) BAL_FreeImage( &resImage1 );
        if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to allocate output image...\n", 0 );
      }
      ptrResImage = &resImage;
    }

#define _GREYLEVEL_DESALLOCATIONS {                                    \
  if ( ptrResImage != (bal_image*)NULL ) BAL_FreeImage( &resImage );   \
  if ( ptrResImage1 != (bal_image*)NULL ) BAL_FreeImage( &resImage1 ); \
  if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 ); \
  if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );    \
  if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );    \
  BAL_FreeImage( &theTemplate );                                       \
  BAL_FreeTransformation( &theTrsf );                                  \
}

    if ( par.nimages > 0 ) {

      if ( par.write_extremities == 0 ) {
        imin = 1;
        imax = par.nimages;
      }
      else {
        imin = 0;
        imax = par.nimages+1;
      }

      for ( i=imin; i<=imax; i++ ) {

        t = (float)(i)/(float)(par.nimages+1);

        if ( BAL_InterpolateGreyLevelImages( ptrImage0, ptrImage1,
                                             ptrResImage0, ptrResImage1, ptrResImage,
                                             &theTrsf, t, par.interpolation ) != 1 ) {
          _GREYLEVEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when computing image at t=%f\n", par.index );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error during computation...\n", 0 );
        }

        if ( par.output_image_0[0] != '\0' ) {
          if ( BAL_WriteImage( ptrResImage0, _outputname( tmpname,  par.output_image_0, i, t ) ) != 1 ) {
            _GREYLEVEL_DESALLOCATIONS;
            if ( _verbose_ )
              fprintf( stderr, "... error when writing image '%s'\n", _outputname( tmpname,  par.output_image_0, i, t ) );
            API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing ouput image #0...\n", 0 );
          }
        }

        if ( par.output_image_1[0] != '\0' ) {
          if ( BAL_WriteImage( ptrResImage1, _outputname( tmpname,  par.output_image_1, i, t ) ) != 1 ) {
            _GREYLEVEL_DESALLOCATIONS;
            if ( _verbose_ )
              fprintf( stderr, "... error when writing image '%s'\n", _outputname( tmpname,  par.output_image_1, i, t ) );
            API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output image #1...\n", 0 );
          }
        }

        if ( par.output_image[0] != '\0' ) {
          if ( BAL_WriteImage( ptrResImage, _outputname( tmpname,  par.output_image, i, t ) ) != 1 ) {
            _GREYLEVEL_DESALLOCATIONS;
            if ( _verbose_ )
              fprintf( stderr, "... error when writing image '%s'\n", _outputname( tmpname,  par.output_image, i, t )  );
            API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output image...\n", 0 );
          }
        }

      }

    }
    else if ( 0.0 <= par.index && par.index <= 1.0 ) {

      if ( BAL_InterpolateGreyLevelImages( ptrImage0, ptrImage1,
                                           ptrResImage0, ptrResImage1, ptrResImage,
                                           &theTrsf, par.index, par.interpolation ) != 1 ) {
        _GREYLEVEL_DESALLOCATIONS;
        if ( _verbose_ )
          fprintf( stderr, "... error when computing image at t=%f\n", par.index );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error during computation...\n", 0 );
      }

      if ( par.output_image_0[0] != '\0' ) {
        if ( BAL_WriteImage( ptrResImage0, par.output_image_0 ) != 1 ) {
          _GREYLEVEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when writing image '%s'\n", par.output_image_0 );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing ouput image #0...\n", 0 );
        }
      }

      if ( par.output_image_1[0] != '\0' ) {
        if ( BAL_WriteImage( ptrResImage1, par.output_image_1 ) != 1 ) {
          _GREYLEVEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when writing image '%s'\n", par.output_image_1 );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output image #1...\n", 0 );
        }
      }

      if ( par.output_image[0] != '\0' ) {
        if ( BAL_WriteImage( ptrResImage, par.output_image ) != 1 ) {
          _GREYLEVEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when writing image '%s'\n", par.output_image );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output image...\n", 0 );
        }
      }

    }
    else {
      _GREYLEVEL_DESALLOCATIONS;
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "weird case...\n", 0 );
    }

    if ( ptrResImage != (bal_image*)NULL ) BAL_FreeImage( &resImage );
    if ( ptrResImage1 != (bal_image*)NULL ) BAL_FreeImage( &resImage1 );
    if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 );
    if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
    if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
    ptrImage0 = (bal_image*)NULL;
    ptrImage1 = (bal_image*)NULL;
    ptrResImage0 = (bal_image*)NULL;
    ptrResImage1 = (bal_image*)NULL;
    ptrResImage = (bal_image*)NULL;

    time_inter = _GetTime();
    clock_inter = _GetClock();

    if ( par.print_time >= 2 ) {
      fprintf( stderr, "grey-level interpolation: elapsed (real) time = %f\n", time_inter - time_init );
      fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_inter - clock_init );
      fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_inter - clock_init)/(time_inter - time_init) );
    }
  }



  /************************************************************
   *
   * labels images
   *
   ************************************************************/

  if ( par.input_label_0[0] != '\0' && par.input_label_1[0] != '\0'
       && par.input_label_correspondence[0] != '\0'
       && ( par.output_label_0[0] != '\0' || par.output_label_1[0] != '\0' ) ) {

    if ( par.output_label_0[0] == '\0' && par.output_label_1[0] == '\0' ) {
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        if ( _verbose_ )
          fprintf( stderr, "... there are input label images but no output images\n" );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "weird case...\n", 0 );
    }

    if ( par.input_label_correspondence[0] == '\0' ) {
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "label correspondence file is mandatory...\n", 0 );
    }

    /* reading input images
     */
    if ( BAL_ReadImage( &theImage0, par.input_label_0, 0 ) != 1 ) {
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_label_0 );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read input label image #0...\n", 0 );
    }
    ptrImage0 = &theImage0;

    if ( BAL_ReadImage( &theImage1, par.input_label_1, 0 ) != 1 ) {
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_label_1 );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read input label image #1...\n", 0 );
    }
    ptrImage1 = &theImage1;


    /* reading cell correspondences
     * 1. cell correspondences from files
     * 2. build list of cells
     * 3. add neighbors from cell correspondences to cell lists
     */
    BAL_InitCellCorrespondence( &theCellCorrespondence );

    if ( BAL_ReadCorrespondenceList( &(theCellCorrespondence.cliques),
                                     par.input_label_correspondence ) != 1 ) {
      if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
        fprintf( stderr, "... error when reading '%s'\n", par.input_label_correspondence );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to read correspondence list...\n", 0 );
    }

    if ( BAL_FillCellList( &theImage0, theCellCorrespondence.list0 ) != 1 ) {
      BAL_FreeCellCorrespondence( &theCellCorrespondence );
      if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to build cell list from labels #0 ...\n", 0 );
    }

    if ( BAL_FillCellList( &theImage1, theCellCorrespondence.list1 ) != 1 ) {
      BAL_FreeCellCorrespondence( &theCellCorrespondence );
      if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to build cell list from labels #1 ...\n", 0 );
    }

    if ( BAL_CrossFillCellCorrespondence( &theCellCorrespondence ) != 1 ) {
      BAL_FreeCellCorrespondence( &theCellCorrespondence );
      if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
      if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
      BAL_FreeImage( &theTemplate );
      BAL_FreeTransformation( &theTrsf );
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to build cell correspondences ...\n", 0 );
    }

    if ( par.relabel ) {
      if ( BAL_RelabelCells( &theImage0, &theImage1, &theCellCorrespondence  ) != 1 ) {
        BAL_FreeCellCorrespondence( &theCellCorrespondence );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to relabel cell correspondences ...\n", 0 );
      }
    }

    if ( _verbose_ ) {
      BAL_CheckCellCorrespondence( stderr, &theCellCorrespondence );
    }


    /* change transformation unit if necessary
     */

    if ( theTrsf.transformation_unit == VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToRealUnit( ptrImage0, ptrImage1, &theTrsf, &theTrsf ) != 1 ) {
        BAL_FreeCellCorrespondence( &theCellCorrespondence );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to change transformation unit...\n", 0 );
      }
    }

    if ( ptrTemplate != (bal_image*)NULL )
      ptrTmp = ptrTemplate;
    else
      ptrTmp = &theImage1;


    /* allocating output images
     */
    if ( par.output_label_0[0] != '\0' ) {
      if ( BAL_AllocImageFromImage( &resImage0, (char*)NULL,
                                        ptrTmp, theImage0.type ) != 1 ) {
        BAL_FreeCellCorrespondence( &theCellCorrespondence );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to allocate output label image #0...\n", 0 );
      }
      ptrResImage0 = &resImage0;
    }

    if ( par.output_label_1[0] != '\0' ) {
      if ( BAL_AllocImageFromImage( &resImage1, (char*)NULL,
                                        ptrTmp, theImage1.type ) != 1 ) {
        if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 );
        BAL_FreeCellCorrespondence( &theCellCorrespondence );
        if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
        if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
        BAL_FreeImage( &theTemplate );
        BAL_FreeTransformation( &theTrsf );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "unable to allocate output label image #1...\n", 0 );
      }
      ptrResImage1 = &resImage1;
    }

#define _LABEL_DESALLOCATIONS {                                        \
  if ( ptrResImage1 != (bal_image*)NULL ) BAL_FreeImage( &resImage1 ); \
  if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 ); \
  BAL_FreeCellCorrespondence( &theCellCorrespondence );                \
  if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );    \
  if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );    \
  BAL_FreeImage( &theTemplate );                                       \
  BAL_FreeTransformation( &theTrsf );                                  \
}

    if ( par.nimages > 0 ) {

      if ( par.write_extremities == 0 ) {
        imin = 1;
        imax = par.nimages;
      }
      else {
        imin = 0;
        imax = par.nimages+1;
      }

      for ( i=imin; i<=imax; i++ ) {

        t = (float)(i)/(float)(par.nimages+1);

        if ( BAL_InterpolateLabelImages( ptrImage0, ptrImage1,
                                         ptrResImage0, ptrResImage1,
                                         &theTrsf, &theCellCorrespondence,
                                         t, par.rmax ) != 1 ) {
          _LABEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when computing image at t=%f\n", par.index );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error during computation...\n", 0 );
        }

        if ( par.output_label_0[0] != '\0' ) {
          if ( BAL_WriteImage( ptrResImage0, _outputname( tmpname,  par.output_label_0, i, t ) ) != 1 ) {
            _LABEL_DESALLOCATIONS;
            if ( _verbose_ )
              fprintf( stderr, "... error when writing image '%s'\n", _outputname( tmpname,  par.output_label_0, i, t ) );
            API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing ouput label image #0...\n", 0 );
          }
        }

        if ( par.output_label_1[0] != '\0' ) {
          if ( BAL_WriteImage( ptrResImage1, _outputname( tmpname,  par.output_label_1, i, t ) ) != 1 ) {
            _LABEL_DESALLOCATIONS;
            if ( _verbose_ )
              fprintf( stderr, "... error when writing image '%s'\n", _outputname( tmpname,  par.output_label_1, i, t ) );
            API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output label image #1...\n", 0 );
          }
        }

      }

    }
    else if ( 0.0 <= par.index && par.index <= 1.0 ) {

      if ( BAL_InterpolateLabelImages( ptrImage0, ptrImage1,
                                       ptrResImage0, ptrResImage1,
                                       &theTrsf, &theCellCorrespondence,
                                       par.index, par.rmax ) != 1 ) {
        _LABEL_DESALLOCATIONS;
        if ( _verbose_ )
          fprintf( stderr, "... error when computing image at t=%f\n", par.index );
        API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error during computation...\n", 0 );
      }

      if ( par.output_label_0[0] != '\0' ) {
        if ( BAL_WriteImage( ptrResImage0, par.output_label_0 ) != 1 ) {
          _LABEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when writing image '%s'\n", par.output_label_0 );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing ouput label image #0...\n", 0 );
        }
      }

      if ( par.output_label_1[0] != '\0' ) {
        if ( BAL_WriteImage( ptrResImage1, par.output_label_1 ) != 1 ) {
          _LABEL_DESALLOCATIONS;
          if ( _verbose_ )
            fprintf( stderr, "... error when writing image '%s'\n", par.output_label_1 );
          API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "error when writing output label image #1...\n", 0 );
        }
      }

    }
    else {
      _LABEL_DESALLOCATIONS;
      API_ErrorParse_interpolateImages( _BaseName( argv[0] ), "weird case...\n", 0 );
    }

    if ( ptrResImage1 != (bal_image*)NULL ) BAL_FreeImage( &resImage1 );
    if ( ptrResImage0 != (bal_image*)NULL ) BAL_FreeImage( &resImage0 );
    BAL_FreeCellCorrespondence( &theCellCorrespondence );
    if ( ptrImage1 != (bal_image*)NULL ) BAL_FreeImage( &theImage1 );
    if ( ptrImage0 != (bal_image*)NULL ) BAL_FreeImage( &theImage0 );
    ptrImage0 = (bal_image*)NULL;
    ptrImage1 = (bal_image*)NULL;
    ptrResImage0 = (bal_image*)NULL;
    ptrResImage1 = (bal_image*)NULL;
    ptrResImage = (bal_image*)NULL;

    time_inter2 = _GetTime();
    clock_inter2 = _GetClock();

    if ( time_inter < 0.0 ) time_inter = time_init;
    if ( time_inter < 0.0 ) clock_inter = clock_init;

    if ( par.print_time >= 2 ) {
      fprintf( stderr, "label interpolation: elapsed (real) time = %f\n", time_inter2 - time_inter );
      fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_inter2 - clock_inter );
      fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_inter2 - clock_inter)/(time_inter2 - time_inter) );
    }
  }





  BAL_FreeImage( &theTemplate );
  BAL_FreeTransformation( &theTrsf );





  if ( par.trace_allocations ) {
    fprintfVtMallocTrace( stderr );
    clearVtMalloc();
  }

  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( argv[0] ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 0 );
}





/************************************************************
 *
 * static functions
 *
 ************************************************************/

static char * _outputname( char *tmp, char *format, int i, float t )
{
   char *proc = "_outputname";
   char *p;

   p = strstr( format, "%");

   if ( p == (char*)NULL ) return( format );

   while( (*p >= '0' && *p <= '9')
          || *p == '.'
          || *p == '%'
          || *p == '-' )
     p++;

   switch( *p ) {
   default :
     if ( _verbose_ )
       fprintf( stderr, "%s: format '%s' not recognized\n", proc, format );
     return( format );
   case 'd' :
     sprintf( tmp, format, i );
     return( tmp );
   case 'f' :
     sprintf( tmp, format, t );
     return( tmp );
   }

   if ( _verbose_ )
     fprintf( stderr, "%s: weird, should not reach this point\n", proc );
   return( format );
}






/************************************************************
 *
 * static functions
 *
 ************************************************************/


#ifdef UNUSED
static char *_Array2Str( int argc, char *argv[] )
{
  char *proc = "_Array2Str";
  int i, l;
  char *s, *t;

  if ( argc <= 1 || argv == (char**)NULL ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: no options in argv[]\n", proc );
    return( (char*)NULL );
  }

  /* there are argc-1 strings
   * compute the sum of string lengths from 1 to argc-1
   * + number of interval between successive strings (argc-2)
   * + 1 to add a trailing '\0'
   */
  for ( l=argc-1, i=1; i<argc; i++ ) {
    l += strlen( argv[i] );
  }

  s = (char*)vtmalloc( l * sizeof( char ), "s", proc );
  if ( s == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( (char*)NULL );
  }

  for ( t=s, i=1; i<argc; i++ ) {
    (void)strncpy( t, argv[i], strlen( argv[i] ) );
    t += strlen( argv[i] );
    if ( i < argc-1 ) {
      *t = ' ';
      t++;
    }
    else {
      *t = '\0';
    }
  }

  return( s );
}
#endif


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



static double _GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}



static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
