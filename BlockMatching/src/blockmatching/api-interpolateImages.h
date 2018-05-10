/*************************************************************************
 * api-interpolateImages.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu  2 jul 2015 09:43:10 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _api_interpolateimages_h_
#define _api_interpolateimages_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-image.h>



typedef struct lineCmdParamInterpolateImages {

  /* image names
   * image_0 is the 'floating' image (blockmatching convention)
   * image_1 is the 'reference' image (blockmatching convention)
   * and the transformation goes from the reference to the floating
   */
  char input_image_0[STRINGLENGTH];
  char input_image_1[STRINGLENGTH];
  char input_label_0[STRINGLENGTH];
  char input_label_1[STRINGLENGTH];

  char input_label_correspondence[STRINGLENGTH];

  /* output images (image names in case of one index
   * or format in case of interval)
   */
  char output_image[STRINGLENGTH];
  char output_image_0[STRINGLENGTH];
  char output_image_1[STRINGLENGTH];
  char output_label_0[STRINGLENGTH];
  char output_label_1[STRINGLENGTH];

  ImageType output_type;

  /* transformation
   */
  char input_real_transformation[STRINGLENGTH];
  char input_voxel_transformation[STRINGLENGTH];

  /* specific arguments
   */
  float index;
  int nimages;
  int write_extremities;
  int relabel;

  /* template
   */
  char template_name[STRINGLENGTH];
  bal_integerPoint template_dim;
  bal_doublePoint template_voxel;

  /* specific arguments
   */
  float rmax;
  enumTransformationInterpolation interpolation;

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;
  int trace_allocations;

} lineCmdParamInterpolateImages;










extern char *API_Help_interpolateImages( int h );

extern void API_ErrorParse_interpolateImages( char *program, char *str, int flag );

extern void API_InitParam_interpolateImages( lineCmdParamInterpolateImages *par );

extern void API_PrintParam_interpolateImages( FILE *theFile, char *program,
                                         lineCmdParamInterpolateImages *par,
                                         char *str );

extern void API_ParseParam_interpolateImages( int firstargc, int argc, char *argv[],
                                 lineCmdParamInterpolateImages *p );



#ifdef __cplusplus
}
#endif

#endif
