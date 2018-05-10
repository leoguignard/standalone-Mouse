/*************************************************************************
 * api-cropImage.h -
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

#ifndef _api_cropimage_h_
#define _api_cropimage_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-image.h>



typedef struct lineCmdParamCropImage {

  /* image names and output type
   */
  char input_name[STRINGLENGTH];
  char output_name[STRINGLENGTH];

  char output_real_transformation_name[STRINGLENGTH];
  char output_voxel_transformation_name[STRINGLENGTH];

  /* specific arguments
   */
  char template_name[STRINGLENGTH];

  bal_integerPoint origin;
  bal_integerPoint dim;
  bal_integerPoint slice;

  int originValue;
  int analyzeFiji;

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;

} lineCmdParamCropImage;





/* this API is kept for historical reasons
 * ie compatibility with STracking
 */
extern int cropImage( char *theim_name,
                      char *resim_name,
                      char *real_transformation_name,
                      char *voxel_transformation_name,
                      char *template_image_name, /* template */
                      int analyzeFiji,
                      bal_integerPoint origin,
                      bal_integerPoint dim,
                      bal_integerPoint slice,
                      int isDebug,
                      int isVerbose );

/* this API read the input files and then call the following one
 */
extern int API_INTERMEDIARY_cropImage( char *theim_name,
                                       char *resim_name,
                                       char *real_transformation_name,
                                       char *voxel_transformation_name,
                                       char *template_image_name,
                                       char *param_str_1, char *param_str_2 );




/* this API does the job, ie invert the read input
 * transformation into the allocated output one
 */
extern int API_cropImage( bal_image *image,
                             bal_image *imres,
                             char *param_str_1,
                             char *param_str_2 );



extern char *API_Help_cropImage( int h );

extern void API_ErrorParse_cropImage( char *program, char *str, int flag );

extern void API_InitParam_cropImage( lineCmdParamCropImage *par );

extern void API_PrintParam_cropImage( FILE *theFile, char *program,
                                         lineCmdParamCropImage *par,
                                         char *str );

extern void API_ParseParam_cropImage( int firstargc, int argc, char *argv[],
                                 lineCmdParamCropImage *p );



#ifdef __cplusplus
}
#endif

#endif
