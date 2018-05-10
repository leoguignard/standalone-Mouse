/*************************************************************************
 * api-copyTrsf.h -
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
 * to generate files:
 * sed -e "s/cropImage/execuTable/g" \
 *     -e "s/CropImage/ExecuTable/g" \
 *     -e "s/cropimage/executable/g" \
 *     [api-]cropImage.[c,h] > [api-]execTable.[c,h]
 *
 */

#ifndef _api_copytrsf_h_
#define _api_copytrsf_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-image.h>



typedef struct lineCmdParamCopyTrsf {

  /* image names and output type
   */
  char thetrsf_name[STRINGLENGTH];
  char restrsf_name[STRINGLENGTH];

  enumTypeTransfo transformation_type;

  char template_image_name[STRINGLENGTH];

  bal_integerPoint dim;
  bal_doublePoint voxel;

  enumUnitTransfo thetrsf_unit;
  enumUnitTransfo restrsf_unit;

  char floating_image_name[STRINGLENGTH];

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;

} lineCmdParamCopyTrsf;



/* this one is kept for historical reasons,
 * ie Stracking compilation but should disappear
 */
extern int copyTrsf( char* thetrsf_name,
                     char* restrsf_name,
                     enumUnitTransfo thetrsf_unit,
                     enumUnitTransfo restrsf_unit,
                     char *template_image_name,
                     char *floating_image_name,
                     bal_integerPoint dim,
                     bal_doublePoint voxel,
                     enumTypeTransfo transformation_type,
                     int isDebug,
                     int isVerbose );

extern int API_INTERMEDIARY_copyTrsf( char *thetrsf_name,
                                      char *restrsf_name,
                                      char *template_image_name,
                                      char *floating_image_name,
                                      char *param_str_1, char *param_str_2 );

/* this function should be (re)write to be called from
 * API_INTERMEDIARY_copyTrsf()
 * To be done
 */
extern int API_copyTrsf( bal_image *image,
                             bal_image *imres,
                             char *param_str_1,
                             char *param_str_2 );



extern char *API_Help_copyTrsf( int h );

extern void API_ErrorParse_copyTrsf( char *program, char *str, int flag );

extern void API_InitParam_copyTrsf( lineCmdParamCopyTrsf *par );

extern void API_PrintParam_copyTrsf( FILE *theFile, char *program,
                                         lineCmdParamCopyTrsf *par,
                                         char *str );

extern void API_ParseParam_copyTrsf( int firstargc, int argc, char *argv[],
                                 lineCmdParamCopyTrsf *p );



#ifdef __cplusplus
}
#endif

#endif
