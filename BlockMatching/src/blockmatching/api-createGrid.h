/*************************************************************************
 * api-createGrid.h -
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

#ifndef _api_creategrid_h_
#define _api_creategrid_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-image.h>


typedef enum typeOutputCreateGrid {
  _GRID,
  _MOSAIC
} typeOutputCreateGrid;


typedef struct lineCmdParamCreateGrid {

  /* image names and output type
   */
  char input_name[STRINGLENGTH];
  char output_name[STRINGLENGTH];

  ImageType output_type;

  /* geometry information for image template
   */

  char template_name[STRINGLENGTH];

  bal_integerPoint dim;
  bal_doublePoint voxel;

  /* specific arguments
   */

  typeOutputCreateGrid image;

  bal_integerPoint offset;
  bal_integerPoint spacing;
  float value;

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;

} lineCmdParamCreateGrid;



/* this one is kept for historical reasons,
 * ie Stracking compilation but should disappear
 */
extern int createGrid( char *resim_name,
                       char *template_image_name,
                       bal_integerPoint dim,
                       bal_doublePoint voxel,
                       bal_integerPoint offset,
                       bal_integerPoint spacing );


extern int API_INTERMEDIARY_createGrid( char *theim_name,
                                        char *resim_name,
                                        char *template_image_name,
                                        char *param_str_1, char *param_str_2 );

extern int API_createGrid( bal_image *image,
                             char *param_str_1,
                             char *param_str_2 );



extern char *API_Help_createGrid( int h );

extern void API_ErrorParse_createGrid( char *program, char *str, int flag );

extern void API_InitParam_createGrid( lineCmdParamCreateGrid *par );

extern void API_PrintParam_createGrid( FILE *theFile, char *program,
                                         lineCmdParamCreateGrid *par,
                                         char *str );

extern void API_ParseParam_createGrid( int firstargc, int argc, char *argv[],
                                 lineCmdParamCreateGrid *p );



#ifdef __cplusplus
}
#endif

#endif
