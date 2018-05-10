/*************************************************************************
 * api-applyTrsfToPoints.h -
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

#ifndef _api_applytrsftopoints_h_
#define _api_applytrsftopoints_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-image.h>



typedef struct lineCmdParamApplyTrsfToPoints {

  /* image names and output type
   */
  char input_name[STRINGLENGTH];
  char output_name[STRINGLENGTH];
  char input_real_transformation[STRINGLENGTH];

  /* specific arguments
   */

  /* ... */

  /* general parameters
   */
  int print_lineCmdParam;
  int print_time;

} lineCmdParamApplyTrsfToPoints;



extern int applyTrsfToPoint( bal_doublePoint thePt,
                             bal_doublePoint* resPt,
                             char *real_transformation_name,
                             int verbose,
                             int debug );



extern int API_applyTrsfToPoints( bal_doublePointList *floatingPoints,
                                  bal_doublePointList *referencePoints,
                                  bal_transformation *theTrsf,
                                  char *param_str_1,
                                  char *param_str_2 );



extern char *API_Help_applyTrsfToPoints( int h );

extern void API_ErrorParse_applyTrsfToPoints( char *program, char *str, int flag );

extern void API_InitParam_applyTrsfToPoints( lineCmdParamApplyTrsfToPoints *par );

extern void API_PrintParam_applyTrsfToPoints( FILE *theFile, char *program,
                                         lineCmdParamApplyTrsfToPoints *par,
                                         char *str );

extern void API_ParseParam_applyTrsfToPoints( int firstargc, int argc, char *argv[],
                                 lineCmdParamApplyTrsfToPoints *p );



#ifdef __cplusplus
}
#endif

#endif
