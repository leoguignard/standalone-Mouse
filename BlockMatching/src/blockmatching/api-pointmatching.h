/*************************************************************************
 * api-pointmatching.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 14 sep 2016 18:28:04 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _api_pointmatching_h_
#define _api_pointmatching_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <typedefs.h>

#include <bal-estimator.h>
#include <bal-transformation.h>

typedef struct lineCmdParamPointmatching {

  /* file names
     - images
     - transformations
  */
  char floating_points[STRINGLENGTH];
  char reference_points[STRINGLENGTH];

  char result_real_transformation[STRINGLENGTH];
  char result_voxel_transformation[STRINGLENGTH];

  char result_residual[STRINGLENGTH];

  /* parameters for  matching
   */
  enumTypeTransfo transformation_type;
  bal_estimator estimator;

  /* geometry information for points
   */
  enumUnitTransfo points_unit;
  bal_doublePoint floating_voxel;
  bal_doublePoint reference_voxel;

  /* geometry information for vectorfield estimate
   */

  char template_name[STRINGLENGTH];

  bal_integerPoint dim;
  bal_doublePoint voxel;

  /* ... */

  /* general parameters
   */
  char command_line_file[STRINGLENGTH];

  int print_lineCmdParam;
  int print_time;

} lineCmdParamPointmatching;



extern bal_transformation *API_pointmatching( bal_doublePointList *floatingPoints,
                                              bal_doublePointList *referencePoints,
                                              bal_image *imtemplate,
                                              char *result_residual,
                                              char *param_str_1,
                                              char *param_str_2 );



extern char *API_Help_pointmatching( int h );

extern void API_ErrorParse_pointmatching( char *program, char *str, int flag );

extern void API_InitParam_pointmatching( lineCmdParamPointmatching *par );

extern void API_PrintParam_pointmatching( FILE *theFile, char *program,
                                         lineCmdParamPointmatching *par,
                                         char *str );

extern void API_ParseParam_pointmatching( int firstargc, int argc, char *argv[],
                                 lineCmdParamPointmatching *p );



#ifdef __cplusplus
}
#endif

#endif
