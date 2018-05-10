/*************************************************************************
 * api-pointmatching.c -
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <bal-field.h>
#include <bal-field-tools.h>
#include <bal-point.h>
#include <bal-transformation-tools.h>

#include <api-pointmatching.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_pointmatching( char *str, lineCmdParamPointmatching *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/



bal_transformation *API_pointmatching( bal_doublePointList *floatingPoints,
                                       bal_doublePointList *referencePoints,
                                       bal_image *imtemplate_transformation,
                                       char *result_residual,
                                       char *param_str_1,
                                       char *param_str_2 )
{
  char *proc = "API_pointmatching";
  bal_transformation *resultTransformation = (bal_transformation*)NULL;
  FIELD theField;
  int max_n_data, min_n_data;
  lineCmdParamPointmatching par;



  /* parameter initialization
   */
  API_InitParam_pointmatching( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_pointmatching( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_pointmatching( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_pointmatching( stderr, proc, &par, (char*)NULL );



  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  /* allocate computed transformation
   */
  resultTransformation = (bal_transformation*)vtmalloc( sizeof(bal_transformation),
                                                        "resultTransformation", proc );
  if ( resultTransformation == (bal_transformation*)NULL ) {
    if ( _verbose_ )
        fprintf( stderr, "%s: can not allocate initial result transformation container\n", proc );
    return( (bal_transformation*)NULL );
  }

  BAL_InitTransformation( resultTransformation );

  if ( BAL_AllocTransformation( resultTransformation,
                                par.transformation_type,
                                imtemplate_transformation ) != 1 ) {
    vtfree( resultTransformation );
    if ( _verbose_ )
        fprintf( stderr, "%s: can not allocate initial result transformation\n", proc );
    return( (bal_transformation*)NULL );
  }



  /* read points
   */
  max_n_data = min_n_data = floatingPoints->n_data;

  if ( floatingPoints->n_data != referencePoints->n_data ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: list of points have different lengths\n", proc );
      fprintf( stderr, "   floating list contains %d points\n", floatingPoints->n_data );
      fprintf( stderr, "   reference list contains %d points\n", referencePoints->n_data );
    }
    if ( max_n_data < referencePoints->n_data )
      max_n_data = referencePoints->n_data;
    if ( min_n_data > referencePoints->n_data )
      min_n_data = referencePoints->n_data;
    if ( _verbose_ ) {
      fprintf( stderr, "   uses the first %d of both lists\n", min_n_data );
    }
  }

  if ( BAL_AllocateField( &theField, max_n_data ) != 1 ) {
    BAL_FreeTransformation( resultTransformation );
    vtfree( resultTransformation );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating displacement field\n", proc );
    return( (bal_transformation*)NULL );
  }

  if ( imtemplate_transformation != (bal_image *)NULL ) {
    if ( BAL_CopyImageGeometryToField( imtemplate_transformation, &theField ) != 1 ) {
        BAL_FreeField( &theField );
        BAL_FreeTransformation( resultTransformation );
        vtfree( resultTransformation );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when setting displacement field from image\n", proc );
        return( (bal_transformation*)NULL );
    }
  }
  else {
      if ( BAL_SetFieldVoxelSizes( &theField, referencePoints->vx, referencePoints->vy, referencePoints->vz ) != 1 ) {
          BAL_FreeField( &theField );
          BAL_FreeTransformation( resultTransformation );
          vtfree( resultTransformation );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when setting displacement field voxel sizes\n", proc );
          return( (bal_transformation*)NULL );
      }
  }


  /* won't work in the points are in voxel units, with different
     voxel sizes
     build a displacement field where
     origin = reference point
     displacement = floating point - reference point
  */
  if ( BAL_ComputePairingFieldFromDoublePointList( &theField, floatingPoints,
                                                   referencePoints ) != 1 ) {
    BAL_FreeField( &theField );
    BAL_FreeTransformation( resultTransformation );
    vtfree( resultTransformation );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing displacement field\n", proc );
    return( (bal_transformation*)NULL );
  }



  /* transformation computation procedure assumes that
     the points are in voxel units (in the reference referential)
     Well, this is due to historical reasons, which is definitively not
     a good reason
  */
  if ( theField.unit == REAL_UNIT ) {
    size_t i;
    for ( i=0; i<theField.n_computed_pairs; i++ ) {
      theField.data[i].origin.x /= theField.vx;
      theField.data[i].origin.y /= theField.vy;
      theField.data[i].origin.z /= theField.vz;
      theField.data[i].vector.x /= theField.vx;
      theField.data[i].vector.y /= theField.vy;
      theField.data[i].vector.z /= theField.vz;
    }
    theField.unit = VOXEL_UNIT;
  }



  /* computes transformation
   * the resulting transformation allows
   * - to reasmple a floating image into a reference one,
   * - to transform the reference points into the floating frame
   * it goes from the reference frame to the floating one
   */
  if ( BAL_ComputeIncrementalTransformation( resultTransformation,
                                             &theField, &(par.estimator) ) != 1 ) {
    BAL_FreeField( &theField );
    BAL_FreeTransformation( resultTransformation );
    vtfree( resultTransformation );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing transformation\n", proc );
    return( (bal_transformation*)NULL );
  }



  /* computes residuals and write residuals on disk
   * such an operation should not be in this API
   * (I/O operations are supposed to take place in the calling procedure)
   * however, this is the most convenient.
   */

  if ( result_residual != (char*)NULL && result_residual[0] != '\0' ) {
    if ( BAL_ComputeTransformationResiduals( resultTransformation,
                                             &theField ) != 1 ) {
      BAL_FreeField( &theField );
      BAL_FreeTransformation( resultTransformation );
      vtfree( resultTransformation );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when computing residuals\n", proc );
      return( (bal_transformation*)NULL );
    }
    {
      FILE *f = fopen( result_residual, "w" );
      size_t i;
      double sum = 0.0;
      if ( f == (FILE*)NULL ) {
        BAL_FreeField( &theField );
        BAL_FreeTransformation( resultTransformation );
        vtfree( resultTransformation );
        if ( _verbose_ )
          fprintf( stderr, "%s: error opening residual file '%s'\n", proc, result_residual );
        return( (bal_transformation*)NULL );
      }
      for ( i=0; i<theField.n_computed_pairs; i++ )
        fprintf( f, "%f\n", sqrt( theField.data[i].error ) );
      fprintf( f, "\n" );
      for ( i=0; i<theField.n_selected_pairs; i++ )
        sum += sqrt( theField.pointer[i]->error );
      fprintf( f, "# ----------------------------------------\n" );
      fprintf( f, "# average on %lu points\n", theField.n_selected_pairs );
      fprintf( f, "# ----------------------------------------\n" );
      fprintf( f, "%f\n", sum / (double)theField.n_selected_pairs );
      fclose( f );
    }
  }

  BAL_FreeField( &theField );

  return( resultTransformation );
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
static char *usage = "-reference %s -floating %s\n\
 [-result-transformation|-res-trsf %s]\n\
 [-result-voxel-transformation|-res-voxel-trsf %s]\n\
 [-residuals %s]\n\
 [-transformation-type|-transformation|-trsf-type %s]\n\
 [-estimator-type|-estimator|-es-type wlts|lts|wls|ls]\n\
 [-lts-fraction %lf] [-lts-deviation %f] [-lts-iterations %d]\n\
 [-fluid-sigma|-lts-sigma[-ll|-hl] %lf %lf %lf]\n\
 [-vector-propagation-type|-propagation-type|-propagation direct|skiz]\n\
 [-vector-propagation-distance|-propagation-distance|-pdistance %f]\n\
 [-vector-fading-distance|-fading-distance|-fdistance %f]\n\
 [-points-unit|-unit voxel|real]\n\
 [-points-voxel|-points-vs %f %f [%f]]\n\
 [-floating-voxel %f %f [%f]]\n\
 [-reference-voxel %f %f [%f]]\n\
 [-template|-t|-dims %s]\n\
 [-dim %d %d [%d] | [-x %d] [-y %d] [-z %d]]\n\
 [-voxel | -pixel | -vs %f %f [%f] | [-vx %f] [-vy %f] [-vz %f] ]\n\
 [-command-line %s]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-nv|-noverbose] [-debug|-D] [-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-notime]\n\
 [-help|-h]";




/*---------------- longueur maximale de ligne --------------------------------*/
/*----------------------------------------------------------------------------*/
static char *detail = "\
###\n\
Computes a transformation with the paired points from the reference to the\n\
floating. This is somehow the 'blockmatching' way, the resulting\n\
transformation will allow to resample the floating image (from which the\n\
floating points are drawn) onto the reference one.\n\
### File names ###\n\
-reference|-ref %s  # name of the reference point list (still points)\n\
   There is a line in the file per point of the form 'x y z'\n\
   Lines beginning with either '#' or '\%' are considered as comments\n\
-floating|-flo %s   # name of the point list to be registered (floating points)\n\
   points are assumed to be paired according to their rank\n\
[-result-transformation|-res-trsf %s] # name of the result transformation\n\
  in 'real' coordinates. Goes from 'reference' to 'floating', ie allows to \n\
  resample 'floating' in the geometry of 'reference'.\n\
  If indicated, '-result-voxel-transformation' is ignored.\n\
[-result-voxel-transformation|-res-voxel-trsf %s] # name of the result\n\
  transformation in 'voxel' coordinates.\n\
[-residuals %s]     # name of the output file for the residual values\n\
### transformation type ###\n\
[-transformation-type|-transformation|-trsf-type %s] # transformation type\n\
  translation2D, translation3D, translation-scaling2D, translation-scaling3D,\n\
  rigid2D, rigid3D, rigid, similitude2D, similitude3D, similitude,\n\
  affine2D, affine3D, affine, vectorfield2D, vectorfield3D, vectorfield, vector\n\
### transformation estimation ###\n\
[-estimator-type|-estimator|-es-type %s] # transformation estimator\n\
  wlts: weighted least trimmed squares\n\
  lts: least trimmed squares\n\
  wls: weighted least squares\n\
  ls: least squares\n\
[-lts-fraction %lf] # for trimmed estimations, fraction of pairs that are kept\n\
[-lts-deviation %lf] # for trimmed estimations, defines the threshold to discard\n\
  pairings, ie 'average + this_value * standard_deviation'\n\
[-lts-iterations %d] # for trimmed estimations, the maximal number of iterations\n\
### transformation estimation: vectorfield estimation dedicated parameters ###\n\
[-fluid-sigma|-lts-sigma] %lf %lf %lf] # sigma for fluid regularization,\n\
  ie field interpolation and regularization for pairings (only for vector field)\n\
[-vector-propagation-type|-propagation-type|-propagation direct|skiz] # \n\
  direct: exact propagation (but slow)\n\
  skiz: approximate propagation (but faster)\n\
[-vector-propagation-distance|-propagation-distance|-pdistance %f] # \n\
  distance propagation of initial pairings (ie displacements)\n\
  this implies the same displacement for the spanned sphere\n\
  (only for vectorfield). Distance is in world units (not voxels).\n\
[-vector-fading-distance|-fading-distance|-fdistance %f] # \n\
  area of fading for initial pairings (ie displacements)\n\
  this allows progressive transition towards null displacements\n\
  and thus avoid discontinuites. Distance is in world units (not voxels).\n\
### geometry information for points\n\
[-points-unit|-unit voxel|real] # units of the points (default is REAL)\n\
[-points-voxel|-points-vs %f %f [%f]] #\n\
  allows to specify the voxel sizes when points are given in voxel units\n\
  same voxel sizes is assumed for both the floating and reference points\n\
[-floating-voxel %f %f [%f]] #\n\
[-reference-voxel %f %f [%f]] #\n\
### deformation field geometry / voxel size ###\n\
  specifying the output image/vector field with the above options\n\
  is mandatory for this kind of transformation\n\
[-template|-t|-dims %s] # template image for the dimensions\n\
                         of the output image\n\
[-dim %d %d [%d]]      # output image dimensions\n\
[-x %d]                # X dimension of the ouput image\n\
[-y %d]                # Y dimension of the ouput image\n\
[-z %d]                # Z dimension of the ouput image\n\
[-voxel|-pixel|-vs %f %f [%f]]    # output image voxel sizes\n\
[-vx %f]               # voxel size along X dimension of the ouput image\n\
[-vy %f]               # voxel size along Y dimension of the ouput image\n\
[-vz %f]               # voxel size along Z dimension of the ouput image\n\
### misc writing stuff ###\n\
[-command-line %s]           # write the command line\n\
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





char *API_Help_pointmatching( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_pointmatching( char *program, char *str, int flag )
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



void API_InitParam_pointmatching( lineCmdParamPointmatching *p )
{
    (void)strncpy( p->floating_points, "\0", 1 );
    (void)strncpy( p->reference_points, "\0", 1 );

    (void)strncpy( p->result_real_transformation, "\0", 1 );
    (void)strncpy( p->result_voxel_transformation, "\0", 1 );

    (void)strncpy( p->result_residual, "\0", 1 );

    /* parameters for  matching
     */
    p->transformation_type = RIGID_3D;
    BAL_InitEstimator( &(p->estimator) );
    p->estimator.type = TYPE_LS;

    /* geometry information for points
     */
    p->points_unit = REAL_UNIT;
    p->floating_voxel.x = -1.0;
    p->floating_voxel.y = -1.0;
    p->floating_voxel.z = -1.0;
    p->reference_voxel.x = -1.0;
    p->reference_voxel.y = -1.0;
    p->reference_voxel.z = -1.0;

    /* geometry information for vectorfield estimate
     */

    (void)strncpy( p->template_name, "\0", 1 );

    p->dim.x = 0;
    p->dim.y = 0;
    p->dim.z = 0;

    p->voxel.x = -1.0;
    p->voxel.y = -1.0;
    p->voxel.z = -1.0;

    /* general parameters
     */
    (void)strncpy( p->command_line_file, "\0", 1 );
    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_pointmatching( FILE *theFile, char *program,
                                  lineCmdParamPointmatching *p, char *str )
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


  fprintf( f, "# file names\n" );

  fprintf( f, "- floating points file is " );
  if ( p->floating_points != (char*)NULL && p->floating_points[0] != '\0' )
    fprintf( f, "'%s'\n", p->floating_points );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- reference points file is " );
  if ( p->reference_points != (char*)NULL && p->reference_points[0] != '\0' )
    fprintf( f, "'%s'\n", p->reference_points );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- p->result_real_transformation = " );
  if ( p->result_voxel_transformation[0] != '\0' )
      fprintf( f, "'%s'\n", p->result_voxel_transformation );
  else
      fprintf( f, "NULL\n" );

  fprintf( f, "- p->result_voxel_transformation = " );
  if ( p->result_voxel_transformation[0] != '\0' )
      fprintf( f, "'%s'\n", p->result_voxel_transformation );
  else
      fprintf( f, "NULL\n" );

  fprintf( f, "- p->result_residual = " );
  if ( p->result_residual[0] != '\0' )
      fprintf( f, "'%s'\n", p->result_residual );
  else
      fprintf( f, "NULL\n" );

  fprintf( f, "# parameters for matching\n" );

  BAL_PrintTypeTransformation( f, p->transformation_type, "p->transformation_type = " );
  fprintf( f, "- estimator =" );
  BAL_PrintEstimator(f, &(p->estimator) );

  fprintf( f, "# geometry information for points\n" );

  fprintf( f, "- p->points_unit = " );
  switch( p->points_unit ) {
  case UNDEF_UNIT : fprintf( f, "UNDEF_UNIT\n" ); break;
  case VOXEL_UNIT : fprintf( f, "VOXEL_UNIT\n" ); break;
  case REAL_UNIT :  fprintf( f, "REAL_UNIT\n" ); break;
  }
  BAL_PrintDoublePoint( f, &(p->floating_voxel),  "- p->floating_voxel  = " );
  BAL_PrintDoublePoint( f, &(p->reference_voxel), "- p->reference_voxel = " );

  fprintf( f, "# geometry information for vectorfield estimate\n" );

  fprintf( f, "- p->template_name = " );
  if ( p->template_name[0] != '\0' )
      fprintf( f, "'%s'\n", p->template_name );
  else
      fprintf( f, "NULL\n" );

  BAL_PrintIntegerPoint( f, &(p->dim), "- p->dim = " );
  BAL_PrintDoublePoint( f, &(p->voxel), "- p->voxel = " );


  fprintf( f, "# general parameters\n" );

  fprintf( f, "- p->command_line_file = " );
  if ( p->command_line_file[0] != '\0' )
      fprintf( f, "'%s'\n", p->command_line_file );
  else
      fprintf( f, "NULL\n" );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_pointmatching( char *str, lineCmdParamPointmatching *p )
{
  char *proc = "_API_ParseParam_pointmatching";
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

  API_ParseParam_pointmatching( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_pointmatching( int firstargc, int argc, char *argv[],
                                  lineCmdParamPointmatching *p )
{
  int i;
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

        /* point file names
         */
        if ( (strcmp ( argv[i], "-reference") == 0 && argv[i][10] == '\0')
             || (strcmp ( argv[i], "-ref") == 0 && argv[i][4] == '\0') ) {
           i++;
           if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -reference", 0 );
           (void)strcpy( p->reference_points, argv[i] );
        }

        else if ( (strcmp ( argv[i], "-floating") == 0 && argv[i][9] == '\0')
                  || (strcmp ( argv[i], "-flo") == 0 && argv[i][4] == '\0') ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -floating", 0 );
          (void)strcpy( p->floating_points, argv[i] );
        }

        /* transformation file names
         */
        else if ( strcmp ( argv[i], "-result-transformation" ) == 0
                  || strcmp ( argv[i], "-res-trsf" ) == 0 ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -result-transformation", 0 );
          (void)strcpy( p->result_real_transformation, argv[i] );
        }

        else if ( strcmp ( argv[i], "-result-voxel-transformation" ) == 0
                  || strcmp ( argv[i], "-res-voxel-trsf" ) == 0 ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -result-voxel-transformation", 0 );
          (void)strcpy( p->result_voxel_transformation, argv[i] );
        }

        /* residuals
         */
        else if ( strcmp ( argv[i], "-residuals" ) == 0 ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -residuals", 0 );
          (void)strcpy( p->result_residual, argv[i] );
        }

        /* transformation type
         */
        else  if ( strcmp ( argv[i], "-transformation-type" ) == 0
         || strcmp ( argv[i], "-transformation" ) == 0
         || strcmp ( argv[i], "-trsf-type" ) == 0 ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-transformation-type", 0 );
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
                      || strcmp ( argv[i], "vector" ) == 0 ) {
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
              API_ErrorParse_pointmatching( (char*)NULL, "-transformation-type", 0 );
            }
        }

        /* estimator definition and computation
         */
        else if ( strcmp ( argv[i], "-estimator-type") == 0
                  || strcmp ( argv[i], "-estimator") == 0
                  || strcmp ( argv[i], "-es-type") == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-estimator-type", 0 );
          if ( (strcmp ( argv[i], "ltsw" ) == 0 && argv[i][4] == '\0')
               || (strcmp ( argv[i], "wlts" ) == 0 && argv[i][4] == '\0') ) {
            p->estimator.type = TYPE_WLTS;
          }
          else if ( strcmp ( argv[i], "lts" ) == 0 && argv[i][3] == '\0' ) {
            p->estimator.type = TYPE_LTS;
          }
          else if ( (strcmp ( argv[i], "lsw" ) == 0 && argv[i][3] == '\0')
                    || (strcmp ( argv[i], "wls" ) == 0 && argv[i][3] == '\0') ) {
            p->estimator.type = TYPE_WLS;
          }
          else if ( strcmp ( argv[i], "ls" ) == 0 && argv[i][2] == '\0' ) {
            p->estimator.type = TYPE_LS;
          }
          else {
            fprintf( stderr, "unknown estimator type: '%s'\n", argv[i] );
            API_ErrorParse_pointmatching( (char*)NULL, "-estimator-type", 0 );
          }
        }

        else if ( strcmp ( argv[i], "-lts-fraction" ) == 0
                  || strcmp ( argv[i], "-lts-cut" ) == 0) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-lts-fraction", 0 );
          status = sscanf( argv[i], "%lf", &(p->estimator.retained_fraction) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "-lts-fraction", 0 );
        }

        else if ( strcmp ( argv[i], "-lts-deviation" ) == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-lts-deviation", 0 );
          status = sscanf( argv[i], "%lf", &(p->estimator.standard_deviation_threshold) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "-lts-deviation", 0 );
        }

        else if ( strcmp ( argv[i], "-lts-iterations" ) == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-lts-iterations", 0 );
          status = sscanf( argv[i], "%d", &(p->estimator.max_iterations) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "-lts-iterations", 0 );
        }

        else if ( (strcmp (argv[i], "-fluid-sigma" ) == 0 && argv[i][12] == '\0')
                  || (strcmp (argv[i], "-lts-sigma" ) == 0 && argv[i][10] == '\0') ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -lts-sigma %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->estimator.sigma.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -lts-sigma %lf", 0 );
          i ++;
          if ( i >= argc) {
            p->estimator.sigma.y = p->estimator.sigma.x;
            p->estimator.sigma.z = p->estimator.sigma.x;
          }
          else {
            status = sscanf( argv[i], "%lf", &(p->estimator.sigma.y) );
            if ( status <= 0 ) {
              i--;
              p->estimator.sigma.y = p->estimator.sigma.x;
              p->estimator.sigma.z = p->estimator.sigma.x;
            }
            else {
              i ++;
              if ( i >= argc) p->estimator.sigma.z = 0;
              else {
                status = sscanf( argv[i], "%lf", &(p->estimator.sigma.z) );
                if ( status <= 0 ) {
                  i--;
                  p->estimator.sigma.z = 0;
                }
              }
            }
          }
        }

        else if ( strcmp (argv[i], "-vector-propagation-type" ) == 0
                   || strcmp (argv[i], "-propagation-type" ) == 0
                   || ( strcmp (argv[i], "-propagation" ) == 0 && argv[i][12] == '\0') ) {
          i++;
          if ( strcmp (argv[i], "direct" ) == 0 ) {
            p->estimator.propagation.type = TYPE_DIRECT_PROPAGATION;
          }
          else if ( strcmp (argv[i], "skiz" ) == 0 ) {
            p->estimator.propagation.type = TYPE_SKIZ_PROPAGATION;
          }
          else {
            API_ErrorParse_pointmatching( (char*)NULL, "unknown propagation type for '-vector-propagation-type'", 0 );
          }
        }

        else if ( strcmp (argv[i], "-vector-propagation-distance" ) == 0
                   || strcmp (argv[i], "-propagation-distance" ) == 0
                   || (strcmp (argv[i], "-pdistance" ) == 0 && argv[i][10] == '\0') ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-vector-propagation-distance", 0 );
          status = sscanf( argv[i], "%f", &(p->estimator.propagation.constant) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "-vector-propagation-distance", 0 );
        }
        else if ( strcmp (argv[i], "-vector-fading-distance" ) == 0
                   || strcmp (argv[i], "-fading-distance" ) == 0
                   || (strcmp (argv[i], "-fdistance" ) == 0 && argv[i][10] == '\0') ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-vector-fading-distance", 0 );
          status = sscanf( argv[i], "%f", &(p->estimator.propagation.fading) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "-vector-fading-distance", 0 );
        }

        /* geometry information for points
         */
        else if ( strcmp ( argv[i], "-points-unit") == 0
                  || strcmp ( argv[i], "-units") == 0
                  || strcmp ( argv[i], "-unit") == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "-units", 0 );
          if ( (strcmp ( argv[i], "real" ) == 0 && argv[i][4] == '\0') ) {
            p->points_unit = REAL_UNIT;
          }
          else if ( (strcmp ( argv[i], "voxel" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "pixel" ) == 0 && argv[i][5] == '\0') ) {
            p->points_unit = VOXEL_UNIT;
          }
          else {
            fprintf( stderr, "unknown unis: '%s'\n", argv[i] );
            API_ErrorParse_pointmatching( (char*)NULL, "-points-unit", 0 );
          }
        }

        else if ( strcmp ( argv[i], "-points-voxel" ) == 0
                  || strcmp ( argv[i], "-points-pixel" ) == 0
                  || (strcmp ( argv[i], "-points-vs" ) == 0  && argv[i][3] == '\0')
                  || strcmp ( argv[i], "-floating-voxel" ) == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -points-voxel %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->floating_voxel.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -points-voxel %lf", 0 );
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -points-voxel %lf %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->floating_voxel.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -points-voxel %lf %lf", 0 );
          i ++;
          if ( i >= argc) p->floating_voxel.z = 1;
          else {
            status = sscanf( argv[i], "%lf", &(p->floating_voxel.z) );
            if ( status <= 0 ) {
              i--;
              p->floating_voxel.z = 1;
            }
          }
          if ( strcmp ( argv[i], "-points-voxel" ) == 0
               || strcmp ( argv[i], "-points-pixel" ) == 0
               || (strcmp ( argv[i], "-points-vs" ) == 0  && argv[i][3] == '\0') ) {
              p->reference_voxel.x = p->floating_voxel.x;
              p->reference_voxel.x = p->floating_voxel.x;
              p->reference_voxel.x = p->floating_voxel.x;
          }
        }

        else if ( strcmp ( argv[i], "-reference-voxel" ) == 0 ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -reference-voxel %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->reference_voxel.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -reference-voxel %lf", 0 );
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -reference-voxel %lf %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->reference_voxel.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -reference-voxel %lf %lf", 0 );
          i ++;
          if ( i >= argc) p->reference_voxel.z = 1;
          else {
            status = sscanf( argv[i], "%lf", &(p->reference_voxel.z) );
            if ( status <= 0 ) {
              i--;
              p->reference_voxel.z = 1;
            }
          }
        }

        /* geometry information for vectorfield estimate
         */
        else if ( strcmp ( argv[i], "-template") == 0
             || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
             || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -template", 0 );
          (void)strcpy( p->template_name, argv[i] );
        }

        else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -dim %d", 0 );
          status = sscanf( argv[i], "%d", &(p->dim.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -dim %d", 0 );
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -dim %d %d", 0 );
          status = sscanf( argv[i], "%d", &(p->dim.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -dim %d %d", 0 );
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

        else if ( strcmp ( argv[i], "-x") == 0 && argv[i][2] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -x...\n", 0 );
          status = sscanf( argv[i],"%d",&(p->dim.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -x...\n", 0 );
        }
        else if ( strcmp ( argv[i], "-y" ) == 0  && argv[i][2] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -y...\n", 0 );
          status = sscanf( argv[i],"%d",&(p->dim.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -y...\n", 0 );
        }
        else if ( strcmp ( argv[i], "-z" ) == 0  && argv[i][2] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -z...\n", 0 );
          status = sscanf( argv[i],"%d",&(p->dim.z) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -z...\n", 0 );
        }

        else if ( strcmp (argv[i], "-voxel" ) == 0
                  || strcmp ( argv[i], "-pixel" ) == 0
                  || (strcmp ( argv[i], "-vs" ) == 0  && argv[i][3] == '\0') ) {
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -voxel %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->voxel.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -voxel %lf", 0 );
          i ++;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -voxel %lf %lf", 0 );
          status = sscanf( argv[i], "%lf", &(p->voxel.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -voxel %lf %lf", 0 );
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

        else if ( strcmp ( argv[i], "-vx" ) == 0   && argv[i][3] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -vx...\n", 0 );
          status = sscanf( argv[i], "%lf", &(p->voxel.x) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -vx...\n", 0 );
        }
        else if ( strcmp ( argv[i], "-vy" ) == 0   && argv[i][3] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -vy...\n", 0 );
          status = sscanf( argv[i], "%lf", &(p->voxel.y) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -vy...\n", 0 );
        }
        else if ( strcmp ( argv[i], "-vz" ) == 0   && argv[i][3] == '\0' ) {
          i += 1;
          if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -vz...\n", 0 );
          status = sscanf( argv[i], "%lf", &(p->voxel.z) );
          if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -vz...\n", 0 );
        }

        /* misc writing stuff
         */
        else if ( strcmp ( argv[i], "-command-line") == 0 ) {
          i++;
          if ( i >= argc) API_ErrorParse_pointmatching( (char*)NULL, "parsing -command-line", 0 );
          (void)strcpy( p->command_line_file, argv[i] );
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
           if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
             API_ErrorParse_pointmatching( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
           }
        }

        else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
           i ++;
           if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -max-chunks ...\n", 0 );
           status = sscanf( argv[i], "%d", &maxchunks );
           if ( status <= 0 ) API_ErrorParse_pointmatching( (char*)NULL, "parsing -max-chunks ...\n", 0 );
           if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
        }

        else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                 ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
           i ++;
           if ( i >= argc)    API_ErrorParse_pointmatching( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
             API_ErrorParse_pointmatching( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
           }
        }

        /* general parameters
         */
        else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                  || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
           API_ErrorParse_pointmatching( (char*)NULL, (char*)NULL, 1);
        }
        else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                  || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
           API_ErrorParse_pointmatching( (char*)NULL, (char*)NULL, 0);
        }
        else if ( strcmp ( argv[i], "-verbose" ) == 0
                  || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
          if ( _n_call_parse_ == 1 ) {
            if ( _verbose_ <= 0 ) _verbose_ = 1;
            else                  _verbose_ ++;
            BAL_IncrementVerboseInBalPoint();
          }
        }
        else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                  || strcmp ( argv[i], "-noverbose" ) == 0
                  || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
            _verbose_ = 0;
            BAL_SetVerboseInBalPoint( 0 );
        }
        else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                  || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
          if ( _n_call_parse_ == 1 ) {
            if ( _debug_ <= 0 ) _debug_ = 1;
            else                _debug_ ++;
            BAL_IncrementDebugInBalPoint();
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
            API_ErrorParse_pointmatching( (char*)NULL, text, 0);
        }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          fprintf( stderr, "... parsing '%s'\n", argv[i] );
          API_ErrorParse_pointmatching( (char*)NULL, "too many file names ...\n", 0 );
       }
  }



}
