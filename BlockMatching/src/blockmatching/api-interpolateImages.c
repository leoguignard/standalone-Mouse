/*************************************************************************
 * api-interpolateImages.c -
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
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <reech4x4.h>
#include <reech4x4-coeff.h>
#include <reech-def.h>
#include <vtmalloc.h>

#include <bal-interpolation.h>
#include <bal-transformation-inversion.h>

#include <api-interpolateImages.h>






static int _verbose_ = 1;
static int _debug_ = 0;


#ifdef UNUSED
static void _API_ParseParam_interpolateImages( char *str, lineCmdParamInterpolateImages *p );
#endif


/************************************************************
 *
 * main API
 *
 ************************************************************/









/************************************************************
 *
 * static functions
 *
 ************************************************************/


#ifdef UNUSED
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
#endif





/************************************************************
 *
 * help / documentation
 *
 ************************************************************/



static char *usage = "[-input0|-floating|-flo] %s\n\
 [-input1|-reference|-ref] %s\n\
 [-input0-label|-floating-label|-flo-label %s]\n\
 [-input1-label|-reference-label|-ref-label %s]\n\
 [-label-correspondence|-label-corr %s]\n\
 [-result[-format]|-res|-output[-format]|-format] %s\n\
 [-output0[-format] %s] [-output1[-format] %s]\n\
 [-output0-label[-format] %s] [-output1-label[-format] %s]\n\
 [-transformation |-trsf %s]\n\
 [-voxel-transformation |-voxel-trsf %s]\n\
 [-i|-index %f] [-n|-nimages %d]\n\
 [-write-extremities|-we] [-do-not-write-extremities|-nwe]\n\
 [-relabel]\n\
 [-template-image|-template|-dims|-t %s]\n\
 [-template-dimension[s]|-template-dim|-dimension[s]|-dim %d %d [%d]]\n\
 [-x %d] [-y %d] [-z %d]\n\
 [-template-voxel|-voxel-size|-voxel|-pixel|-vs %lf %lf [%lf]]\n\
 [-r|-rmax %f]\n\
 [-nearest|-linear|-cspline] [-interpolation nearest|linear|cspline]\n\
 [-inversion-error %lf] [-inversion-iteration %d]\n\
 [-inversion-derivation-sigma %lf]\n\
 [-inversion-initialization zero|forward] [-inversion-forward-sigma %lf]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallelism-type|-parallel-type default|none|openmp|omp|pthread|thread]\n\
 [-omp-scheduling|-omps default|static|dynamic-one|dynamic|guided]\n\
 [output-image-type | -type s8|u8|s16|u16...]\n\
 [-verbose|-v] [-no-verbose|-noverbose|-nv]\n\
 [-debug|-D] [-no-debug|-nodebug]\n\
 [-print-parameters|-param]\n\
 [-print-time|-time] [-no-time|-notime]\n\
 [-trace-memory|-memory] [-no-memory|-nomemory]\n\
 [-help|-h]";



static char *detail = "\
# input images / data\n\
 [-input0|-floating|-flo] %s  # first/left greylevel input image\n\
 [-input1|-reference|-ref] %s # second/right greylevel input image\n\
 [-input0-label|-floating-label|-flo-label %s]  # \n\
    first/left label input image\n\
 [-input1-label|-reference-label|-ref-label %s] # \n\
    second/right label input image\n\
 [-label-correspondence|-label-corr %s] # correspondence file\n\
    it contains lines of the forme 'l0 l1 - r0 r1 r2'\n\
    where the li are a set of labels from the left image\n\
    that correspond to the set of rj labels of the right image\n\
# output images/formats\n\
 when '-index' is used, the output strings correspond to file names\n\
 while when '-nimages' is used, they correspond to file name formats\n\
 ie they contains a '%d' to describe a generic name a la printf\n\
 [-result[-format]|-res|-output[-format]|-format] %s # result image (format)\n\
    linear combination of deformed both left and right images\n\
 [-output0[-format] %s] # deformed left greylevel image (format)\n\
 [-output1[-format] %s] # deformed right greylevel image (format) \n\
 [-output0-label[-format] %s] # deformed left label image (format)\n\
 [-output1-label[-format] %s] # deformed right label image (format)\n\
# specific parameters\n\
 [-index|-i %f] # interpolate at 'index' position (between 0 and 1)\n\
    near 0 is close to left image while near 1 is close to right image\n\
    result is only one image\n\
 [-nimages|-n %d] # interpolate 'n' intermediary images.\n\
    Output image names are assumed to be in 'a la printf' format\n\
    indexes begins at 0 (with '-write-extremities') or 1\n\
    Eg: specifying '9' comes to divide the interval into 10 intervals\n\
 [-write-extremities|-we] # interpolate at extremities\n\
 [-do-not-write-extremities|-nwe] # \n\
 [-relabel] # relabel label image\n\
   in case of one-to-one mapping, ensure that cells have the same label\n\
# transformations\n\
 -transformation|-trsf %s # transformation to be applied\n\
    in 'real' coordinates. Goes from 'right image' to 'left image'.\n\
    If indicated, '-voxel-transformation' is ignored.\n\
 -voxel-transformation|-voxel-trsf %s:  transformation to be applied\n\
    in 'voxel' coordinates.\n\
# template image for output image geometry. If no information is given,\n\
    input image geometry is used\n\
 -template-image|-template|-dims|-t %s: template image\n\
    to set image geometry for result image\n\
 -template-dimension[s]|-template-dim|-dimension[s]|-dim %d %d [%d]: dimensions\n\
    of the result image\n\
 -x %d: dimension along X of the result image\n\
 -y %d: dimension along Y of the result image\n\
 -z %d: dimension along Z of the result image\n\
 -template-voxel|-voxel-size|-voxel|-pixel|-vs %lf %lf [%lf]:\n\
    voxel sizes of the result image\n\
# specific parameters\n\
 -r|-rmax %f: maximal radius for neighbor search\n\
   must be > 0.\n\
   specifying a very small rmax (<<1) allows to select only\n\
   compatible nearest interpolated points\n\
 -interpolation nearest|linear|cspline: selection of interpolation mode\n\
 -nearest: nearest neighor interpolation mode (for binary or lable images)\n\
 -linear: bi- or tri-linear interpolation\n\
 -cspline: cubic spline\n\
# inversion of vector field based transformations\n\
 -inversion-error %lf: absolute error (in real world unit) to determine convergence\n\
 -inversion-iteration %d: maximal number of iterations to reach convergence\n\
 -inversion-derivation-sigma %lf: standard deviation of the gaussian used to compute derivatives\n\
 -inversion-error-image %s: create an image of convergence defaults\n\
 -inversion-initialization %s:\n\
    zero:\n\
    forward: forward interpolation from input vector field\n\
 -inversion-forward-sigma %lf: gaussian kernel for interpolation\n\
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
  -no-verbose|-noverbose|-nv: no verboseness\n\
  -debug|-D: increase debug level\n\
  -no-debug|-nodebug: no debug indication\n\
  -print-parameters|-param:\n\
  -print-time|-time:\n\
  -no-time|-notime:\n\
  -trace-memory|-memory:\n\
  -no-memory|-nomemory:\n\
  -h: print option list\n\
  -help: print option list + details\n\
###########################################################\n\
 example of use: if a transformation is computed with 'blockmatching'\n\
 with 'blockmatching -flo image-0 -ref image-1 -res-trsf trsf -res result ...'\n\
 then 'interpolateImages [-flo] image-0 [-ref] image-1 [-res] result -trsf trsf ...'\n\
 allows to interpolate between image-0 and image-1 \n\
";





char *API_Help_interpolateImages( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_interpolateImages( char *program, char *str, int flag )
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



void API_InitParam_interpolateImages( lineCmdParamInterpolateImages *p )
{
    (void)strncpy( p->input_image_0, "\0", 1 );
    (void)strncpy( p->input_image_1, "\0", 1 );
    (void)strncpy( p->input_label_0, "\0", 1 );
    (void)strncpy( p->input_label_1, "\0", 1 );

    (void)strncpy( p->input_label_correspondence, "\0", 1 );

    (void)strncpy( p->output_image, "\0", 1 );
    (void)strncpy( p->output_image_0, "\0", 1 );
    (void)strncpy( p->output_image_1, "\0", 1 );
    (void)strncpy( p->output_label_0, "\0", 1 );
    (void)strncpy( p->output_label_1, "\0", 1 );

    p->output_type = TYPE_UNKNOWN;

    (void)strncpy( p->input_real_transformation, "\0", 1 );
    (void)strncpy( p->input_voxel_transformation, "\0", 1 );

    p->index = 0.5;
    p->nimages = 0;
    p->write_extremities = 0;
    p->relabel = 0;

    (void)strncpy( p->template_name, "\0", 1 );

    p->template_dim.x = 0;
    p->template_dim.y = 0;
    p->template_dim.z = 0;

    p->template_voxel.x = -1.0;
    p->template_voxel.y = -1.0;
    p->template_voxel.z = -1.0;

    p->rmax = -1.0;
    p->interpolation = LINEAR;

    p->print_lineCmdParam = 0;
    p->print_time = 0;
    p->trace_allocations = 0;
}





void API_PrintParam_interpolateImages( FILE *theFile, char *program,
                                  lineCmdParamInterpolateImages *p, char *str )
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


  fprintf( f, "# image names\n" );

  fprintf( f, "- input image #0 is " );
  if ( p->input_image_0 != (char*)NULL && p->input_image_0[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_image_0 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input image #1 is " );
  if ( p->input_image_1 != (char*)NULL && p->input_image_1[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_image_1 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input label image #0 is " );
  if ( p->input_label_0 != (char*)NULL && p->input_label_0[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_label_0 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input label image #1 is " );
  if ( p->input_label_1 != (char*)NULL && p->input_label_1[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_label_1 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- input label correspondence is " );
  if ( p->input_label_correspondence != (char*)NULL && p->input_label_correspondence[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_label_correspondence );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image (format) is " );
  if ( p->output_image != (char*)NULL && p->output_image[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_image );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image #0 (format) is " );
  if ( p->output_image_0 != (char*)NULL && p->output_image_0[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_image_0 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image #1 (format) is " );
  if ( p->output_image_1 != (char*)NULL && p->output_image_1[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_image_1 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output label #0 (format) is " );
  if ( p->output_label_0 != (char*)NULL && p->output_label_0[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_label_0 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output label #1 (format) is " );
  if ( p->output_label_1 != (char*)NULL && p->output_label_1[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_label_1 );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# transformation names\n" );

  fprintf( f, "- transformation to be applied (real units) is " );
  if ( p->input_real_transformation != (char*)NULL && p->input_real_transformation[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_real_transformation );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- transformation to be applied (voxel units) is " );
  if ( p->input_voxel_transformation != (char*)NULL && p->input_voxel_transformation[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_voxel_transformation );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "# specific parameter\n" );

  fprintf( f, "- intermediary image index = %f\n", p->index );
  fprintf( f, "- intermediary image number = %d\n", p->nimages );
  fprintf( f, "- interpolate at extremities = %d\n", p->write_extremities );
  fprintf( f, "- relabel label images = %d\n", p->relabel );

  fprintf( f, "# template for output image geometry\n" );

  fprintf( f, "- template image is " );
  if ( p->template_name != (char*)NULL && p->template_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->template_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- template image dimensions are [%d %d %d]\n",
           p->template_dim.x, p->template_dim.y, p->template_dim.z );

  fprintf( f, "- template image voxel sizes are [%f %f %f]\n",
           p->template_voxel.x, p->template_voxel.y, p->template_voxel.z );

  fprintf( f, "# specific parameter\n" );

  fprintf( f, "- rmax = %f\n", p->rmax );
  fprintf( f, "- interpolation mode = " );
  switch ( p->interpolation ) {
  default :      fprintf( f, "unknown\n" ); break;
  case NEAREST : fprintf( f, "nearest point\n" ); break;
  case LINEAR :  fprintf( f, "bi- or tri-linear\n" ); break;
  case CSPLINE : fprintf( f, "cubic spline\n" ); break;
  }

  fprintf( f, "# general image related parameters\n" );

  fprintf( f, "- output image type = " );
  switch ( p->output_type ) {
  default :     fprintf( f, "TYPE_UNKNOWN\n" ); break;
  case SCHAR :  fprintf( f, "SCHAR\n" ); break;
  case UCHAR :  fprintf( f, "UCHAR\n" ); break;
  case SSHORT : fprintf( f, "SSHORT\n" ); break;
  case USHORT : fprintf( f, "USHORT\n" ); break;
  case UINT :   fprintf( f, "UINT\n" ); break;
  case SINT :   fprintf( f, "INT\n" ); break;
  case ULINT :  fprintf( f, "ULINT\n" ); break;
  case FLOAT :  fprintf( f, "FLOAT\n" ); break;
  case DOUBLE : fprintf( f, "DOUBLE\n" ); break;
  }

  fprintf( f, "# misc\n" );
  fprintf( f, "- p->print_lineCmdParam =  %d\n", p->print_lineCmdParam );
  fprintf( f, "- p->print_time =  %d\n", p->print_time );
  fprintf( f, "- p->trace_allocations =  %d\n", p->trace_allocations );
  fprintf( f, "\n" );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/


#ifdef UNUSED
static void _API_ParseParam_interpolateImages( char *str, lineCmdParamInterpolateImages *p )
{
  char *proc = "_API_ParseParam_interpolateImages";
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

  API_ParseParam_interpolateImages( 0, argc, argv, p );

  vtfree( argv );
}
#endif




static int _n_call_parse_ = 0;

void API_ParseParam_interpolateImages( int firstargc, int argc, char *argv[],
                                  lineCmdParamInterpolateImages *p )
{
  int i;
  int input0isread = 0;
  int input1isread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  int status;
  int maxchunks;
  int o=0, s=0, r=0;
  int iterations;
  double error, sigma;

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
            API_ErrorParse_interpolateImages( (char*)NULL, "'-' is not handled...\n", 0 );
          }

          /* input image names
           */
          else if ( (strcmp ( argv[i], "-input0") == 0 && argv[i][7] == '\0')
                    || (strcmp ( argv[i], "-floating") == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-flo") == 0 && argv[i][4] == '\0') ) {
            if ( input0isread == 1 )
              API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input0, input #0 already read", 0 );
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input0", 0 );
            (void)strcpy( p->input_image_0, argv[i] );
            input0isread = 1;
          }
          else if ( (strcmp ( argv[i], "-input1") == 0 && argv[i][7] == '\0')
                    || (strcmp ( argv[i], "-reference") == 0 && argv[i][10] == '\0')
                    || (strcmp ( argv[i], "-ref") == 0 && argv[i][4] == '\0') ) {
            if ( input1isread == 1 )
              API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input1, input #1 already read", 0 );
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input1", 0 );
            (void)strcpy( p->input_image_1, argv[i] );
            input1isread = 1;
          }
          else if ( (strcmp ( argv[i], "-input0-label") == 0 && argv[i][13] == '\0')
                    || (strcmp ( argv[i], "-floating-label") == 0 && argv[i][15] == '\0')
                    || (strcmp ( argv[i], "-flo-label") == 0 && argv[i][10] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input0-label", 0 );
            (void)strcpy( p->input_label_0, argv[i] );
          }
          else if ( (strcmp ( argv[i], "-input1-label") == 0 && argv[i][13] == '\0')
                    || (strcmp ( argv[i], "-reference-label") == 0 && argv[i][16] == '\0')
                    || (strcmp ( argv[i], "-ref-label") == 0 && argv[i][10] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -input1-label", 0 );
            (void)strcpy( p->input_label_1, argv[i] );
          }
          else if ( (strcmp ( argv[i], "-label-correspondence") == 0 )
                    || (strcmp ( argv[i], "-label-corr") == 0 && argv[i][11] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -label-correspondence", 0 );
            (void)strcpy( p->input_label_correspondence, argv[i] );
          }


          /* output image names
           */
          else if ( (strcmp ( argv[i], "-result") == 0 && argv[i][7] == '\0')
                    || (strcmp ( argv[i], "-res") == 0 && argv[i][4] == '\0')
                    || (strcmp ( argv[i], "-result-format") == 0 && argv[i][14] == '\0')
                    || (strcmp ( argv[i], "-format") == 0 && argv[i][7] == '\0')
                    || (strcmp ( argv[i], "-output") == 0 && argv[i][7] == '\0')
                    || (strcmp ( argv[i], "-output-format") == 0 && argv[i][14] == '\0') ) {
            if ( outputisread == 1 )
              API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output(-format), output already read", 0 );
            i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output(-format)", 0 );
            (void)strcpy( p->output_image, argv[i] );
            outputisread = 1;
          }
          else if ( (strcmp ( argv[i], "-output0") == 0 && argv[i][8] == '\0')
                     || (strcmp ( argv[i], "-output0-format") == 0 && argv[i][15] == '\0') ) {             i++;
             if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output0(-format)", 0 );
             (void)strcpy( p->output_image_0, argv[i] );
          }
          else if ( (strcmp ( argv[i], "-output1") == 0 && argv[i][8] == '\0')
                     || (strcmp ( argv[i], "-output1-format") == 0 && argv[i][15] == '\0') ) {             i++;
             if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output1(-format)", 0 );
             (void)strcpy( p->output_image_1, argv[i] );
          }
          else if ( (strcmp ( argv[i], "-output0-label") == 0 && argv[i][14] == '\0')
                     || (strcmp ( argv[i], "-output0-label-format") == 0 && argv[i][21] == '\0') ) {             i++;
             if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output0-label(-format)", 0 );
             (void)strcpy( p->output_label_0, argv[i] );
          }
          else if ( (strcmp ( argv[i], "-output1-label") == 0 && argv[i][14] == '\0')
                    || (strcmp ( argv[i], "-output1-label-format") == 0 && argv[i][21] == '\0') ) {             i++;
            if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -output1-label(-format)", 0 );
            (void)strcpy( p->output_label_1, argv[i] );
          }


          /* transformation names
           */

          else if ( strcmp ( argv[i], "-transformation") == 0
                    || (strcmp ( argv[i], "-trsf") == 0 && argv[i][5] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -transformation", 0 );
                 (void)strcpy( p->input_real_transformation, argv[i] );
          }
          else if ( strcmp ( argv[i], "-voxel-transformation") == 0
                    || (strcmp ( argv[i], "-voxel-trsf") == 0 && argv[i][11] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -voxel-transformation", 0 );
                 (void)strcpy( p->input_voxel_transformation, argv[i] );
          }

          /* specific arguments
           */

          else if ( (strcmp ( argv[i], "-i") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-index") == 0 && argv[i][6] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -index", 0 );
                 status = sscanf( argv[i], "%f", &p->index );
                 if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -index ...\n", 0 );
          }
          else if ( (strcmp ( argv[i], "-n") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-nimages") == 0 && argv[i][8] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -nimages", 0 );
                 status = sscanf( argv[i], "%d", &p->nimages );
                 if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -nimages ...\n", 0 );
          }
          else if ( (strcmp ( argv[i], "-we") == 0 && argv[i][3] == '\0')
                    || strcmp ( argv[i], "-write-extremities") == 0 ) {
            p->write_extremities = 1;
          }
          else if ( (strcmp ( argv[i], "-nwe") == 0 && argv[i][4] == '\0')
                    || strcmp ( argv[i], "-do-not-write-extremities") == 0 ) {
            p->write_extremities = 0;
          }
          else if ( (strcmp ( argv[i], "-relabel") == 0 && argv[i][8] == '\0') ) {
            p->relabel = 1;
          }

          /* template
           */

          else if ( strcmp ( argv[i], "-template-image") == 0
                    || (strcmp ( argv[i], "-template") == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-image", 0 );
                 (void)strcpy( p->template_name, argv[i] );
          }
          else if ( strcmp ( argv[i], "-template-dimensions") == 0
                    || strcmp ( argv[i], "-template-dimension") == 0
                    || strcmp ( argv[i], "-template-dim") == 0
                    || strcmp ( argv[i], "-dimensions") == 0
                    || strcmp ( argv[i], "-dimension") == 0
                    || (strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-dimensions %d", 0 );
            status = sscanf( argv[i], "%d", &(p->template_dim.x) );
            if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-dimensions %d", 0 );
            i ++;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-dimensions %d %d", 0 );
            status = sscanf( argv[i], "%d", &(p->template_dim.y) );
            if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-dimensions %d %d", 0 );
            i ++;
            if ( i >= argc) p->template_dim.z = 1;
            else {
              status = sscanf( argv[i], "%d", &(p->template_dim.z) );
              if ( status <= 0 ) {
                i--;
                p->template_dim.z = 1;
              }
            }
          }
          else if ( strcmp ( argv[i], "-x") == 0 && argv[i][2] == '\0' ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -x ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->template_dim.x) );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -x ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-y") == 0 && argv[i][2] == '\0' ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -y ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->template_dim.y) );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -y ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-z") == 0 && argv[i][2] == '\0' ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -z ...\n", 0 );
              status = sscanf( argv[i], "%d", &(p->template_dim.z) );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -z ...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-template-voxel") == 0
                    || strcmp ( argv[i], "-voxel-size") == 0
                    || (strcmp (argv[i], "-voxel" ) == 0 && argv[i][6] == '\0')
                    || (strcmp (argv[i], "-pixel" ) == 0 && argv[i][6] == '\0')
                    || (strcmp (argv[i], "-vs" ) == 0 && argv[i][3] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-voxel %lf", 0 );
            status = sscanf( argv[i], "%lf", &(p->template_voxel.x) );
            if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-voxel %lf", 0 );
            i ++;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-voxel %lf %lf", 0 );
            status = sscanf( argv[i], "%lf", &(p->template_voxel.y) );
            if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -template-voxel %lf %lf", 0 );
            i ++;
            if ( i >= argc) p->template_voxel.z = 1;
            else {
              status = sscanf( argv[i], "%lf", &(p->template_voxel.z) );
              if ( status <= 0 ) {
                i--;
                p->template_voxel.z = 1;
              }
            }
          }

          /* specific parameters
           */

          else if ( (strcmp ( argv[i], "-r") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-rmax") == 0 && argv[i][5] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -rmax", 0 );
                 status = sscanf( argv[i], "%f", &p->rmax );
                 if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -rmax ...\n", 0 );
          }

          else if ( strcmp ( argv[i], "-interpolation" ) == 0 ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -interpolation...\n", 0 );
            if ( strcmp ( argv[i], "nearest" ) == 0 ) {
               p->interpolation = NEAREST;
            }
            else if ( strcmp ( argv[i], "linear" ) == 0 ) {
              p->interpolation = LINEAR;
            }
            else if ( strcmp ( argv[i], "cspline" ) == 0 ) {
                p->interpolation = CSPLINE;
            }
            else {
              fprintf( stderr, "unknown interpolation mode: '%s'\n", argv[i] );
              API_ErrorParse_interpolateImages( (char*)NULL, "parsing -interpolation ...\n", 0 );
            }
          }

          else if ( strcmp ( argv[i], "-nearest" ) == 0 ) {
            p->interpolation = NEAREST;
          }
          else if ( strcmp ( argv[i], "-linear" ) == 0 ) {
            p->interpolation = LINEAR;
          }
          else if ( strcmp ( argv[i], "-cspline" ) == 0 ) {
              p->interpolation = CSPLINE;
          }

          /* vector field transformation, inversion monitoring
           */
          else if ( strcmp ( argv[i], "-inversion-error") == 0 && argv[i][16] == '\0') {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-error ...\n", 0 );
              status = sscanf( argv[i], "%lf", &error );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-error ...\n", 0 );
              if ( error > 0.0 ) BAL_SetErrorMaxForVectorFieldInversionInBalTransformationInversion( error );
          }
          else if ( strcmp ( argv[i], "-inversion-iteration") == 0
                    || strcmp ( argv[i], "-inversion-iterations") == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-iterations ...\n", 0 );
              status = sscanf( argv[i], "%d", &iterations );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-iterations ...\n", 0 );
              if ( iterations >= 0 ) BAL_SetIterationsMaxForVectorFieldInversionInBalTransformationInversion( iterations );
          }
          else if ( strcmp ( argv[i], "-inversion-derivation-sigma") == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-derivation-sigma ...\n", 0 );
              status = sscanf( argv[i], "%lf", &sigma );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-derivation-sigma ...\n", 0 );
              if ( sigma > 0.0 ) BAL_SetDerivationSigmaForVectorFieldInversionInBalTransformationInversion( sigma );
          }
          else if ( strcmp ( argv[i], "-inversion-initialization" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-initialization ...\n", 0 );
             if ( strcmp ( argv[i], "zero" ) == 0 ) {
                 BAL_SetInitializationForVectorFieldInversionInBalTransformationInversion( ZERO );
             }
             else if ( strcmp ( argv[i], "forward" ) == 0 ) {
                 BAL_SetInitializationForVectorFieldInversionInBalTransformationInversion( FORWARD_INTERPOLATION );
             }
             else {
               fprintf( stderr, "unknown initialization type: '%s'\n", argv[i] );
               API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-initialization ...\n", 0 );
             }
          }
          else if ( strcmp ( argv[i], "-inversion-forward-sigma") == 0 ) {
              i++;
              if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-forward-sigma ...\n", 0 );
              status = sscanf( argv[i], "%lf", &sigma );
              if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -inversion-forward-sigma ...\n", 0 );
              if ( sigma > 0.0 ) BAL_SetForwardSigmaForVectorFieldInversionInBalTransformationInversion( sigma );
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
             if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
               API_ErrorParse_interpolateImages( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
               API_ErrorParse_interpolateImages( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
             }
          }

          /* general image related parameters
           */

          else if ( strcmp ( argv[i], "-r" ) == 0 && argv[i][2] == '\0' ) {
             r = 1;
          }
          else if ( strcmp ( argv[i], "-s" ) == 0 && argv[i][2] == '\0' ) {
             s = 1;
          }
          else if ( strcmp ( argv[i], "-o" ) == 0 && argv[i][2] == '\0' ) {
             i += 1;
             if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -o...\n", 0 );
             status = sscanf( argv[i],"%d",&o );
             if ( status <= 0 ) API_ErrorParse_interpolateImages( (char*)NULL, "parsing -o...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_interpolateImages( (char*)NULL, "parsing -type...\n", 0 );
            if ( strcmp ( argv[i], "s8" ) == 0 && argv[i][2] == '\0' ) {
               p->output_type = SCHAR;
            }
            else if ( strcmp ( argv[i], "u8" ) == 0 && argv[i][2] == '\0' ) {
               p->output_type = UCHAR;
            }
            else if ( strcmp ( argv[i], "s16" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SSHORT;
            }
            else if ( strcmp ( argv[i], "u16" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = USHORT;
            }
            else if ( strcmp ( argv[i], "s32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SINT;
            }
            else if ( strcmp ( argv[i], "u32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = UINT;
            }
            else if ( strcmp ( argv[i], "s64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = SLINT;
            }
            else if ( strcmp ( argv[i], "u64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = ULINT;
            }
            else if ( strcmp ( argv[i], "r32" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = FLOAT;
            }
            else if ( strcmp ( argv[i], "r64" ) == 0 && argv[i][3] == '\0' ) {
              p->output_type = DOUBLE;
            }
            else {
              API_ErrorParse_interpolateImages( (char*)NULL, "parsing -type...\n", 0 );
            }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_interpolateImages( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_interpolateImages( (char*)NULL, (char*)NULL, 0);
          }
          else if ( strcmp ( argv[i], "-verbose" ) == 0
                    || (strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _verbose_ <= 0 ) _verbose_ = 1;
              else                  _verbose_ ++;
              if ( 0 ) {
                incrementVerboseInReech4x4();
                incrementVerboseInReech4x4Coeff();
                incrementVerboseInReechDef();
              }
              BAL_IncrementVerboseInBalTransformation();
              BAL_IncrementVerboseInBalTransformationInversion();
              BAL_IncrementVerboseInBalInterpolation();
            }
          }
          else if ( strcmp ( argv[i], "-no-verbose" ) == 0
                    || strcmp ( argv[i], "-noverbose" ) == 0
                    || (strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0') ) {
              _verbose_ = 0;
              setVerboseInReech4x4( 0 );
              setVerboseInReech4x4Coeff( 0 );
              setVerboseInReechDef( 0 );
              BAL_SetVerboseInBalTransformation( 0 );
              BAL_SetVerboseInBalTransformationInversion( 0 );
              BAL_SetVerboseInBalInterpolation( 0 );
          }
          else if ( (strcmp ( argv[i], "-debug" ) == 0 && argv[i][6] == '\0')
                    || (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
            if ( _n_call_parse_ == 1 ) {
              if ( _debug_ <= 0 ) _debug_ = 1;
              else                _debug_ ++;
              BAL_IncrementDebugInBalInterpolation();
            }
          }
          else if ( (strcmp ( argv[i], "-no-debug" ) == 0 && argv[i][9] == '\0')
                    || (strcmp ( argv[i], "-nodebug" ) == 0 && argv[i][8] == '\0') ) {
              _debug_ = 0;
              BAL_SetDebugInBalInterpolation( 0 );
          }

          else if ( strcmp ( argv[i], "-print-parameters" ) == 0
                    || (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
             p->print_lineCmdParam = 1;
          }

          else if ( strcmp ( argv[i], "-print-time" ) == 0
                     || (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
             if ( p->print_time <= 0 )
               p->print_time = 1;
             else
               p->print_time ++;
          }
          else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')
                      || (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
             p->print_time = 0;
          }

          else if ( strcmp ( argv[i], "-trace-memory" ) == 0
                     || (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
             if ( _n_call_parse_ == 1 ) {
               incrementTraceInVtMalloc( );
               if ( p->trace_allocations  <= 0 ) p->trace_allocations  = 1;
               else                              p->trace_allocations  ++;
             }
             if ( 0 ) setParallelism( _NO_PARALLELISM_ );
          }
          else if ( (strcmp ( argv[i], "-nomemory" ) == 0 && argv[i][9] == '\0')
                      || (strcmp ( argv[i], "-no-memory" ) == 0 && argv[i][10] == '\0') ) {
             setTraceInVtMalloc( 0 );
          }

          /* unknown option
           */
          else {
              sprintf(text,"unknown option %s\n",argv[i]);
              API_ErrorParse_interpolateImages( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_interpolateImages( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( input0isread == 0 ) {
              (void)strcpy( p->input_image_0, argv[i] );
              input0isread = 1;
          }
          else if ( input1isread == 0 ) {
              (void)strcpy( p->input_image_1, argv[i] );
              input1isread = 1;
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_image, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_interpolateImages( (char*)NULL, "too many file names ...\n", 0 );
          }
      }
  }

  /* output image type
   */
  if ( (o != 0) || (s != 0) || (r != 0) ) {
    if ( (o == 1) && (s == 1) && (r == 0) ) p->output_type = SCHAR;
    else if ( (o == 1) && (s == 0) && (r == 0) ) p->output_type = UCHAR;
    else if ( (o == 2) && (s == 0) && (r == 0) ) p->output_type = USHORT;
    else if ( (o == 2) && (s == 1) && (r == 0) ) p->output_type = SSHORT;
    else if ( (o == 4) && (s == 1) && (r == 0) ) p->output_type = SINT;
    else if ( (o == 0) && (s == 0) && (r == 1) ) p->output_type = FLOAT;
    else {
        API_ErrorParse_interpolateImages( (char*)NULL, "unable to determine uotput image type ...\n", 0 );
    }
  }

}
