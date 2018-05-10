/*************************************************************************
 * api-createGrid.c -
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
 * sed -e "s/createGrid/execuTable/g" \
 *     -e "s/CreateGrid/ExecuTable/g" \
 *     -e "s/creategrid/executable/g" \
 *     [api-]createGrid.[c,h] > [api-]execTable.[c,h]
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <bal-image-tools.h>

#include <api-createGrid.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_createGrid( char *str, lineCmdParamCreateGrid *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/

/* this one is kept for historical reasons,
 * ie Stracking compilation but should disappear
 */
int createGrid( char *resim_name,
                char *template_image_name,
                bal_integerPoint dim,
                bal_doublePoint voxel,
                bal_integerPoint offset,
                bal_integerPoint spacing )
{
  char *proc = "createGrid";
  char str[256], *s;

  if ( _verbose_ )
      fprintf( stderr, "Warning, '%s' is obsolete\n", proc );

  s = str;

  if ( dim.x >= 1 && dim.y >= 1 ) {
    sprintf( s, "-dim %d %d ", dim.x, dim.y );
    s = &(str[strlen(str)]);
    if ( dim.z >= 1 ) {
      sprintf( s, "%d ", dim.z );
      s = &(str[strlen(str)]);
    }
  }

  if ( voxel.x > 0.0 && voxel.y > 0.0 ) {
    sprintf( s, "-voxel %f %f ", voxel.x, voxel.y );
    s = &(str[strlen(str)]);
    if ( voxel.z > 0.0 ) {
      sprintf( s, "%f ", voxel.z );
      s = &(str[strlen(str)]);
    }
  }

  sprintf( s, "-offset %d %d %d", offset.x, offset.y, offset.z );
  s = &(str[strlen(str)]);

  sprintf( s, "-spacing %d %d %d", spacing.x, spacing.y, spacing.z );
  s = &(str[strlen(str)]);

  if ( API_INTERMEDIARY_createGrid( (char*)NULL,
                                   resim_name,
                                   template_image_name,
                                   s, (char*)NULL ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: some error occurs\n", proc );
      return( -1 );
  }

  return( 0 );
}





int API_INTERMEDIARY_createGrid( char *theim_name,
                                 char *resim_name,
                                 char *template_image_name,
                                 char *param_str_1, char *param_str_2 )
{
  char *proc = "API_INTERMEDIARY_createGrid";
  int v;
  bal_image theim;
  bal_image tmpim;
  bal_image *ptrim = (bal_image*)NULL;
  char *resname = (char*)NULL;
  lineCmdParamCreateGrid par;



  /* parameter initialization
   */
  API_InitParam_createGrid( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_createGrid( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_createGrid( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_createGrid( stderr, proc, &par, (char*)NULL );



  /* if there is an input image name
   * try to read it
   * unable verbose in bal-image.c to prevent error message
   */
  if ( theim_name != (char*)NULL && theim_name[0] != '\0' ) {
    v = BAL_GetVerboseInBalImage();
    BAL_SetVerboseInBalImage( 0 );
    if ( BAL_ReadImage( &theim, theim_name, 0 ) != 1 ) {
        if ( _verbose_ >= 3 )
          fprintf( stderr, "%s: can not read input image '%s'\n", proc, theim_name );
    }
    else {
      ptrim = &theim;
    }
    BAL_SetVerboseInBalImage( v );
  }

  /* do nothing if the input image was read
   * else, this may be the output image name
   * then initialize the resut image
   */
  if ( ptrim != (bal_image*)NULL ) {
      ;
  }
  else if ( template_image_name != (char*)NULL && template_image_name[0] != '\0' ) {
    if ( BAL_ReadImage( &tmpim, template_image_name, 0 ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: can not read template image '%s'\n", proc, template_image_name );
        return( -1 );
    }
    if ( tmpim.type == par.output_type ) {
      if ( BAL_FillImage( &tmpim, 0.0 ) != 1 ) {
        BAL_FreeImage( &tmpim );
        if ( _verbose_ )
          fprintf( stderr, "%s: can not initialize result image\n", proc );
        return( -1 );
      }
      ptrim = &tmpim;
    }
    else {
      if ( BAL_AllocImageFromImage( &theim, (char *)NULL, &tmpim, par.output_type ) != 1 ) {
        BAL_FreeImage( &tmpim );
        if ( _verbose_ )
          fprintf( stderr, "%s: can not allocate result image\n", proc );
        return( -1 );
      }
      BAL_FreeImage( &tmpim );
      ptrim = &theim;
    }
  }

  else if ( par.dim.x > 0 && par.dim.y > 0 ) {
    if ( par.dim.z > 0 ) {
      if ( BAL_AllocFullImage( &theim, (char*)NULL, par.dim.x, par.dim.y, par.dim.z, 1, 1.0, 1.0, 1.0, par.output_type ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate result image\n", proc );
        return( -1 );
      }
    }
    else {
      if ( BAL_AllocFullImage( &theim, (char*)NULL, par.dim.x, par.dim.y, 1, 1, 1.0, 1.0, 1.0, par.output_type ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate result image (dimz=1) \n", proc );
        return( -1 );
      }
    }
    if ( par.voxel.x > 0.0 ) theim.vx = par.voxel.x;
    if ( par.voxel.y > 0.0 ) theim.vy = par.voxel.y;
    if ( par.voxel.z > 0.0 ) theim.vz = par.voxel.z;
    if ( BAL_SetImageVoxelSizes( &theim, theim.vx, theim.vy, theim.vz ) != 1 ) {
      BAL_FreeImage( &theim );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to initialize result image geometry\n", proc );
      return( -1 );
    }
    ptrim = &theim;
  }

  else {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, unable to specify result image geometry\n", proc );
    return( -1 );
  }


  /* API call
   */
  if ( API_createGrid( ptrim, param_str_1, param_str_2 ) !=1 ) {
    BAL_FreeImage( ptrim );
    if ( _verbose_ )
      fprintf( stderr, "%s: can not draw grid \n", proc );
    return( -1 );
  }

  /* choose the output name
   */
  if ( theim_name != (char*)NULL && theim_name[0] != '\0' ) {
    if  (resim_name != (char*)NULL && resim_name[0] != '\0' ) {
      resname = resim_name;
    }
    else {
      resname = theim_name;
    }
  }
  else {
    resname = resim_name;
  }

  if ( BAL_WriteImage( ptrim, resname ) != 1 ) {
    BAL_FreeImage( ptrim );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write result image '%s'\n", proc, resname );
    return( -1 );
  }

  BAL_FreeImage( ptrim );
  return( 1 );
}





int API_createGrid( bal_image *image, char *param_str_1, char *param_str_2 )
{
  char *proc = "API_createGrid";
  lineCmdParamCreateGrid par;



  /* parameter initialization
   */
  API_InitParam_createGrid( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_createGrid( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_createGrid( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_createGrid( stderr, proc, &par, (char*)NULL );

  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  switch( par.image ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such output image type not handled yet\n", proc );
    return( -1 );
  case _GRID :
    if ( BAL_DrawGrid( image, &(par.offset), &(par.spacing), par.value ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: can not draw grid \n", proc );
      return( -1 );
    }
    break;
  case _MOSAIC :
    if ( BAL_DrawMosaic( image, &(par.spacing) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: can not draw mosaic \n", proc );
      return( -1 );
    }
    break;
  }

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
static char *usage = "[image-in] image-out\n\
 [-template|-t|-dims %s]\n\
 [-dim %d %d [%d]]\n\
 [-voxel | -pixel | -vs %f %f [%f]]\n\
 [-type-output|-type grid|mosaic]\n\
 [-offset %d %d [%d]]\n\
 [-spacing %d %d [%d]]\n\
 [-value %f]\n\
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
# ...\n\
### template for image creation ###\n\
[-template|-t|-dims %s] # template image for the dimensions\n\
                          of the output image\n\
[-dim %d %d [%d]]      # output image dimensions\n\
[-voxel|-pixel|-vs %f %f [%f]]    # output image voxel sizes\n\
###\n\
[-type-output|-type grid|mosaic]   #\n\
### grid parameters\n\
[-offset %d %d [%d]]  #\n\
[-spacing %d %d [%d]] #\n\
[-value %f]           #\n\
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





char *API_Help_createGrid( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_createGrid( char *program, char *str, int flag )
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



void API_InitParam_createGrid( lineCmdParamCreateGrid *p )
{
    (void)strncpy( p->input_name, "\0", 1 );
    (void)strncpy( p->output_name, "\0", 1 );
    p->output_type = UCHAR;

    (void)strncpy( p->template_name, "\0", 1 );

    p->dim.x = 256;
    p->dim.y = 256;
    p->dim.z = 1;

    p->voxel.x = 1.0;
    p->voxel.y = 1.0;
    p->voxel.z = 1.0;

    p->image = _GRID;

    p->offset.x = 0;
    p->offset.y = 0;
    p->offset.z = 0;

    p->spacing.x = 10;
    p->spacing.y = 10;
    p->spacing.z = 10;

    p->value = 255;

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_createGrid( FILE *theFile, char *program,
                                  lineCmdParamCreateGrid *p, char *str )
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

  fprintf( f, "- input image is " );
  if ( p->input_name != (char*)NULL && p->input_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_name );
  else
    fprintf( f, "'NULL'\n" );

  fprintf( f, "- output image is " );
  if ( p->output_name != (char*)NULL && p->output_name[0] != '\0' )
    fprintf( f, "'%s'\n", p->output_name );
  else
    fprintf( f, "'NULL'\n" );

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

  fprintf( f, "# geometry information for template image\n" );

  fprintf( f, "- p->template_name = " );
  if ( p->template_name[0] != '\0' )
      fprintf( f, "'%s'\n", p->template_name );
  else
      fprintf( f, "NULL\n" );

  BAL_PrintIntegerPoint( f, &(p->dim), "- p->dim = " );
  BAL_PrintDoublePoint( f, &(p->voxel), "- p->voxel = " );

  fprintf( f, "# grid parameter\n" );

  BAL_PrintIntegerPoint( f, &(p->offset), "- p->offset = " );
  BAL_PrintIntegerPoint( f, &(p->spacing), "- p->spacing = " );
  fprintf( f, "- p->value = %f\n", p->value );

  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_createGrid( char *str, lineCmdParamCreateGrid *p )
{
  char *proc = "_API_ParseParam_createGrid";
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

  API_ParseParam_createGrid( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_createGrid( int firstargc, int argc, char *argv[],
                                  lineCmdParamCreateGrid *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
  char text[STRINGLENGTH];
  int status;
  int maxchunks;
  int o=0, s=0, r=0;

  _n_call_parse_ ++;

  /* option line parsing
   */
  for ( i=firstargc; i<argc; i++ ) {

      /* strings beginning with '-'
       */
      if ( argv[i][0] == '-' ) {
          if ( argv[i][1] == '\0' ) {
            if ( inputisread == 0 ) {
              (void)strcpy( p->input_name,  "<" );  /* standart input */
              inputisread = 1;
            }
            else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name,  ">" );  /* standart output */
              outputisread = 1;
            }
            else {
              API_ErrorParse_createGrid( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
            }
          }

          /* template image for image creation
           */
          else if ( strcmp ( argv[i], "-template") == 0
                    || (strcmp ( argv[i], "-t") == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "-dims") == 0 && argv[i][5] == '\0') ) {
            i++;
            if ( i >= argc) API_ErrorParse_createGrid( (char*)NULL, "parsing -template", 0 );
            (void)strcpy( p->template_name, argv[i] );
          }

          else if ( strcmp (argv[i], "-dim" ) == 0 && argv[i][4] == '\0' ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -dim %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.x) );
            if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -dim %d", 0 );
            i ++;
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -dim %d %d", 0 );
            status = sscanf( argv[i], "%d", &(p->dim.y) );
            if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -dim %d %d", 0 );
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
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -voxel %lf", 0 );
            status = sscanf( argv[i], "%lf", &(p->voxel.x) );
            if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -voxel %lf", 0 );
            i ++;
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -voxel %lf %lf", 0 );
            status = sscanf( argv[i], "%lf", &(p->voxel.y) );
            if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -voxel %lf %lf", 0 );
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

          /* image output
           */

          else if ( strcmp (argv[i], "-type-output" ) == 0
                    || (strcmp (argv[i], "-type" ) == 0 && argv[i][5] == '\0') ) {
            i ++;
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -type-output", 0 );
            if ( strcmp (argv[i], "grid" ) == 0 ) {
              p->image = _GRID;
            }
            else if ( strcmp (argv[i], "mosaic" ) == 0 ) {
                p->image = _MOSAIC;
            }
            else {
              if ( _verbose_ )
                fprintf( stderr, "unknown output image tyoe '%s'\n", argv[i] );
              API_ErrorParse_createGrid( (char*)NULL, "parsing -type-output", 0 );
            }
          }

          /* grid parameters
           */

          else if ( strcmp (argv[i], "-offset" ) == 0 ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -offset %d", 0 );
              status = sscanf( argv[i], "%d", &(p->offset.x) );
              if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -offset %d", 0 );
              i ++;
              if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -offset %d %d", 0 );
              status = sscanf( argv[i], "%d", &(p->offset.y) );
              if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -offset %d %d", 0 );
              i ++;
              if ( i < argc) {
                status = sscanf( argv[i], "%d", &(p->offset.z) );
                if ( status <= 0 ) i--;
              }
          }
          else if ( strcmp (argv[i], "-spacing" ) == 0 ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -spacing %d", 0 );
              status = sscanf( argv[i], "%d", &(p->spacing.x) );
              if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -spacing %d", 0 );
              i ++;
              if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -spacing %d %d", 0 );
              status = sscanf( argv[i], "%d", &(p->spacing.y) );
              if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -spacing %d %d", 0 );
              i ++;
              if ( i < argc) {
                status = sscanf( argv[i], "%d", &(p->spacing.z) );
                if ( status <= 0 ) i--;
              }
          }

          else if ( strcmp (argv[i], "-value" ) == 0 ) {
              i ++;
              if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -value", 0 );
              status = sscanf( argv[i], "%f", &(p->value) );
              if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -value", 0 );
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
             if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
               API_ErrorParse_createGrid( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
               API_ErrorParse_createGrid( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
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
             if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -o...\n", 0 );
             status = sscanf( argv[i],"%d",&o );
             if ( status <= 0 ) API_ErrorParse_createGrid( (char*)NULL, "parsing -o...\n", 0 );
          }
          else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
            i += 1;
            if ( i >= argc)    API_ErrorParse_createGrid( (char*)NULL, "parsing -type...\n", 0 );
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
              API_ErrorParse_createGrid( (char*)NULL, "parsing -type...\n", 0 );
            }
          }

          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_createGrid( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_createGrid( (char*)NULL, (char*)NULL, 0);
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
              API_ErrorParse_createGrid( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_createGrid( (char*)NULL, "too long file name ...\n", 0 );
          }
          else if ( inputisread == 0 ) {
              (void)strcpy( p->input_name, argv[i] );
              inputisread = 1;
          }
          else if ( outputisread == 0 ) {
              (void)strcpy( p->output_name, argv[i] );
              outputisread = 1;
          }
          else {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_createGrid( (char*)NULL, "too many file names ...\n", 0 );
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
        API_ErrorParse_createGrid( (char*)NULL, "unable to determine uotput image type ...\n", 0 );
    }
  }

}
