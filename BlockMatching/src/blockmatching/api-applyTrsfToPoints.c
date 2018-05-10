/*************************************************************************
 * api-applyTrsfToPoints.c -
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
 * sed -e "s/applyTrsfToPoints/execuTable/g" \
 *     -e "s/ApplyTrsfToPoints/ExecuTable/g" \
 *     -e "s/applytrsftopoints/executable/g" \
 *     [api-]applyTrsfToPoints.[c,h] > [api-]execTable.[c,h]
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <chunks.h>
#include <vtmalloc.h>

#include <bal-point.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

#include <api-applyTrsfToPoints.h>






static int _verbose_ = 1;
static int _debug_ = 0;


static void _API_ParseParam_applyTrsfToPoints( char *str, lineCmdParamApplyTrsfToPoints *p );



/************************************************************
 *
 * main API
 *
 ************************************************************/

int applyTrsfToPoint( bal_doublePoint thePt,
                      bal_doublePoint* resPt,
                      char *real_transformation_name,
                      int verbose,
                      int debug )
{
  char *proc = "applyTrsfToPoint";
  bal_transformation theTrsf;

  if ( _verbose_ )
      fprintf( stderr, "Warning, '%s' is obsolete\n", proc );

  _verbose_ = verbose;
  _debug_ = debug;

  BAL_InitTransformation( &theTrsf );

  if ( real_transformation_name != NULL &&  real_transformation_name[0] != '\0' ) {
    if ( BAL_ReadTransformation( &theTrsf, real_transformation_name ) != 1 ) {
      if ( _verbose_ )
          fprintf( stderr, "%s: unable to read 'real' transformation '%s'\n", proc, real_transformation_name );
      return( -1 );
    }
  }
  else {
    if ( _verbose_ )
      fprintf( stderr, "%s: no transformation \n", proc );
    return -1;
  }

  if ( BAL_TransformDoublePoint( &thePt, resPt, &theTrsf ) != 1 ) {
    BAL_FreeTransformation( &theTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to transform point (%f,%f,%f)\n", proc, thePt.x, thePt.y, thePt.z );
    return( -1 );
  }

  if ( _debug_ )
    fprintf( stderr, "\t write %lf %lf %lf\n", resPt->x, resPt->y, resPt->z );

  BAL_FreeTransformation( &theTrsf );

  return( 0 );
}





int API_applyTrsfToPoints( bal_doublePointList *thePoints,
                           bal_doublePointList *resPoints,
                           bal_transformation *theTrsf,
                           char *param_str_1,
                           char *param_str_2 )
{
  char *proc = "API_applyTrsfToPoints";
  lineCmdParamApplyTrsfToPoints par;
  int i;



  /* parameter initialization
   */
  API_InitParam_applyTrsfToPoints( &par );

  /* parameter parsing
   */
  if ( param_str_1 != (char*)NULL )
      _API_ParseParam_applyTrsfToPoints( param_str_1, &par );
  if ( param_str_2 != (char*)NULL )
      _API_ParseParam_applyTrsfToPoints( param_str_2, &par );

  if ( par.print_lineCmdParam )
      API_PrintParam_applyTrsfToPoints( stderr, proc, &par, (char*)NULL );

  if ( thePoints->unit != theTrsf->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: points to be transformed and transformation are in different units\n", proc );
    return( -1 );
  }

  if ( thePoints->n_data > resPoints->n_allocated_data ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: reference point list too small, job will be done partially\n", proc );
  }


  /************************************************************
   *
   *  here is the stuff
   *
   ************************************************************/

  if ( theTrsf == (bal_transformation*)NULL ) {
    for (i=0; i<thePoints->n_data && i<resPoints->n_allocated_data; i++ ) {
       resPoints->data[i] = thePoints->data[i];
    }
  }
  else {
    for (i=0; i<thePoints->n_data && i<resPoints->n_allocated_data; i++ ) {
      if ( BAL_TransformDoublePoint( &(thePoints->data[i]),
                                     &(resPoints->data[i]), theTrsf ) != 1 ) {
         if ( _verbose_ )
           fprintf( stderr, "%s: unable to transform point #%d\n", proc, i );
         return( -1 );
      }
    }
  }

  if ( thePoints->n_data > resPoints->n_allocated_data )
    resPoints->n_data = resPoints->n_allocated_data;
  else
    resPoints->n_data = thePoints->n_data;

  resPoints->unit = thePoints->unit;

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
static char *usage = "%s %s\n\
 [-transformation |-trsf %s]\n\
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
 transform a list of points from a file with a given transformation.\n\
 In the file, each line consists in the 3 coordinates 'x y z' of each\n\
 point.\n\
# ...\n\
 -transformation|-trsf %s # transformation to be applied\n\
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





char *API_Help_applyTrsfToPoints( int h )
{
    if ( h == 0 )
        return( usage );
    return( detail );
}





void API_ErrorParse_applyTrsfToPoints( char *program, char *str, int flag )
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



void API_InitParam_applyTrsfToPoints( lineCmdParamApplyTrsfToPoints *p )
{
    (void)strncpy( p->input_name, "\0", 1 );
    (void)strncpy( p->output_name, "\0", 1 );
    (void)strncpy( p->input_real_transformation, "\0", 1 );

    p->print_lineCmdParam = 0;
    p->print_time = 0;
}





void API_PrintParam_applyTrsfToPoints( FILE *theFile, char *program,
                                  lineCmdParamApplyTrsfToPoints *p, char *str )
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

  fprintf( f, "- transformation file name is " );
  if ( p->input_real_transformation != (char*)NULL && p->input_real_transformation[0] != '\0' )
    fprintf( f, "'%s'\n", p->input_real_transformation );
  else
    fprintf( f, "'NULL'\n" );


  fprintf( f, "==================================================\n" );
}





/************************************************************
 *
 * parameters parsing
 *
 ************************************************************/



static void _API_ParseParam_applyTrsfToPoints( char *str, lineCmdParamApplyTrsfToPoints *p )
{
  char *proc = "_API_ParseParam_applyTrsfToPoints";
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

  API_ParseParam_applyTrsfToPoints( 0, argc, argv, p );

  vtfree( argv );
}





static int _n_call_parse_ = 0;

void API_ParseParam_applyTrsfToPoints( int firstargc, int argc, char *argv[],
                                  lineCmdParamApplyTrsfToPoints *p )
{
  int i;
  int inputisread = 0;
  int outputisread = 0;
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
              API_ErrorParse_applyTrsfToPoints( (char*)NULL, "too many file names, parsing '-' ...\n", 0 );
            }
          }

          else if ( strcmp ( argv[i], "-transformation") == 0
                    || (strcmp ( argv[i], "-trsf") == 0 && argv[i][5] == '\0') ) {
                 i++;
                 if ( i >= argc) API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -transformation", 0 );
                 (void)strcpy( p->input_real_transformation, argv[i] );
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
             if ( i >= argc)    API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
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
               API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -parallelism-type ...\n", 0 );
             }
          }

          else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             status = sscanf( argv[i], "%d", &maxchunks );
             if ( status <= 0 ) API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -max-chunks ...\n", 0 );
             if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
          }

          else if ( strcmp ( argv[i], "-omp-scheduling" ) == 0 ||
                   ( strcmp ( argv[i], "-omps" ) == 0 && argv[i][5] == '\0') ) {
             i ++;
             if ( i >= argc)    API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -omp-scheduling, no argument\n", 0 );
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
               API_ErrorParse_applyTrsfToPoints( (char*)NULL, "parsing -omp-scheduling ...\n", 0 );
             }
          }


          /* general parameters
           */
          else if ( (strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0')
                    || (strcmp ( argv[i], "--help" ) == 0 && argv[i][6] == '\0') ) {
             API_ErrorParse_applyTrsfToPoints( (char*)NULL, (char*)NULL, 1);
          }
          else if ( (strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0')
                    || (strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0') ) {
             API_ErrorParse_applyTrsfToPoints( (char*)NULL, (char*)NULL, 0);
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
              API_ErrorParse_applyTrsfToPoints( (char*)NULL, text, 0);
          }
      }

      /* strings beginning with a character different from '-'
       */
      else {
          if ( strlen( argv[i] ) >= STRINGLENGTH ) {
              fprintf( stderr, "... parsing '%s'\n", argv[i] );
              API_ErrorParse_applyTrsfToPoints( (char*)NULL, "too long file name ...\n", 0 );
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
              API_ErrorParse_applyTrsfToPoints( (char*)NULL, "too many file names ...\n", 0 );
          }
      }
  }

  /* if not enough file names
   */
  if ( inputisread == 0 ) {
    (void)strcpy( p->input_name,  "<" );  /* standart input */
    inputisread = 1;
  }
  if ( outputisread == 0 ) {
    (void)strcpy( p->output_name,  ">" );  /* standart output */
    outputisread = 1;
  }

}
