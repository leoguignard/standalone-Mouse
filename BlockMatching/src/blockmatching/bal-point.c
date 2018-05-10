/*************************************************************************
 * bal-point.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 23 sep 2013 16:54:01 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vtmalloc.h>

#include <bal-point.h>

#define STRLENGTH 512

static int _verbose_ = 1;
static int _debug_ = 0;


void BAL_SetVerboseInBalPoint( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalPoint(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalPoint(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}


void BAL_SetDebugInBalPoint( int v )
{
  _debug_ = v;
}

void BAL_IncrementDebugInBalPoint(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalPoint(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}





static int _size_to_be_allocated_ = 100;



/************************************************************
 *
 * MANAGEMENT (bal_typeFieldPointList)
 *
 ************************************************************/

void BAL_InitTypeFieldPointList( bal_typeFieldPointList *l )
{
  l->data = (bal_typeFieldPoint*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->unit =  VOXEL_UNIT;
  l->vx = 1.0;
  l->vy = 1.0;
  l->vz = 1.0;
}



void BAL_FreeTypeFieldPointList( bal_typeFieldPointList *l )
{
  if ( l->data != (bal_typeFieldPoint*)NULL ) {
    vtfree( l->data );
  } 
  BAL_InitTypeFieldPointList( l );
} 



void BAL_PrintTypeFieldPointList( FILE *f, bal_typeFieldPointList *l )
{
  int i;
  
  for ( i=0; i<l->n_data; i++ )
    BAL_PrintTypeFieldPoint( f, &(l->data[i]), (char*)NULL );
}





static int BAL_AddTypeFieldPointToTypeFieldPointList( bal_typeFieldPointList *l, 
						      bal_typeFieldPoint *c )
{
  char *proc = "BAL_AddTypeFieldPointToTypeFieldPointList";
  int s =  l->n_allocated_data;
  bal_typeFieldPoint *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (bal_typeFieldPoint*)vtmalloc( s * sizeof(bal_typeFieldPoint), "data", proc );
    if ( data == (bal_typeFieldPoint*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(bal_typeFieldPoint) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = *c;
  l->n_data ++;

  return( 1 );
} 





int BAL_ReadTypeFieldPointList( bal_typeFieldPointList *l, char *filename )
{
  char *proc = "BAL_ReadTypeFieldPointList";
  FILE *f;
  int i, j, r;
  float x, y, z;
  bal_typeFieldPoint c;
  char str[STRLENGTH];
  char *s;

  if ( filename != (char*)NULL
       && filename[0] != '\0'
       && filename[0] != '<' ) {
    f = fopen( filename, "r" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
      return( -1 );
    }
  }
  else {
    f = stdin;
  }

  i = 0;
  while ( fgets( str, STRLENGTH, f ) != NULL ) {

    s = str;

    /* skip trailing blanks
     */
    for ( j=0; j<STRLENGTH && (*s == ' ' || *s == '\t'); j++ )
      s++;

    /* skip empty lines
     */
    if ( *s == '\0' || *s == '\n' )
      continue;

    /* skip comments
     */
    if ( *s == '#' || *s == '%' || (s[0] =='/' && s[1] == '/') )
      continue;

    r=sscanf( str, "%f %f %f", &x, &y, &z );

    switch( r ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: error when translating '%s'\n", proc, s );
      continue;
      break;
    case 3 :
      break;
    case 2 :
      z = 0.0;
      break;
    }

    if ( _debug_ )
      fprintf( stderr, "#%04d: read %lf %lf %lf\n", i, x, y, z );

    c.x = x;
    c.y = y;
    c.z = z;

    if ( BAL_AddTypeFieldPointToTypeFieldPointList( l, &c ) != 1 ) {
      if ( f != stdin ) fclose( f );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
                 proc, i, x, y, z );
      return( -1 );
    }
    i ++;
  }

  if ( 0 && _verbose_ )
    fprintf( stderr, "%s: read %d values\n", proc, i );

  if ( f != stdin ) fclose( f );

  return( 1 );
}





/************************************************************
 *
 * MANAGEMENT (bal_doublePointList)
 *
 ************************************************************/

void BAL_InitDoublePointList( bal_doublePointList *l )
{
  l->data = (bal_doublePoint*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->unit = REAL_UNIT;
  l->vx = 1.0;
  l->vy = 1.0;
  l->vz = 1.0;
}



int BAL_AllocDoublePointList( bal_doublePointList *l, int s )
{
  char *proc = "BAL_AllocDoublePointList";
  bal_doublePoint *data;

  if ( s <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative or null number of points\n", proc );
    return( -1 );
  }

  data = (bal_doublePoint*)vtmalloc( s * sizeof(bal_doublePoint), "data", proc );
  if ( data == (bal_doublePoint*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  l->n_data = 0;
  l->n_allocated_data = s;
  l->data = data;

  return( 1 );
}



void BAL_FreeDoublePointList( bal_doublePointList *l )
{
  if ( l->data != (bal_doublePoint*)NULL ) {
    vtfree( l->data );
  }
  BAL_InitDoublePointList( l );
}



void BAL_PrintDoublePointList( FILE *f, bal_doublePointList *l )
{
  int i;

  fprintf( f, "%d points, %d allocated points\n", l->n_data, l->n_allocated_data );

  for ( i=0; i<l->n_data; i++ )
    BAL_PrintDoublePoint( f, &(l->data[i]), (char*)NULL );
}





static int BAL_AddDoublePointToDoublePointList( bal_doublePointList *l,
                  bal_doublePoint *c )
{
  char *proc = "BAL_AddDoublePointToDoublePointList";
  int s =  l->n_allocated_data;
  bal_doublePoint *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (bal_doublePoint*)vtmalloc( s * sizeof(bal_doublePoint), "data", proc );
    if ( data == (bal_doublePoint*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(bal_doublePoint) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = *c;
  l->n_data ++;

  return( 1 );
}





int BAL_ReadDoublePointList( bal_doublePointList *l, char *filename )
{
  char *proc = "BAL_ReadDoublePointList";
  FILE *f;
  int i, j, r;
  float x, y, z;
  bal_doublePoint c;
  char str[STRLENGTH];
  char *s;

  if ( filename != (char*)NULL
       && filename[0] != '\0'
       && filename[0] != '<' ) {
    f = fopen( filename, "r" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
      return( -1 );
    }
  }
  else {
    f = stdin;
  }

  i = 0;
  while ( fgets( str, STRLENGTH, f ) != NULL ) {

    s = str;

    /* skip trailing blanks
     */
    for ( j=0; j<STRLENGTH && (*s == ' ' || *s == '\t'); j++ )
      s++;

    /* skip empty lines
     */
    if ( *s == '\0' || *s == '\n' )
      continue;

    /* skip comments
     */
    if ( *s == '#' || *s == '%' || (s[0] =='/' && s[1] == '/') )
      continue;

    r=sscanf( str, "%f %f %f", &x, &y, &z );

    switch( r ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: error when translating '%s'\n", proc, s );
      continue;
      break;
    case 3 :
      break;
    case 2 :
      z = 0.0;
      break;
    }

    if ( _debug_ )
      fprintf( stderr, "#%04d: read %lf %lf %lf\n", i, x, y, z );

    c.x = x;
    c.y = y;
    c.z = z;

    if ( BAL_AddDoublePointToDoublePointList( l, &c ) != 1 ) {
      if ( f != stdin ) fclose( f );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add measure #%d (=%f %f %f) to list\n",
                 proc, i, x, y, z );
      return( -1 );
    }
    i ++;
  }

  if ( 0 && _verbose_ )
    fprintf( stderr, "%s: read %d values\n", proc, i );

  if ( f != stdin ) fclose( f );

  return( 1 );
}





int BAL_WriteDoublePointList( bal_doublePointList *l, char *filename )
{
  char *proc = "BAL_WriteDoublePointList";
  FILE *f;
  int i;

  if ( filename != (char*)NULL
       && filename[0] != '\0'
       && filename[0] != '>' ) {
    f = fopen( filename, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
      return( -1 );
    }
  }
  else {
    f = stdout;
  }


  for ( i=0; i<l->n_data; i++ ) {
     fprintf( f, "%f %f %f\n", l->data[i].x, l->data[i].y, l->data[i].z );
  }

  if ( f != stdout ) fclose( f );

  return( 1 );
}

