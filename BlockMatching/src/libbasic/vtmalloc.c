/*************************************************************************
 * malloc.c -
 *
 * Copyright (c) INRIA 2016
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 16 oct 2016 21:00:33 CEST
 *
 *
 * ADDITIONS, CHANGES
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <vtmalloc.h>

static int _verbose_ = 1;


static int _trace_allocations_ = 0;

void setTraceInVtMalloc( int t )
{
  _trace_allocations_ = t;
}

void incrementTraceInVtMalloc( )
{
  _trace_allocations_++;
}

/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/


#define MALLOCSTRINGLENGTH 41
typedef struct _allocation {
  void *ptr;
  size_t size;
  size_t sum_size;
  int in_use;
  char var[MALLOCSTRINGLENGTH];
  char from[MALLOCSTRINGLENGTH];
} _allocation;



static void _initAllocation( _allocation *a )
{
  int i;
  a->ptr = (void*)NULL;
  a->size = 0;
  a->sum_size = 0;
  a->in_use = 0;
  for ( i=0; i<MALLOCSTRINGLENGTH; i++ )
    a->var[i] = '\0';
  for ( i=0; i<MALLOCSTRINGLENGTH; i++ )
    a->from[i] = '0';
}



static void _fprintfAllocation( FILE *f, _allocation *a )
{
    if ( a->in_use == 0 )
      fprintf( f, "x " );
    else
      fprintf( f, "  " );

    fprintf( f, "ptr=%14p, size=%lu, total_size=%lu, var='%s', from='%s'\n",
             a->ptr, a->size,
             a->sum_size,  a->var,
             a->from );
}



typedef struct _allocationList {
   _allocation *data;
   int n_data;
   int n_allocated_data;
   size_t max_size;
   size_t sum_size;
   int n_other_free;
} _allocationList;


static int allocationListIsInitialized = 0;
static _allocationList allocationList;



static void _initAllocationList( _allocationList *l )
{
  l->data = (_allocation*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->max_size = 0;
  l->sum_size = 0;
  l->n_other_free = 0;
}



static void _freeAllocationList( _allocationList *l )
{
  if ( allocationListIsInitialized == 0 ) {
      _initAllocationList( &allocationList );
      allocationListIsInitialized = 1;
  }
  if ( l->data != (_allocation*)NULL )
    free( l->data );
  _initAllocationList( l );
}





static int _allocations_to_be_allocated_ = 10000;



void setAllocationsInVtMalloc( int a )
{
  if ( a > 1 )
    _allocations_to_be_allocated_ = a;
}





static int _addAllocationToList( _allocationList *l, _allocation *a )
{
  char *proc = "_addAllocationToList";
  int i, s =  l->n_allocated_data;
  _allocation *tmp = (_allocation*)NULL;
  _allocation *data;
  _allocation n;

  if ( allocationListIsInitialized == 0 ) {
      _initAllocationList( &allocationList );
      allocationListIsInitialized = 1;
  }

  if ( l->n_data == l->n_allocated_data ) {
    s += _allocations_to_be_allocated_;
    data = (_allocation*)malloc( s * sizeof(_allocation) );
    if ( data == (_allocation*)NULL ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: allocation error", proc );
        fprintf( stderr, "var='data' from='%s' of size=%lu\n", proc, s * sizeof(_allocation) );
        fprintf( stderr, "\t current allocation size is %lu\n", l->sum_size );
      }
      return( -1 );
    }
    for ( i=l->n_allocated_data; i<s; i++ )
      _initAllocation( &(data[i]) );
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_allocation) );
      tmp = l->data;
    }
    l->n_allocated_data = s;
    l->data = data;
    if ( tmp != (_allocation*)NULL ) vtfree( tmp );

    _initAllocation( &n );
    n.ptr = (void*)data;
    n.size = s * sizeof(_allocation);
    n.in_use = 1;
    (void)strncpy( n.var, "l->data", MALLOCSTRINGLENGTH-1 );
    (void)strncpy( n.from, proc, MALLOCSTRINGLENGTH-1 );

    (void)memcpy( &(l->data[l->n_data]), &n, sizeof(_allocation) );
    l->sum_size += n.size;
    l->data[l->n_data].sum_size = l->sum_size;
    if ( l->max_size < l->sum_size ) {
      l->max_size = l->sum_size;
      if ( _trace_allocations_ >= 3 )
        fprintf( stderr, "\t new maximum current allocation size: %lu\n", l->max_size );
    }
    l->n_data ++;

  }

  (void)memcpy( &(l->data[l->n_data]), a, sizeof(_allocation) );
  l->sum_size += a->size;
  l->data[l->n_data].sum_size = l->sum_size;
  if ( l->max_size < l->sum_size ) {
    l->max_size = l->sum_size;
    if ( _trace_allocations_ >= 3 )
      fprintf( stderr, "\t new maximum current allocation size: %lu\n", l->max_size );
  }
  l->n_data ++;

  return( 1 );

}





static void _fprintfAllocationSize( FILE *f, char *str, size_t size )
{
  size_t s = size;
  size_t i;
  if ( s < 1024 ) {
    fprintf( f, "%lu bytes", s );
    return;
  }
  if ( _trace_allocations_ >= 2 ) {
     fprintf( f, "%lu bytes\n", s );
     for ( i=0; i<strlen(str); i++ )
       fprintf( f, " " );
  }

  s = s >> 10;
  if ( s < 1024 ) {
    fprintf( f, "%lu kbytes", s );
    return;
  }
  if ( _trace_allocations_ >= 2 ) {
     fprintf( f, "%lu bytes\n", s );
     for ( i=0; i<strlen(str); i++ )
       fprintf( f, " " );
  }

  s = s >> 10;
  if ( s < 1024 ) {
    fprintf( f, "%lu mbytes", s );
    return;
  }
  if ( _trace_allocations_ >= 2 ) {
     fprintf( f, "%lu bytes\n", s );
     for ( i=0; i<strlen(str); i++ )
       fprintf( f, " " );
  }

  s = s >> 10;
  if ( s < 1024 ) {
    fprintf( f, "%lu gbytes", s );
    return;
  }
  if ( _trace_allocations_ >= 2 ) {
     fprintf( f, "%lu bytes\n", s );
     for ( i=0; i<strlen(str); i++ )
       fprintf( f, " " );
  }

  s = s >> 10;
  fprintf( f, "%lu tbytes", s );
}





static void _fprintfAllocationList( FILE *f, _allocationList *l )
{
  int i;
  if ( allocationListIsInitialized == 0 ) {
      _initAllocationList( &allocationList );
      allocationListIsInitialized = 1;
  }
  if ( _trace_allocations_ >= 2 ) {
    for ( i=0; i<l->n_data; i++ ) {
      fprintf( f, "alloc[%4d]: ", i );
      _fprintfAllocation( f, &(l->data[i]) );
    }
  }
  fprintf( f, "- current allocation size = " );
  _fprintfAllocationSize( f, "- current allocation size = ", l->sum_size );
  fprintf( f, "\n" );
  fprintf( f, "- maximum allocation size = " );
  _fprintfAllocationSize( f, "- maximum allocation size = ", l->max_size );
  fprintf( f, "\n" );

  fprintf( f, "- pointers freed but not traced = %d\n", l->n_other_free );
  fprintf( f, "warning, some allocations may have not been traced\n" );
}





/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/



void vtfree( void *ptr )
{
  int i, f;
  if ( _trace_allocations_ ) {
    for ( i=0, f=0; i<allocationList.n_data && f==0; i++ ) {
      if ( allocationList.data[i].in_use == 0 ) continue;
      if ( allocationList.data[i].ptr != ptr ) continue;
      allocationList.data[i].in_use = 0;
      allocationList.sum_size -= allocationList.data[i].size;
      f = 1;
    }
    if ( f == 0 ) allocationList.n_other_free ++;
  }
  free( ptr );
}



void *vtmalloc( size_t size, char *var, char *from )
{
  char *proc = "vtmalloc";
  _allocation a;
  void *ptr;

  if ( size <= 0 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: allocation size negative or null\n", proc );
      fprintf( stderr, "\t var='%s' from='%s' of size=%lu\n",
               var, from, size );
    }
    return( (void*)NULL );
  }

  ptr = malloc( size );
  if ( ptr == (void*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: allocation failed\n", proc );
      fprintf( stderr, "\t var='%s' from='%s' of size=%lu\n",
               var, from, size );
      if ( _trace_allocations_ )
        fprintf( stderr, "\t current allocation size is %lu\n", allocationList.sum_size );
    }
    if ( _trace_allocations_ >= 2 )
      fprintfVtMallocTrace( stderr );
    return( (void*)NULL );
  }

  if ( _trace_allocations_ ) {
    _initAllocation( &a );
    a.ptr = ptr;
    a.size = size;
    a.in_use = 1;
    if ( var != (char*)NULL && var[0] != '\0' ) {
      (void)strncpy( a.var, var, MALLOCSTRINGLENGTH );
      a.var[MALLOCSTRINGLENGTH-1] = '\0';
    }
    if ( from != (char*)NULL && from[0] != '\0' ) {
      (void)strncpy( a.from, from, MALLOCSTRINGLENGTH );
      a.from[MALLOCSTRINGLENGTH-1] = '\0';
    }
    if ( _addAllocationToList( &allocationList, &a ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to add allocation to list\n", proc );
      }
    }
  }

  return( ptr );
}



void clearVtMalloc( )
{
  _freeAllocationList( &allocationList );
}


void fprintfVtMallocTrace( FILE *f )
{
   _fprintfAllocationList( f, &allocationList );
}

