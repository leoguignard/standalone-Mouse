/*************************************************************************
 * malloc.h -
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
 *
 */


#ifndef _vtmalloc_h_
#define _vtmalloc_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

extern void setTraceInVtMalloc( int t );
extern void incrementTraceInVtMalloc( );
extern void setAllocationsInVtMalloc( int a );

extern void vtfree( void *ptr );
extern void *vtmalloc( size_t size, char *var, char *from );

extern void clearVtMalloc( );
void fprintfVtMallocTrace( FILE *f );

#ifdef __cplusplus
}
#endif

#endif
