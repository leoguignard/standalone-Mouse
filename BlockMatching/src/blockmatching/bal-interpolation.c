/*************************************************************************
 * bal-interpolation.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun  4 avr 2016 22:09:07 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vtmalloc.h>

#include <bal-transformation-copy.h>
#include <bal-transformation-inversion.h>

#include <bal-interpolation.h>

static int _verbose_ = 1;
static int _debug_ = 0;

void BAL_SetVerboseInBalInterpolation( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalInterpolation(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalInterpolation(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void BAL_SetDebugInBalInterpolation( int d )
{
  _debug_ = d;
}

void BAL_IncrementDebugInBalInterpolation(  )
{
  _debug_ ++;
}

void BAL_DecrementDebugInBalInterpolation(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




/**********************************************************************
 *
 * labels
 *
 **********************************************************************/




static void _initLabelList( typeLabelList *l )
{
  l->data = (int*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeLabelList( typeLabelList *l )
{
  if ( l->data != (int*)NULL )
    vtfree( l->data );
  _initLabelList( l );
}



static int _labels_to_be_allocated_ = 3;

static int _addLabelToLabelList( typeLabelList *l, int n )
{
  char *proc = "_addLabelToLabelList";
  int s =  l->n_allocated_data;
  int *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _labels_to_be_allocated_;
    data = (int*)vtmalloc( s * sizeof(int), "data", proc );
    if ( data == (int*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(int) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = n;
  l->n_data ++;

  return( 1 );
}



static void _fprintfLabelList( FILE *f, typeLabelList *l )
{
  int i;
  fprintf( f, "[" );
  for ( i=0; i<l->n_data; i++ ) {
    fprintf( f, "%d", l->data[i] );
    if ( i < l->n_data-1 )
      fprintf( f, "," );
  }
  fprintf( f, "]" );
}





/**********************************************************************
 *
 * cells
 *
 **********************************************************************/



static void _initCell( typeCell *c )
{
  c->left_corner.x = -1;
  c->left_corner.y = -1;
  c->left_corner.z = -1;

  c->right_corner.x = -1;
  c->right_corner.y = -1;
  c->right_corner.z = -1;

  c->npoints = 0;

  c->new_label = -1;

  _initLabelList( &(c->corresp) );
}



static void _freeCell( typeCell *c )
{
  _freeLabelList( &(c->corresp) );
  _initCell( c );
}



static void _fprintfCell( FILE *f, typeCell *c )
{
  fprintf( f, "(%3d,%3d,%3d) x (%3d,%3d,%3d)",
           c->left_corner.x, c->left_corner.y, c->left_corner.z,
           c->right_corner.x, c->right_corner.y, c->right_corner.z );
  fprintf( f, " - " );
  fprintf( f, " # = %5d ", c->npoints );
  fprintf( f, " - " );
  _fprintfLabelList( f, &(c->corresp) );
  fprintf( f, "\n" );
}



static void _initCellList( typeCellList *l )
{
  l->data = (typeCell*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeCellList( typeCellList *l )
{
  int i;
  if ( l->data != (typeCell*)NULL ) {
    for ( i=0; i<l->n_data; i++ )
      _freeCell( &(l->data[i]) );
    vtfree( l->data );
  }
  _initCellList( l );
}



static int _addLabelToCellList( typeCellList *l, int n )
{
  char *proc = "_addLabelToCellList";
  int i;
  typeCell *data;

  if ( n < 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with negative labels\n", proc );
    return( -1 );
  }

  if ( n < l->n_allocated_data )
    return( 1 );

  data = (typeCell *)vtmalloc( (n+1) * sizeof(typeCell), "data", proc );
  if ( l->n_allocated_data > 0 ) {
    (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeCell) );
    vtfree( l->data );
  }
  for ( i=l->n_allocated_data; i<=n; i++ )
    _initCell( &(data[i]) );

  l->n_allocated_data = n+1;
  l->n_data = n+1;
  l->data = data;

  return( 1 );
}



static void _fprintfCellList( FILE *f, typeCellList *l, char *d )
{
  int i;

  if ( d != (char*)NULL ) {
    fprintf( f, "--- '%s' ---\n", d );
  }
  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i].left_corner.x >=0 || l->data[i].left_corner.y >= 0
         || l->data[i].left_corner.z >= 0
         || l->data[i].corresp.n_data > 0 ) {
      if ( d != (char*)NULL )
        fprintf( f, "%s", d );
      else
        fprintf( f, "cell" );
      fprintf( f, "[%4d] = ", i );
      _fprintfCell( f, &(l->data[i]) );
    }
  }
  fprintf( f, "\n" );
}



int BAL_FillCellList( bal_image *image, typeCellList *l )
{
  char *proc = "BAL_FillCellList";
  int min, max;
  size_t i, v;
  int x, y, z;
  typeCell *c;

  v = image->vdim * image->nrows * image->ncols * image->nplanes;

#define _CELLLISTMINMAX( TYPE ) {     \
  TYPE *buf = (TYPE*)image->data;     \
  min = max = buf[0];                 \
  for ( i=1; i<v; i++ ) {             \
    if ( min > buf[i] ) min = buf[i]; \
    if ( max < buf[i] ) max = buf[i]; \
  }                                   \
}
  switch( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _CELLLISTMINMAX( u8 );
    break;
  case SSHORT :
    _CELLLISTMINMAX( s16 );
    break;
  case USHORT :
    _CELLLISTMINMAX( u16 );
    break;
  }

  if ( min < 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: minimal value of labels is negative\n", proc );
    return( -1 );
  }

  if ( _debug_ )
    fprintf( stderr, "%s: label values are in [%d,%d]\n", proc, min, max );

  l->data = (typeCell*)vtmalloc( (max+1) * sizeof(typeCell), "l->data", proc );
  if ( l->data == (typeCell*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  l->n_data = l->n_allocated_data = max+1;

  for ( i=0; i<=(unsigned int)max; i++ )
    _initCell( &(l->data[i]) );

#define _CELLLISTFILL( TYPE ) {                           \
  TYPE ***buf = (TYPE***)image->array;                    \
  for ( z=0; z<(int)image->nplanes; z++ )                 \
  for ( y=0; y<(int)image->nrows; y++ )                   \
  for ( x=0; x<(int)image->ncols; x++ ) {                 \
    c = &(l->data[buf[z][y][x]]);                         \
    c->npoints ++;                                        \
    if ( c->left_corner.x < 0 || c->left_corner.y < 0 || c->left_corner.z < 0 ) { \
      c->left_corner.x = c->right_corner.x = x;           \
      c->left_corner.y = c->right_corner.y = y;           \
      c->left_corner.z = c->right_corner.z = z;           \
    }                                                     \
    else {                                                \
      if ( c->left_corner.x > x ) c->left_corner.x = x;   \
      if ( c->left_corner.y > y ) c->left_corner.y = y;   \
      if ( c->left_corner.z > z ) c->left_corner.z = z;   \
      if ( c->right_corner.x < x ) c->right_corner.x = x; \
      if ( c->right_corner.y < y ) c->right_corner.y = y; \
      if ( c->right_corner.z < z ) c->right_corner.z = z; \
    }                                                     \
  }                                                       \
}

  switch( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _CELLLISTFILL( u8 );
    break;
  case SSHORT :
    _CELLLISTFILL( s16 );
    break;
  case USHORT :
    _CELLLISTFILL( u16 );
    break;
  }

  return( 1 );
}





static int BAL_RelabelCellImage( bal_image *image, typeCellList *l )
{
  char *proc = "BAL_RelabelCellImage";
  size_t i, v;

  v = image->vdim * image->nrows * image->ncols * image->nplanes;

#define _RELABELCELLIMAGE( TYPE ) {     \
  TYPE *buf = (TYPE*)image->data;       \
  for ( i=0; i<v; i++ ) {               \
    buf[i] = l->data[buf[i]].new_label; \
  }                                     \
}

  switch( image->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _RELABELCELLIMAGE( u8 );
    break;
  case SSHORT :
    _RELABELCELLIMAGE( s16 );
    break;
  case USHORT :
    _RELABELCELLIMAGE( u16 );
    break;
  }

  return( 1 );
}





/**********************************************************************
 *
 * correspondences
 *
 **********************************************************************/



static void _initCorrespondence( typeCorrespondence *c )
{
  c->label0 = &(c->allocatedLabel0);
  c->label1 = &(c->allocatedLabel1);
  _initLabelList( &(c->allocatedLabel0) );
  _initLabelList( &(c->allocatedLabel1) );
}



static void _freeCorrespondence( typeCorrespondence *c )
{
  _freeLabelList( &(c->allocatedLabel0) );
  _freeLabelList( &(c->allocatedLabel1) );
  _initCorrespondence( c );
}



static void _printCorrespondence( FILE *f, typeCorrespondence *c )
{
  int i;

  if ( c->label0->n_data <= 0 && c->label1->n_data <= 0 )
      return;

  if ( c->label0->n_data > 0 ) {
    for ( i=0; i<c->label0->n_data; i++ )
      fprintf( f, "%d ", c->label0->data[i] );
  }
  else {
    fprintf( f, " " );
  }

  fprintf( f, "-" );

  if ( c->label1->n_data > 0 ) {
    for ( i=0; i<c->label1->n_data; i++ )
      fprintf( f, " %d", c->label1->data[i] );
    /* keep a trailing space to be strictly
     * equivalent to writing procedure of Gael
     */
    fprintf( f, " " );
  }
  else {
    fprintf( f, " " );
  }

  fprintf( f, "\n" );
}





/**********************************************************************
 *
 * correspondence list
 *
 **********************************************************************/



static int _correspondences_to_be_allocated_ = 10;

static int _addCorrespondenceToCorrespondenceList( typeCorrespondenceList *l,
                                                   typeCorrespondence *c)
{
  char *proc = "_addCorrespondenceToCorrespondenceList";
  int i, s = l->n_allocated_data;
  typeCorrespondence *data;
  typeCorrespondence *t;

  if ( l->n_data == l->n_allocated_data ) {

    s += _correspondences_to_be_allocated_;
    data = (typeCorrespondence*)vtmalloc( s * sizeof(typeCorrespondence),
                                          "data", proc );
    if ( data == (typeCorrespondence*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error (array)\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typeCorrespondence) );
      vtfree( l->data );
      /* do not forget to set again the pointers !
       */
      for ( i=0; i<l->n_allocated_data; i++ ) {
        data[i].label0 = &(data[i].allocatedLabel0);
        data[i].label1 = &(data[i].allocatedLabel1);
      }
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  t = &(l->data[l->n_data]);
  _initCorrespondence( t );

  if ( c->allocatedLabel0.n_data > 0 ) {
      t->allocatedLabel0.data = (int*)vtmalloc( c->allocatedLabel0.n_data * sizeof(int),
                                                "t->allocatedLabel0.data", proc );
      if ( t->allocatedLabel0.data == (int*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: allocation error (label0) \n", proc );
          return( -1 );
      }
      (void)memcpy( t->allocatedLabel0.data, c->allocatedLabel0.data, c->allocatedLabel0.n_data *sizeof(int) );
      t->allocatedLabel0.n_data = c->allocatedLabel0.n_data;
      t->allocatedLabel0.n_allocated_data = c->allocatedLabel0.n_data;
  }

  if ( c->allocatedLabel1.n_data > 0 ) {
      t->allocatedLabel1.data = (int*)vtmalloc( c->allocatedLabel1.n_data * sizeof(int),
                                                "t->allocatedLabel1.data", proc );
      if ( t->allocatedLabel1.data == (int*)NULL ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: allocation error (label1) \n", proc );
          return( -1 );
      }
      (void)memcpy( t->allocatedLabel1.data, c->allocatedLabel1.data, c->allocatedLabel1.n_data *sizeof(int) );
      t->allocatedLabel1.n_data = c->allocatedLabel1.n_data;
      t->allocatedLabel1.n_allocated_data = c->allocatedLabel1.n_data;
  }

  l->n_data ++;

  return( 1 );
}






static void _initCorrespondenceList( typeCorrespondenceList *l )
{
  l->data = (typeCorrespondence*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeCorrespondenceList( typeCorrespondenceList *l )
{
  int c;
  if ( l->data != (typeCorrespondence*)NULL ) {
    for ( c=0; c<l->n_data; c++ ) {
      _freeCorrespondence( &(l->data[c]) );
    }
    vtfree( l->data );
  }
}



static void _maxLabelsCorrespondenceList( typeCorrespondenceList *l,
                                          int *max0, int *max1 )
{
  int m0, m1;
  int i, j;

  m0 = m1 = -1;

  for ( i=0; i<l->n_data; i++ ) {
    if ( l->data[i].label0->n_data > 0 ) {
      for ( j=0; j<l->data[i].label0->n_data; j++ ) {
        if ( m0 == -1 ) {
          m0 = l->data[i].label0->data[j];
        }
        else {
          if ( m0 < l->data[i].label0->data[j] )
              m0 = l->data[i].label0->data[j];
        }
      }
    }
    if ( l->data[i].label1->n_data > 0 ) {
      for ( j=0; j<l->data[i].label1->n_data; j++ ) {
         if ( m1 == -1 ) {
           m1 = l->data[i].label1->data[j];
         }
         else {
           if ( m1 < l->data[i].label1->data[j] )
             m1 = l->data[i].label1->data[j];
         }
      }
    }
  }
  *max0 = m0;
  *max1 = m1;
}





/**********************************************************************
 *
 * correspondence list / I/O procedures
 *
 **********************************************************************/



static void BAL_FprintfCorrespondenceList( FILE *f, typeCorrespondenceList *l )
{
  int c;
  for ( c=0; c<l->n_data; c++ ) {
      _printCorrespondence( f, &(l->data[c]) );
  }
}



static int BAL_FscanfCorrespondenceList( FILE *f, typeCorrespondenceList *l )
{
  char *proc = "BAL_FscanfCorrespondenceList";
  typeCorrespondence c;
  typeLabelList *label;
  char line[512], *r;
  int n, stop;
  int ndata = l->n_data;

  _initCorrespondence( &c );

  while ( fgets( line, 512, f ) != (char*)NULL ) {

      c.label0->n_data = 0;
      c.label1->n_data = 0;
      label = c.label0;

      r = line;

      /* comments
       */
      if ( *r == '#' || *r == '%' || (r[0] =='/' && r[0] =='/') )
          continue;

      for ( stop=0; stop==0; ) {
        switch( *r ) {
        default :
            _freeCorrespondence( &c );
            if ( _verbose_ )
                fprintf( stderr, "%s: unexpected reading, line '%s'\n", proc, r );
            return( -1 );
            break;
        case ' ' :
        case '\t' :
            r++;
            break;
        case '\n' :
        case '\0' :
            stop = 1;
            break;
        case '1' :
        case '2' :
        case '3' :
        case '4' :
        case '5' :
        case '6' :
        case '7' :
        case '8' :
        case '9' :
            if ( sscanf( r, "%d", &n ) != 1 ) {
                _freeCorrespondence( &c );
                if ( _verbose_ )
                    fprintf( stderr, "%s: reading of '%s' failed\n", proc, r );
                return( -1 );
            }
            if ( _addLabelToLabelList( label, n ) != 1 ) {
                _freeCorrespondence( &c );
                if ( _verbose_ )
                    fprintf( stderr, "%s: can not add label %d from line '%s'\n", proc, n, line );
                return( -1 );
            }
            while ( *r >= '0' && *r <= '9' )
                r++;
            break;
        case '-' :
            if ( label == c.label0 ) {
                label = c.label1;
            }
            else {
                _freeCorrespondence( &c );
                if ( _verbose_ )
                    fprintf( stderr, "%s: '-' found twice in 'line '%s'\n", proc, line );
                return( -1 );
            }
            r++;
            break;
        }
      }


      if ( c.label0->n_data > 0 || c.label1->n_data > 0 ) {
          if ( _addCorrespondenceToCorrespondenceList( l, &c ) != 1 ) {
              _freeCorrespondence( &c );
              if ( _verbose_ )
                  fprintf( stderr, "%s: can not add correspondence to list\n", proc );
              return( -1 );
          }
      }

  }

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "%s: has read %d correspondences\n", proc, l->n_data - ndata );
  }

  _freeCorrespondence( &c );
  return( 1 );
}



int BAL_ReadCorrespondenceList( typeCorrespondenceList *l, char *name )
{
  char *proc = "BAL_ReadCorrespondenceList";
  FILE *f;

  f = fopen( name, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
    return( -1 );
  }

  if ( BAL_FscanfCorrespondenceList( f, l ) != 1 ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading correspondences\n", proc );
    return( -1 );
  }

  fclose( f );
  return( 1 );
}



int BAL_WriteCorrespondenceList( typeCorrespondenceList *l, char *name )
{
  char *proc = "BAL_WriteCorrespondenceList";
  FILE *f;

  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
    return( -1 );
  }

  BAL_FprintfCorrespondenceList( f, l );

  fclose( f );
  return( 1 );
}










/**********************************************************************
 *
 * cell correspondence
 *
 **********************************************************************/



void BAL_InitCellCorrespondence( typeCellCorrespondence *c )
{
  _initCorrespondenceList( &(c->cliques) );
  c->list0 = &(c->allocatedList0);
  c->list1 = &(c->allocatedList1);
  _initCellList( &(c->allocatedList0) );
  _initCellList( &(c->allocatedList1) );
}



void BAL_FreeCellCorrespondence( typeCellCorrespondence *c )
{
  _freeCorrespondenceList( &(c->cliques) );
  _freeCellList( &(c->allocatedList0) );
  _freeCellList( &(c->allocatedList1) );
  BAL_InitCellCorrespondence( c );
}





/**********************************************************************
 *
 *
 *
 **********************************************************************/



static void _invertCellCorrespondence( typeCellCorrespondence *cc )
{
  int i;
  typeCellList *cl;
  typeLabelList *ll;

  for ( i=0; i<cc->cliques.n_data; i++ ) {
    ll = cc->cliques.data[i].label0;
    cc->cliques.data[i].label0 = cc->cliques.data[i].label1;
    cc->cliques.data[i].label1 = ll;
  }

  cl = cc->list0;
  cc->list0 = cc->list1;
  cc->list1 = cl;
}



static int _background_label_ = 1;

int BAL_CrossFillCellCorrespondence( typeCellCorrespondence *cc )
{
  char *proc = "BAL_CrossFillCellCorrespondence";
  int max0, max1;

  typeCorrespondenceList *cl;

  int i0, i1, k;
  int l0, l1;

  /* labels max
   */
  _maxLabelsCorrespondenceList( &(cc->cliques), &max0, &max1 );

  if ( 0 )
    fprintf( stderr, "%s: max A = %d - max B = %d\n", proc, max0, max1 );

  if ( max0 >= cc->list0->n_allocated_data ) {
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: enlarge list #0 from %d to %d\n",
               proc, cc->list0->n_allocated_data, max0+1 );
    }
    if ( _addLabelToCellList( cc->list0, max0 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when adding label %d to list #0\n", proc, max0 );
      return( -1 );
    }
  }

  if ( max1 >= cc->list1->n_allocated_data ) {
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: enlarge list #1 from %d to %d\n",
               proc, cc->list1->n_allocated_data, max1+1 );
    }
    if ( _addLabelToCellList( cc->list1, max1 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when adding label %d to list #1\n", proc, max1 );
      return( -1 );
    }
  }

  /* fill list
   */
  cl = &(cc->cliques);

  for ( k=0; k<cl->n_data; k++ ) {

    /* right part of the clique is not empty
     */
    if ( cl->data[k].label0->n_data > 0 ) {
      for ( i0=0; i0<cl->data[k].label0->n_data; i0++ ) {
        l0 = cl->data[k].label0->data[i0];
        if ( cl->data[k].label1->n_data > 0 ) {
          for ( i1=0; i1<cl->data[k].label1->n_data; i1++ ) {
            l1 = cl->data[k].label1->data[i1];
            if ( _addLabelToLabelList( &(cc->list0->data[l0].corresp), l1 ) != 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: error when adding %d to list0[%d]\n", proc, l1, l0 );
              return( -1 );
            }
            if ( _addLabelToLabelList( &(cc->list1->data[l1].corresp), l0 ) != 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: error when adding %d to list1[%d]\n", proc, l0, l1 );
              return( -1 );
            }
          }
        }
        else {
          if ( _addLabelToLabelList( &(cc->list0->data[l0].corresp), _background_label_ ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when adding %d to list0[%d]\n", proc, _background_label_, l0 );
            return( -1 );
          }
          if ( _addLabelToLabelList( &(cc->list1->data[_background_label_].corresp), l0 ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when adding %d to list1[%d]\n", proc, l0, _background_label_ );
            return( -1 );
          }
        }
      }
    }

    /* right part of the clique is empty
     */
    else {
      if ( cl->data[k].label1->n_data > 0 ) {
        for ( i1=0; i1<cl->data[k].label1->n_data; i1++ ) {
          l1 = cl->data[k].label1->data[i1];
          if ( _addLabelToLabelList( &(cc->list0->data[_background_label_].corresp), l1 ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when adding %d to list0[%d]\n", proc, l1, _background_label_ );
            return( -1 );
          }
          if ( _addLabelToLabelList( &(cc->list1->data[l1].corresp), _background_label_ ) != 1 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: error when adding %d to list1[%d]\n", proc, _background_label_, l1 );
            return( -1 );
          }
        }
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: weird, empty correspondences\n", proc );
        return( -1 );
      }
    }
  }


  /* add background to background correspondences
   */
  if ( _addLabelToLabelList( &(cc->list0->data[_background_label_].corresp),
                             _background_label_ ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when adding %d (background) to list0[%d]\n", proc, _background_label_, _background_label_ );
    return( -1 );
  }
  if ( _addLabelToLabelList( &(cc->list1->data[_background_label_].corresp),
                             _background_label_ ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when adding %d (background) to list1[%d]\n", proc, _background_label_, _background_label_ );
    return( -1 );
  }


  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "==========   correspondences   ==========\n" );
    BAL_FprintfCorrespondenceList( stderr, &(cc->cliques) );
    _fprintfCellList( stderr, cc->list0, "cell#0" );
    _fprintfCellList( stderr, cc->list1, "cell#1" );
    fprintf( stderr, "=========================================\n" );
    fprintf( stderr, "\n" );
  }

  return( 1 );
}


void _CheckCellList( FILE*f, typeCellList *l,
                     char *desc )
{
  int i, j, len, dlen;
  int n00 = 0;
  int n10 = 0;
  int n01 = 0;
  int n11 = 0;
  typeCell *c;

  for ( i=1; i<l->n_data; i++ ) {
    c = &(l->data[i]);
    if ( c->left_corner.x >= 0 || c->left_corner.y >= 0 || c->left_corner.z >= 0 ) {
      if ( c->corresp.n_data > 0 )
        n11++;
      else
        n10++;
    }
    else {
      if ( c->corresp.n_data > 0 )
        n01++;
      else
        n00++;
    }
  }

  if ( desc != (char*)NULL )
    fprintf( f, "--- %s ---", desc );
  else
    fprintf( f, "--- list ---" );

  fprintf( stderr, " max label = %d\n", l->n_data-1 );

  if ( n00 > 0 )
    fprintf( f, "    labels neither in image nor in cliques: %5d/%5d\n",
             n00, l->n_data-1 );
  if ( n11 > 0 )
    fprintf( f, "    labels both in image and in cliques: %5d/%5d\n",
             n11, l->n_data-1 );
  if ( n01 > 0 )
    fprintf( f, "    labels not in image but in cliques: %5d/%5d\n",
             n01, l->n_data-1 );
  if ( n10 > 0 )
    fprintf( f, "    labels in image but not in cliques: %5d/%5d\n",
             n10, l->n_data-1 );

  if ( n10 > 0 ) {
    fprintf( f, "    [" );
    for ( len=5,j=0, i=1; i<l->n_data; i++ ) {
     c = &(l->data[i]);
     if ( (c->left_corner.x >= 0 || c->left_corner.y >= 0
           || c->left_corner.z >= 0) && c->corresp.n_data <= 0 ) {
       if ( j>0 ) {
         fprintf( f, "," );
         len++;
       }
       if ( i< 10 ) dlen = 1;
       else if ( i < 100 ) dlen = 2;
       else if ( i < 1000 ) dlen = 3;
       else if ( i < 10000 ) dlen = 4;
       else dlen = 5;
       if ( len + dlen > 80 ) {
         fprintf( f, "\n     " );
         len = 5;
       }
       fprintf( f, "%d", i );
       j++;
       len += dlen;
     }
    }
    fprintf( stderr, "]\n" );
  }
}



void BAL_CheckCellCorrespondence( FILE *f, typeCellCorrespondence *cc )
{
  _CheckCellList( f, cc->list0, "list #0" );
  _CheckCellList( f, cc->list1, "list #1" );
}





/**********************************************************************
 *
 * cell relabeling
 *
 **********************************************************************/


int BAL_RelabelCells( bal_image *image0, bal_image *image1,
                     typeCellCorrespondence *cc )
{
  char *proc = "BAL_RelabelCells";
  int i, j, l;
  int iteration, changes;
  int n0, n1;
  int newlabel0 = 1;
  int newlabel1 = 1;
  typeLabelList *label;


  /* initialization
   */
  for ( i=0; i<cc->list0->n_data; i++ )
    cc->list0->data[i].new_label = -1;
  for ( i=0; i<cc->list1->n_data; i++ )
    cc->list1->data[i].new_label = -1;

  cc->list0->data[0].new_label = 0;
  cc->list0->data[_background_label_].new_label = _background_label_;
  cc->list1->data[0].new_label = 0;
  cc->list1->data[_background_label_].new_label = _background_label_;

  /* changes the labels in cell lists
   */
  for ( changes=1, iteration=0; changes != 0; iteration ++ ) {
    switch( iteration ) {
    case 0 : n0 = 1; n1 = 1; break;
    case 1 : n0 = 1; n1 = 2; break;
    case 2 : n0 = 1; n1 = -1; break;
    case 3 : n0 = 2; n1 = 1; break;
    default :
      n0 = -1; n1 = -1; break;
    }

    for ( changes=0, i=0; i<cc->cliques.n_data; i++ ) {
      if ( n0 > 0 && cc->cliques.data[i].label0->n_data != n0 )
        continue;
      if ( n1 > 0 && cc->cliques.data[i].label1->n_data != n1 )
        continue;

      /* we've got a clique #n0 - #n1
       */
      if ( newlabel0 == _background_label_ ) newlabel0 ++;
      if ( newlabel1 == _background_label_ ) newlabel1 ++;

      label = cc->cliques.data[i].label0;
      for ( j=0; j<label->n_data; j++ ) {
        l = label->data[j];
        if ( cc->list0->data[l].new_label != -1 ) continue;
        cc->list0->data[l].new_label = newlabel0 ++;
        changes ++;
      }
      label = cc->cliques.data[i].label1;
      for ( j=0; j<label->n_data; j++ ) {
        l = label->data[j];
        if ( cc->list1->data[l].new_label != -1 ) continue;
        cc->list1->data[l].new_label = newlabel1 ++;
        changes ++;
      }
    }

  }

  for ( changes=0, i=0; i<cc->list0->n_data; i++ ) {
    if ( cc->list0->data[i].new_label == -1 ) {
      cc->list0->data[i].new_label = newlabel0 ++;
      changes ++;
    }
  }

  for ( changes=0, i=0; i<cc->list1->n_data; i++ ) {
    if ( cc->list1->data[i].new_label == -1 ) {
      cc->list1->data[i].new_label = newlabel1 ++;
      changes ++;
    }
  }


  /* change the labels in images
   */
  if ( BAL_RelabelCellImage( image0, cc->list0 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to relabel image #0\n", proc );
    return( -1 );
  }
  if ( BAL_RelabelCellImage( image1, cc->list1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to relabel image #1\n", proc );
    return( -1 );
  }

  /* change the labels in correspondences
   */
  for ( i=0; i<cc->cliques.n_data; i++ ) {
    label = cc->cliques.data[i].label0;
    for ( j=0; j<label->n_data; j++ ) {
      l = label->data[j];
      label->data[j] = cc->list0->data[l].new_label;
    }
    label = cc->cliques.data[i].label1;
    for ( j=0; j<label->n_data; j++ ) {
      l = label->data[j];
      label->data[j] = cc->list1->data[l].new_label;
    }
  }


  /* recompute cell lists
   */
  _freeCellList( &(cc->allocatedList0) );
  _freeCellList( &(cc->allocatedList1) );

  if ( BAL_FillCellList( image0, cc->list0 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build cell list from image #0\n", proc );
    return( -1 );
  }

  if ( BAL_FillCellList( image1, cc->list1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build cell list from image #1\n", proc );
    return( -1 );
  }

  if ( BAL_CrossFillCellCorrespondence( cc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build cell correspondences\n", proc );
    return( -1 );
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "==========   correspondences   ==========\n" );
    BAL_FprintfCorrespondenceList( stderr, &(cc->cliques) );
    _fprintfCellList( stderr, cc->list0, "cell#0" );
    _fprintfCellList( stderr, cc->list1, "cell#1" );
    fprintf( stderr, "=========================================\n" );
    fprintf( stderr, "\n" );
  }

  return( 1 );
}






/**********************************************************************
 *
 * transformation related stuff
 *
 **********************************************************************/



static int BAL_MultiplyTrsf( bal_transformation *trsf, bal_transformation *res, float coef )
{
  char *proc= "BAL_MultiplyTrsf";
  size_t i, j, k;

  float ***arrayIn=NULL, ***arrayOut=NULL;
  double *m=NULL, *M=NULL;

  switch (trsf->type) {
  default :
      if ( _verbose_ )
          fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
      return( -1 );
  case VECTORFIELD_3D :
      arrayIn=(float***)trsf->vz.array;
      arrayOut=(float***)res->vz.array;
      for (k = 0 ; k < trsf->vz.nplanes ; k++ )
      for (j = 0 ; j < trsf->vz.nrows ; j++ )
      for (i = 0 ; i < trsf->vz.ncols ; i++ )
          arrayOut[k][j][i]=coef*arrayIn[k][j][i];
  case VECTORFIELD_2D:
      arrayIn=(float***)trsf->vx.array;
      arrayOut=(float***)res->vx.array;
      for (k = 0 ; k < trsf->vx.nplanes ; k++ )
      for (j = 0 ; j < trsf->vx.nrows ; j++ )
      for (i = 0 ; i < trsf->vx.ncols ; i++ )
          arrayOut[k][j][i]=coef*arrayIn[k][j][i];
      arrayIn=(float***)trsf->vy.array;
      arrayOut=(float***)res->vy.array;
      for (k = 0 ; k < trsf->vy.nplanes ; k++ )
      for (j = 0 ; j < trsf->vy.nrows ; j++ )
      for (i = 0 ; i < trsf->vy.ncols ; i++ )
          arrayOut[k][j][i]=coef*arrayIn[k][j][i];
      break;
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
      m=(double*)trsf->mat.m;
      M=(double*)res->mat.m;
      for (i = 0 ; i < (size_t)(trsf->mat.l * trsf->mat.c - 1); i++ )
          M[i] = coef * m[i];
      M[0] += 1-coef;
      M[5] += 1-coef;
      M[10] += 1-coef;
      break;
  }
  return( 1 );
}





typedef enum enumIntermediaryTrsfDestination {
  _TOWARDS_0,
  _TOWARDS_1
} enumIntermediaryTrsfDestination;


/* the input transformation is assumed to be from 1 to 0,
 * ie {0 <- 1},
 * if flag = _TOWARDS_1, we compute T_{1<-t}
 * if flag = _TOWARDS_0, we compute T_{0<-t}
 */
static int BAL_ComputeTrsfAtIntermediaryTime( bal_transformation *theTrsf,
                                       bal_transformation *resTrsf,
                                       float t,
                                       enumIntermediaryTrsfDestination towards )
{
  char *proc = "BAL_ComputeTrsfAtIntermediaryTime";
  bal_transformation tmpTrsf;

  BAL_InitTransformation( &tmpTrsf );

  switch( theTrsf->type ) {
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    }
    return( -1 );
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    if ( BAL_AllocTransformation( &tmpTrsf, theTrsf->type, (bal_image*)NULL ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to allocate linear transformation\n", proc );
      }
      return( -1 );
    }
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    if ( BAL_AllocTransformation( &tmpTrsf, theTrsf->type, &(theTrsf->vx) ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to allocate non-linear transformation\n", proc );
      }
      return( -1 );
    }
    break;
  }

  /* Multiply Transformation : T(t<-1)
   */
  if ( BAL_MultiplyTrsf( theTrsf, &tmpTrsf, 1.0-t ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to multiply transformation\n", proc );
    }
    return( -1 );
  }

  /* Invert Transformation : T(1<-t)
   */
  if ( BAL_InverseTransformation( &tmpTrsf, resTrsf ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert transformation\n", proc );
    return( -1 );
  }

  if ( towards == _TOWARDS_0 ) {
    /* Multiply Transformation : T(0<-t)
     */
    if ( BAL_MultiplyTrsf( resTrsf, resTrsf, -t/(1.0-t) ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to multiply transformation\n", proc );
      }
      return( -1 );
    }
  }

  BAL_FreeTransformation( &tmpTrsf );

  return( 1 );
}





/**********************************************************************
 *
 * image related stuff
 *
 **********************************************************************/


static int BAL_ImageLinearCombination( float c0, bal_image *I0,
                                float c1, bal_image *I1,
                                bal_image *imres )
{
  char *proc = "BAL_ImageLinearCombination";
  size_t i, v;
  float val;

  v = I0->ncols * I0->nrows * I0->nplanes * I0->vdim;

#define _LINEAR_RES_POSINT( TYPE, MAX ) {             \
  TYPE *resbuf = (TYPE*)imres->data;                  \
  for ( i=0; i<v; i++, resbuf++, buf0++, buf1++ ) {   \
     val = c0 * (float)(*buf0) + c1 * (float)(*buf1); \
     if ( val < 0.0 )                                 \
       *resbuf = (TYPE)0;                             \
     else if ( val > MAX )                            \
       *resbuf = (TYPE)MAX;                           \
     else                                             \
       *resbuf = (TYPE)(val + 0.5);                   \
  }                                                   \
}

#define _LINEAR_RES_NEGINT( TYPE, MIN, MAX ) {        \
  TYPE *resbuf = (TYPE*)imres->data;                  \
  for ( i=0; i<v; i++, resbuf++, buf0++, buf1++ ) {   \
     val = c0 * (float)(*buf0) + c1 * (float)(*buf1); \
     if ( val < MIN )                                 \
       *resbuf = (TYPE)MIN;                           \
     else if ( val < 0.0 )                            \
       *resbuf = (TYPE)(val - 0.5);                   \
     else if ( val > MAX )                            \
       *resbuf = (TYPE)MAX;                           \
     else                                             \
       *resbuf = (TYPE)(val + 0.5);                   \
  }                                                   \
}

#define _LINEAR_RES_OTHER( TYPE ) {                   \
  TYPE *resbuf = (TYPE*)imres->data;                  \
  for ( i=0; i<v; i++, resbuf++, buf0++, buf1++ ) {   \
     val = c0 * (float)(*buf0) + c1 * (float)(*buf1); \
       *resbuf = (TYPE)(val);                         \
  }                                                   \
}

#define _LINEAR_SWITCH_RES {                  \
  switch( imres->type ) {                     \
  default :                                   \
    if ( _verbose_ )                          \
      fprintf( stderr, "%s: such output image type not handled yet\n", proc ); \
    return( -1 );                             \
  case SCHAR :                                \
    _LINEAR_RES_NEGINT( s8, -128, 127 );      \
    break;                                    \
  case UCHAR :                                \
    _LINEAR_RES_POSINT( u8, 255 );            \
    break;                                    \
  case SSHORT :                               \
    _LINEAR_RES_NEGINT( s16, -32768, 32767 ); \
    break;                                    \
  case USHORT :                               \
    _LINEAR_RES_POSINT( u16, 65535 );         \
    break;                                    \
  case FLOAT :                                \
    _LINEAR_RES_OTHER( r32 );                 \
    break;                                    \
  }                                           \
}

#define _LINEAR_SCNDINPUT( TYPE ) { \
  TYPE *buf1 = (TYPE*)I1->data;     \
  _LINEAR_SWITCH_RES;               \
}

#define _LINEAR_SWITCH_SCNDINPUT { \
  switch( I1->type ) {             \
  default :                        \
    if ( _verbose_ )               \
      fprintf( stderr, "%s: such input image #1 type not handled yet\n", proc ); \
    return( -1 );                  \
  case SCHAR :                     \
    _LINEAR_SCNDINPUT( s8 );       \
    break;                         \
  case UCHAR :                     \
    _LINEAR_SCNDINPUT( u8 );       \
    break;                         \
  case SSHORT :                    \
    _LINEAR_SCNDINPUT( s16 );      \
    break;                         \
  case USHORT :                    \
    _LINEAR_SCNDINPUT( u16 );      \
    break;                         \
  case FLOAT :                     \
    _LINEAR_SCNDINPUT( r32 );      \
    break;                         \
  }                                \
}

  switch( I0->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such input image #0 type not handled yet\n", proc );
    return( -1 );
  case SCHAR :
    {
      s8 *buf0 = (s8*)I0->data;
      _LINEAR_SWITCH_SCNDINPUT;
    }
    break;
  case UCHAR :
    {
      u8 *buf0 = (u8*)I0->data;
      _LINEAR_SWITCH_SCNDINPUT;
    }
    break;
  case SSHORT :
    {
      s16 *buf0 = (s16*)I0->data;
      _LINEAR_SWITCH_SCNDINPUT;
    }
    break;
  case USHORT :
    {
      u16 *buf0 = (u16*)I0->data;
      _LINEAR_SWITCH_SCNDINPUT;
    }
    break;
  case FLOAT :
    {
      r32 *buf0 = (r32*)I0->data;
      _LINEAR_SWITCH_SCNDINPUT;
    }
    break;
  }

  return( 1 );
}





/**********************************************************************
 *
 * label and couple management stuff
 *
 **********************************************************************/



typedef struct _wlabel {
  int i;
  float d;
} _wlabel;



typedef struct _wlabelList {
  _wlabel *data;
  int n_data;
  int n_allocated_data;
} _wlabelList;



static void _initWLabelList( _wlabelList *l )
{
  l->data = (_wlabel*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeWLabelList( _wlabelList *l )
{
  if ( l->data != (_wlabel*)NULL )
    vtfree( l->data );
  _initWLabelList( l );
}



static int _wlabels_to_be_allocated_ = 3;

static int _addWLabelToWLabelList( _wlabelList *l, int label, float d )
{
  char *proc = "_addWLabelToWLabelList";
  int s =  l->n_allocated_data;
  int i;
  _wlabel *data;

  if ( l->n_data > 0 ) {
    for ( i=0; i<l->n_data; i++ ) {
      if ( l->data[i].i == label ) {
        if ( l->data[i].d > d ) l->data[i].d = d;
        return( 1 );
      }
    }
  }

  if ( l->n_data == l->n_allocated_data ) {
    s += _wlabels_to_be_allocated_;
    data = (_wlabel*)vtmalloc( s * sizeof(_wlabel), "data", proc );
    if ( data == (_wlabel*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_wlabel) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data].i = label;
  l->data[l->n_data].d = d;
  l->n_data ++;

  return( 1 );
}



static int _compareLabels( const void * a, const void * b )
{
  _wlabel *la = (_wlabel *)a;
  _wlabel *lb = (_wlabel *)b;
  if ( la->d < lb->d ) return( -1 );
  if ( la->d > lb->d ) return( 1 );
  return( 0 );
}



static void _fprintfWLabelList( FILE *f, _wlabelList *l, char *d )
{
  int i;
  if ( d != (char*)NULL )
    fprintf( f, "%s=", d );
  fprintf( f, "[" );
  for ( i=0; i<l->n_data; i++ ) {
    fprintf( f, "%d", l->data[i].i );
    if ( i < l->n_data-1 )
      fprintf( stderr, "," );
  }
  fprintf( f, "]" );
}



typedef struct _wcouple {
  int l0;
  int l1;
  float d0;
  float d1;
  float c;
} _wcouple;



typedef struct _wcoupleList {
  _wcouple *data;
  int n_data;
  int n_allocated_data;
} _wcoupleList;



static void _initWCoupleList( _wcoupleList *l )
{
  l->data = (_wcouple*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void _freeWCoupleList( _wcoupleList *l )
{
  if ( l->data != (_wcouple*)NULL )
    vtfree( l->data );
  _initWCoupleList( l );
}



static int _wcouples_to_be_allocated_ = 3;

static int _addWCoupleToWCoupleList( _wcoupleList *l, _wlabel *l0, _wlabel *l1, float t )
{
  char *proc = "_addWCoupleToWCoupleList";
  int s =  l->n_allocated_data;
  _wcouple *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _wcouples_to_be_allocated_;
    data = (_wcouple*)vtmalloc( s * sizeof(_wcouple), "data", proc );
    if ( data == (_wcouple*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(_wcouple) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }


  l->data[l->n_data].l0 = l0->i;
  l->data[l->n_data].d0 = l0->d;
  l->data[l->n_data].l1 = l1->i;
  l->data[l->n_data].d1 = l1->d;

  if ( t < 0.001 )
    l->data[l->n_data].c = l1->d;
  else if ( t > 0.999 )
    l->data[l->n_data].c = l0->d;
  else
    l->data[l->n_data].c = (1.0 - t) * l0->d + t * l1->d;
  l->n_data ++;

  return( 1 );
}



static int _compareCouples( const void * a, const void * b )
{
  _wcouple *ca = (_wcouple *)a;
  _wcouple *cb = (_wcouple *)b;
  if ( ca->c < cb->c ) return( -1 );
  if ( ca->c > cb->c ) return( 1 );
  if ( ca->l0 < cb->l0 ) return( -1 );
  if ( ca->l0 > cb->l0 ) return( 1 );
  return( 0 );
}





/**********************************************************************
 *
 * statistics stuff
 *
 **********************************************************************/



typedef struct _typeStats {
  int *iterations;
  int niterations;
  int nunset;
  float rmax0;
  float rmax1;
} _typeStats;



static void _initStats( _typeStats *s )
{
  s->iterations = (int*)NULL;
  s->niterations = 0;
  s->nunset = 0;
  s->rmax0 = -1.0;
  s->rmax1 = -1.0;
}



static void _freeStats( _typeStats *s )
{
  if ( s->iterations != (int*)NULL )
     vtfree( s->iterations );
}


static int _addStats( _typeStats *s, int i, float rmax0, float rmax1 )
{
  char *proc = "_addStats";
  int *data;
  int j;

  if ( i >= s->niterations ) {
    data = (int*)vtmalloc( (i+1)*sizeof(int), "data", proc );
    if ( data == (int*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    for ( j=0; j<=i; j++ ) data[j] = 0;
    if ( s->iterations != (int*)NULL ) {
      for ( j=0; j<s->niterations; j++ )
        data[j] = s->iterations[j];
      vtfree( s->iterations );
    }
    s->iterations = data;
    s->niterations = i+1;
  }

  s->iterations[i] ++;
  if ( i>0 ) s->iterations[i-1] --;

  if ( s->rmax0 < rmax0 ) s->rmax0 = rmax0;
  if ( s->rmax1 < rmax1 ) s->rmax1 = rmax1;

  return( 1 );
}



static void _fprintfStats( FILE *f, _typeStats *s )
{
  int i;
  int sum = 0;

  fprintf( f, "--- stats ---\n" );

  if ( s->iterations != (int*)NULL ) {
    for ( i=0; i<s->niterations; i++ ) {
      if ( s->iterations[i] > 0 ) {
        sum += s->iterations[i] ;
      }
    }

    if ( s->nunset > 0 ) s->iterations[s->niterations-1] -= s->nunset;

    for ( i=0; i<s->niterations; i++ ) {
      if ( s->iterations[i] > 0 ) {
        fprintf( f, "   - #[iterations = %2d] = %5.2f = %8d/%8d\n",
                 i, 100.0 * (float)s->iterations[i]/(float)sum,
                 s->iterations[i], sum );
      }
    }
  }
  if ( s->nunset > 0 )
    fprintf( f, "   - no label set       = %5.2f = %8d/%8d points\n",
             100.0 * (float)s->nunset/(float)sum, s->nunset, sum );
  fprintf( f, "   - maximum radius in image #0 = %f\n", s->rmax0 );
  fprintf( f, "   - maximum radius in image #1 = %f\n", s->rmax1 );
}



/**********************************************************************
 *
 * label interpolation stuff
 *
 **********************************************************************/



static int _transformPoint( bal_integerPoint *ipt, bal_floatPoint *thePt,
                     bal_floatPoint *thePt0, bal_floatPoint *thePt1,
                     bal_image *theLabel0, bal_image *theLabel1,
                     bal_transformation *theTrsf, bal_transformation *auxTrsf,
                     float t )
{
  char *proc = "_transformPoint";
  float x, y, z = 0.0;

  if ( BAL_IsTransformationLinear( theTrsf ) ) {
     if ( BAL_TransformFloatPoint( thePt, thePt1, theTrsf ) != 1 ) {
       if ( _verbose_ )
         fprintf( stderr, "%s: unable to transform point towards #1\n", proc );
       return( -1 );
     }
     if ( BAL_TransformFloatPoint( thePt, thePt0, auxTrsf ) != 1 ) {
       if ( _verbose_ )
         fprintf( stderr, "%s: unable to transform point towards #0\n", proc );
       return( -1 );
     }
     thePt0->x /= theLabel0->vx;
     thePt0->y /= theLabel0->vy;
     thePt0->z /= theLabel0->vz;
     thePt1->x /= theLabel1->vx;
     thePt1->y /= theLabel1->vy;
     thePt1->z /= theLabel1->vz;
     return( 1 );
  }

  switch( theTrsf->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
  case VECTORFIELD_3D :
    switch( theTrsf->vz.type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such z coordinate type not handled yet\n", proc );
      return( -1 );
    case FLOAT :
      {
        r32 ***theZ = (r32 ***)theTrsf->vz.array;
        z = theZ[ipt->z][ipt->y][ipt->x];
      }
      break;
    }
  case VECTORFIELD_2D :
    switch( theTrsf->vy.type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such y coordinate type not handled yet\n", proc );
      return( -1 );
    case FLOAT :
      {
        r32 ***theY = (r32 ***)theTrsf->vy.array;
        y = theY[ipt->z][ipt->y][ipt->x];
      }
      break;
    }
    switch( theTrsf->vx.type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such x coordinate type not handled yet\n", proc );
      return( -1 );
    case FLOAT :
      {
        r32 ***theX = (r32 ***)theTrsf->vx.array;
        x = theX[ipt->z][ipt->y][ipt->x];
      }
      break;
    }
    break;
  }


  switch( theTrsf->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
  case VECTORFIELD_3D :
    thePt1->x = thePt->x + x;
    thePt1->y = thePt->y + y;
    thePt1->z = thePt->z + z;
    thePt0->x = thePt->x - t/(1.0-t) * x;
    thePt0->y = thePt->y - t/(1.0-t) * y;
    thePt0->z = thePt->z - t/(1.0-t) * z;
    break;
  case VECTORFIELD_2D :
    thePt1->x = thePt->x + x;
    thePt1->y = thePt->y + y;
    thePt1->z = thePt->z;
    thePt0->x = thePt->x - t/(1.0-t) * x;
    thePt0->y = thePt->y - t/(1.0-t) * y;
    thePt0->z = thePt->z;
    break;
  }

  thePt0->x /= theLabel0->vx;
  thePt0->y /= theLabel0->vy;
  thePt0->z /= theLabel0->vz;
  thePt1->x /= theLabel1->vx;
  thePt1->y /= theLabel1->vy;
  thePt1->z /= theLabel1->vz;

  return( 1 );
}



/* return the number of found labels,
 * or -1 in case of error
 */
static int _getLabels( _wlabelList *list,
                       bal_image *theLabel, bal_floatPoint *thePt, float rmax )
{
  char *proc = "_getLabels";
  int xmin, xmax, ymin, ymax, zmin, zmax;
  int i, j, k;
  float d;

  list->n_data = 0;

  /* get nearest point if any
   * else put 0 if the nearest is outside
   */

#define _GETNEARESTLABEL( TYPE ) { \
  TYPE ***buf = (TYPE***)theLabel->array; \
  if ( _addWLabelToWLabelList( list, buf[k][j][i], 0 ) != 1 ) { \
    if ( _verbose_ )             \
      fprintf( stderr, "%s: unable to add label\n", proc ); \
    return( -1 );                \
  }                              \
}

  i = (int)( thePt->x + 0.5 );
  j = (int)( thePt->y + 0.5 );
  k = (int)( thePt->z + 0.5 );

  if ( i >= 0 && i < (int)theLabel->ncols
       && j >= 0 && j < (int)theLabel->nrows
       && k >= 0 && k < (int)theLabel->nplanes ) {
    switch( theLabel->type ) {
    default :
      if ( _verbose_ )
         fprintf( stderr, "%s: such label image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _GETNEARESTLABEL( u8 );
      break;
    case SSHORT :
      _GETNEARESTLABEL( s16 );
      break;
    case USHORT :
      _GETNEARESTLABEL( u16 );
      break;
    }
  }
  else {
    if ( _addWLabelToWLabelList( list, 0, 0.0 ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add label\n", proc );
      return( -1 );
    }
  }

  /* stop if no neighborhood search
   */
  if ( rmax <= 0.0 )
    return( list->n_data );

  /* search over neighborhood
   */

  xmin = (int)( thePt->x - rmax + 0.5 );
  xmax = (int)( thePt->x + rmax + 0.5 );
  if ( xmin < 0 ) {
    if ( xmax < 0 ) return( list->n_data );
    xmin = 0;
    d = ((-1)-thePt->x)*((-1)-thePt->x);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }
  if ( xmax >= (int)theLabel->ncols ) {
    if ( xmin >= (int)theLabel->ncols ) return( list->n_data );
    xmax = theLabel->ncols-1;
    d = (theLabel->ncols-thePt->x)*(theLabel->ncols-thePt->x);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }

  ymin = (int)( thePt->y - rmax + 0.5 );
  ymax = (int)( thePt->y + rmax + 0.5 );
  if ( ymin < 0 ) {
    if ( ymax < 0 ) return( list->n_data );
    ymin = 0;
    d = ((-1)-thePt->y)*((-1)-thePt->y);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }
  if ( ymax >= (int)theLabel->nrows ) {
    if ( ymin >= (int)theLabel->nrows ) return( list->n_data );
    ymax = theLabel->nrows-1;
    d = (theLabel->nrows-thePt->y)*(theLabel->nrows-thePt->y);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }

  zmin = (int)( thePt->z - rmax + 0.5 );
  zmax = (int)( thePt->z + rmax + 0.5 );
  if ( zmin < 0 ) {
    if ( zmax < 0 ) return( list->n_data );
    zmin = 0;
    d = ((-1)-thePt->z)*((-1)-thePt->z);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }
  if ( zmax >= (int)theLabel->nplanes ) {
    if ( zmin >= (int)theLabel->nplanes ) return( list->n_data );
    zmax = theLabel->nplanes-1;
    d = (theLabel->nplanes-thePt->z)*(theLabel->nplanes-thePt->z);
    if ( d <= rmax*rmax ) {
      if ( _addWLabelToWLabelList( list, 0, sqrt(d) ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add label\n", proc );
        return( -1 );
      }
    }
  }

#define _GETLABELS( TYPE ) {       \
  TYPE ***buf = (TYPE***)theLabel->array; \
  for ( k=zmin; k<=zmax; k++ )     \
  for ( j=ymin; j<=ymax; j++ )     \
  for ( i=xmin; i<=xmax; i++ ) {   \
    d = (k-thePt->z)*(k-thePt->z) + (j-thePt->y)*(j-thePt->y) + (i-thePt->x)*(i-thePt->x); \
    if ( d > rmax*rmax ) continue; \
    if ( _addWLabelToWLabelList( list, buf[k][j][i], sqrt(d) ) != 1 ) { \
      if ( _verbose_ )             \
        fprintf( stderr, "%s: unable to add label\n", proc ); \
      return( -1 );                \
    }                              \
  }                                \
}

  switch( theLabel->type ) {
  default :
    if ( _verbose_ )
       fprintf( stderr, "%s: such label image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _GETLABELS( u8 );
    break;
  case SSHORT :
    _GETLABELS( s16 );
    break;
  case USHORT :
    _GETLABELS( u16 );
    break;
  }

  if ( list->n_data > 1 )
    qsort( (void*)(list->data), list->n_data, sizeof(_wlabel), &_compareLabels );

  return( list->n_data );
}



static int _buildCouples( _wcoupleList *l,
                          _wlabelList *list0,
                          _wlabelList *list1,
                          float t )
{
  char *proc = "_buildCouples";
  int i, j;

  l->n_data = 0;

  if ( list0->n_data > 0 && list1->n_data > 0 ) {
    for ( i=0; i<list0->n_data; i++ )
    for ( j=0; j<list1->n_data; j++ ) {
      if ( _addWCoupleToWCoupleList( l, &(list0->data[i]), &(list1->data[j]), t ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add couple\n", proc );
        return( -1 );
      }
    }
  }
  else {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird, this points should not be reached\n", proc );
    return( -1 );
  }

  if ( l->n_data > 1 )
    qsort( (void*)(l->data), l->n_data, sizeof(_wcouple), &_compareCouples );

  return( 1 );
}



static _wcouple *_getCouple( _wcoupleList *l, typeCellCorrespondence *cc )
{
  int i, j;
  typeLabelList *ll;

  for ( i=0; i<l->n_data; i++ ) {
    /* if one of the label is outside (ie 0)
     * it is a valid couple
     */
    if ( l->data[i].l0 == 0 || l->data[i].l1 == 0 ) {
      return( &(l->data[i]) );
    }
    /* label list for label #0 l0
     */
    ll = &(cc->list0->data[ l->data[i].l0 ].corresp);
    for ( j=0; j<ll->n_data; j++ ) {
      if ( ll->data[j] == l->data[i].l1 ) {
        return( &(l->data[i]) );
      }
    }
  }

  return( (_wcouple *)NULL );
}



static int _setLabels( _wcouple *c, bal_integerPoint *ipt,
                       bal_image *resLabel0, bal_image *resLabel1 )
{
  char *proc = "_setLabels";

#define _SETLABEL( TYPE, LABEL, L ) { \
  TYPE ***buf = (TYPE***)LABEL->array; \
  buf[ipt->z][ipt->y][ipt->x] = L;    \
}

  switch( resLabel0->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such label image not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _SETLABEL( u8, resLabel0, c->l0 );
    break;
  case SSHORT :
    _SETLABEL( s16, resLabel0, c->l0 );
    break;
  case USHORT :
    _SETLABEL( u16, resLabel0, c->l0 );
    break;
  }

  switch( resLabel1->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such label image not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _SETLABEL( u8, resLabel1, c->l1 );
    break;
  case SSHORT :
    _SETLABEL( s16, resLabel1, c->l1 );
    break;
  case USHORT :
    _SETLABEL( u16, resLabel1, c->l1 );
    break;
  }

  return( 1 );
}



static void _setWarning( _wlabelList *label, typeCellList *l, char *desc )
{
  int i;
  typeCell *c;

  for ( i=0; i<label->n_data; i++ ) {
    c = &(l->data[ label->data[i].i ]);
    if ( (c->left_corner.x >= 0 || c->left_corner.y >= 0
          || c->left_corner.z >= 0) && c->corresp.n_data <= 0 ) {
        fprintf( stderr, "[pb, ");
        if ( desc != (char*)NULL )
          fprintf( stderr, "'%s' ", desc );
        fprintf( stderr, "= %d]", label->data[i].i );
    }
  }
  return;
}



static int _processPoint( bal_image *theLabel0, bal_image *theLabel1,
                          typeCellCorrespondence *cc,
                          bal_image *resLabel0, bal_image *resLabel1,
                          bal_integerPoint *ipt,
                          bal_floatPoint *thePt0, bal_floatPoint *thePt1,
                          _wlabelList *label0, _wlabelList *label1,
                          _wcoupleList *couples,
                          float rmax0, float rmax1, float t )
{
  char *proc = "_processPoint";
  _wcouple *c;

  /* get labels from both label images
   * there is at least one label,
   * 0 (outside) is in the labels if the nearest point was outside
   */
  if ( _getLabels( label0, theLabel0, thePt0, rmax0 ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when getting labels from image #0\n", proc );
    return( -1 );
  }

  if ( _getLabels( label1, theLabel1, thePt1, rmax1 ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when getting labels from image #0\n", proc );
    return( -1 );
  }

  if ( 0 ) {
    _setWarning( label0, cc->list0, "#0" );
    _setWarning( label1, cc->list1, "#1" );
  }

  if ( 0 ) {
    _fprintfWLabelList( stderr, label0, "l0" );
    _fprintfWLabelList( stderr, label1, "l1" );
  }

  /* building couples
   */
  if ( _buildCouples( couples, label0, label1, t ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build couples of labels\n", proc );
    return( -1 );
  }

  /* is there a valid couple ?
   */
  c = _getCouple( couples, cc );
  if ( c == (_wcouple *)NULL ) return( 0 );

  /* set labels
   */
  if ( _debug_ )
    fprintf( stderr, "found(%d,%d,%d)->[%d-%d]", ipt->x, ipt->y, ipt->z, c->l0, c->l1 );

  if ( _setLabels( c, ipt, resLabel0, resLabel1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to set labels\n", proc );
    return( -1 );
  }

  return( 1 );
}


/* returned values
 * -1 : error
 * 0 : condition not fulfilled, no other candidate
 * 1 : condition fulfilled
 */
static int _testOtherCandidates( int l, float d,
                                 _wlabelList *label0,
                                 _wlabelList *label1,
                                 typeCellList *corr,
                                 float *dmin )
{
  int i, j, k;
  int lj, found;

  *dmin = -1.0;

  for ( i=0; i<label0->n_data; i++ ) {
    if ( label0->data[i].d >= d ) continue;
    if ( label0->data[i].i == l ) continue;
    /* here we have an other label than l with a smaller distance
     */
    for ( j=0; j<corr->data[l].corresp.n_data; j++ ) {
      lj = corr->data[l].corresp.data[j];
      for ( found=0, k=0; k<label1->n_data; k++ ) {
        if ( label1->data[k].i == lj ) {
          found = 1;
        }
      }
      if ( found == 0 ) {
        *dmin = label0->data[i].d;
        return( 1 );
      }
    }
  }

  return( 0 );
}







/* theTrsf is the transformation from t to 1
 * T_{1<-t}
 */
static int BAL_LabelInterpolation( bal_transformation *theTrsf,
                            bal_image *theLabel0, bal_image *theLabel1,
                            typeCellCorrespondence *cc,
                            bal_image *resLabel0, bal_image *resLabel1,
                            float t, float rmax )
{
  char *proc = "BAL_LabelInterpolation";
  bal_transformation auxTrsf;
  bal_transformation *ptrTrsf = (bal_transformation*)NULL;
  bal_integerPoint ipt;
  bal_floatPoint thePt, thePt0, thePt1;
  _wlabelList label0;
  _wlabelList label1;
  _wcoupleList couples;
  int ret;

  _wcouple out;

  float dr0 = 0.0;
  float dr1 = 0.0;
  float r0 = -1.0;
  float r1 = -1.0;
  float rmax0 = -1.0;
  float rmax1 = -1.0;
  float dmin0, dmin1;
  int test0, test1;
  int i;

  _typeStats stats;




  BAL_InitTransformation( &auxTrsf );
  if ( BAL_IsTransformationLinear( theTrsf ) ) {
    if ( BAL_AllocTransformation( &auxTrsf, theTrsf->type, (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate inverse transformation (linear case)\n", proc );
      return( -1 );
    }
    if ( BAL_CopyTransformation( theTrsf, &auxTrsf ) != 1 ) {
      BAL_FreeTransformation( &auxTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to copy transformation (linear case)\n", proc );
      return( -1 );
    }
    if ( BAL_MultiplyTrsf( &auxTrsf, &auxTrsf, -t/(1.0-t) ) != 1 ) {
      BAL_FreeTransformation( &auxTrsf );
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to multiply transformation\n", proc );
      }
      return( -1 );
    }
    ptrTrsf = &auxTrsf;
  }

  /* initialization
   */
  _initWLabelList( &label0 );
  _initWLabelList( &label1 );
  _initWCoupleList( &couples );

  out.l0 = 0;
  out.l1 = 0;
  out.c = 0.0;

  if ( t < 0.001 ) {
    dr0 = 0.0;
    dr1 = 1.0;
  }
  else if ( t > 0.999 ) {
    dr0 = 1.0;
    dr1 = 0.0;
  }
  else {
    if ( t < 0.5 ) {
      dr0 = t / (1.0 - t);
      dr1 = 1.0;
    }
    else {
        dr0 = 1.0;
        dr1 = (1.0 - t) / t;
    }
  }

  _initStats( &stats );

  for ( ipt.z=0; ipt.z<(int)resLabel0->nplanes; ipt.z++ )
  for ( ipt.y=0; ipt.y<(int)resLabel0->nrows; ipt.y++ )
  for ( ipt.x=0; ipt.x<(int)resLabel0->ncols; ipt.x++ ) {


    thePt.x = ipt.x * resLabel0->vx;
    thePt.y = ipt.y * resLabel0->vy;
    thePt.z = ipt.z * resLabel0->vz;

    /* thePt0 and thePt1 are the points in voxel
     * coordinates in both theLabel0 and theLabel1
     */
    if ( _transformPoint( &ipt, &thePt, &thePt0, &thePt1,
                          theLabel0, theLabel1, theTrsf, ptrTrsf, t ) != 1 ) {
      _freeWLabelList( &label0 );
      _freeWLabelList( &label1 );
      _freeWCoupleList( &couples );
      if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to transform point\n", proc );
      return( -1 );
    }


    /* initializations
     */
    r0 = r1 = -1.0;
    rmax0 = rmax1 = -1.0;
    i = 0;


    if ( _debug_ ) {
      fprintf( stderr, "... process (%d,%d,%d)", ipt.x, ipt.y, ipt.z );
    }


    do {

        if ( _addStats( &stats, i, rmax0, rmax1 ) != 1 ) {
          _freeStats( &stats );
          _freeWLabelList( &label0 );
          _freeWLabelList( &label1 );
          _freeWCoupleList( &couples );
          if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to update stats for point (%d,%d,%d)\n", proc, ipt.x, ipt.y, ipt.z );
          return( -1 );
        }
        i++;

        /* returned value
         * -1 = error
         * 0 = no point set
         * try first with the nearest neighbors
         */

        ret = _processPoint( theLabel0, theLabel1, cc,
                             resLabel0, resLabel1, &ipt, &thePt0, & thePt1,
                             &label0, &label1, &couples, r0, r1, t );
        if ( ret == -1 ) {
          _freeStats( &stats );
          _freeWLabelList( &label0 );
          _freeWLabelList( &label1 );
          _freeWCoupleList( &couples );
          if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to set labels for point (%d,%d,%d)\n", proc, ipt.x, ipt.y, ipt.z );
          return( -1 );
        }


        /* first search with nearest neighbors
         * stop if a couple is found
         */
        if ( ret == 1 && r0 < 0.0 && r1 < 0.0 ) {
          if ( _debug_ )
            fprintf( stderr, "nearest=[%d,%d]\n", couples.data[0].l0, couples.data[0].l1 );
          break;
        }

        /* emergency stop
         * with a user-given rmax
         */
        if ( rmax > 0 && (r0 > rmax || r1 > rmax) ) {
          stats.nunset ++;
          if ( _verbose_ >= 2 )
            fprintf( stderr, "%s: no labels set for point (%d,%d,%d)\n", proc, ipt.x, ipt.y, ipt.z );
          if ( _setLabels( &out, &ipt, resLabel0, resLabel1 ) != 1 ) {
            _freeStats( &stats );
            _freeWLabelList( &label0 );
            _freeWLabelList( &label1 );
            _freeWCoupleList( &couples );
            if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to set labels for point (%d,%d,%d)\n", proc, ipt.x, ipt.y, ipt.z );
            return( -1 );
          }
          if ( _debug_ )
            fprintf( stderr, " - emergency stop\n" );
          break;
        }

        /* no couple have been found so far
         * increase radii
         */
        if ( ret == 0 ) {
          r0 = (r0 < 0) ? dr0 : r0 + dr0;
          r1 = (r1 < 0) ? dr1 : r1 + dr1;
          if ( _debug_ )
            fprintf( stderr, "(%f,%f)", r0, r1 );
          continue;
        }

        /* case ret=1
         * a couple has been found
         * test whether other candidates may be found
         */

        test0 =_testOtherCandidates( couples.data[0].l0, couples.data[0].d0,
                                   &label0, &label1, cc->list0, &dmin0 );
        if ( test0 == -1 ) {
          _freeStats( &stats );
          _freeWLabelList( &label0 );
          _freeWLabelList( &label1 );
          _freeWCoupleList( &couples );
          if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to test label #0 %d for point (%d,%d,%d)\n", proc,
                     couples.data[0].l0, ipt.x, ipt.y, ipt.z );
          return( -1 );
        }
        if ( test0 == 1 ) {
          rmax1 = couples.data[0].d1 + (1-t)/t * ( couples.data[0].d0 - dmin0 );
        }

        test1 = _testOtherCandidates( couples.data[0].l1, couples.data[0].d1,
                                   &label1, &label0, cc->list1, &dmin1 );
        if ( test1 == -1 ) {
          _freeStats( &stats );
          _freeWLabelList( &label0 );
          _freeWLabelList( &label1 );
          _freeWCoupleList( &couples );
          if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to test label #1 %d for point (%d,%d,%d)\n", proc,
                     couples.data[0].l1, ipt.x, ipt.y, ipt.z );
          return( -1 );
        }
        if ( test1 == 1 ) {
          rmax0 = couples.data[0].d0 + t/(1-t) * ( couples.data[0].d1 - dmin1 );
        }

        /* no possible candidates
         */
        if ( test0 == 0 && test1 == 0 ) {
          if ( _debug_ )
            fprintf( stderr, "tests failed\n" );
          break;
        }

        /* increase search radii
         * or stop if bounds are raeched
         */

        if ( rmax0 < 0.0 ) {
          if ( rmax1 < 0.0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, this should not be reached\n", proc );
            if ( _debug_ )
              fprintf( stderr, "\n" );
            break;
          }
          else {
            if ( r1 > rmax1 ) {
              if ( _debug_ )
                fprintf( stderr, "\n" );
              break;
            }
            r1 = (r1 < 0) ? dr1 : r1 + dr1;
          }
        }
        else {
          if ( rmax1 < 0.0 ) {
            if ( r0 > rmax0 ) {
              if ( _debug_ )
                fprintf( stderr, "\n" );
              break;
            }
            r0 = (r0 < 0) ? dr0 : r0 + dr0;
          }
          else {
            if ( r0 > rmax0 && r1 > rmax1 ) {
              if ( _debug_ )
                fprintf( stderr, "\n" );
              break;
            }
            if ( r0 <= rmax0 ) r0 = (r0 < 0) ? dr0 : r0 + dr0;
            if ( r1 <= rmax1 ) r1 = (r1 < 0) ? dr1 : r1 + dr1;
          }
        }

        if ( _debug_ )
          fprintf( stderr, "(%f/%f,%f/%f)", r0, rmax0, r1, rmax1 );

    } while( 1 );

  }

  if ( _verbose_ )
    _fprintfStats( stderr, &stats );


  _freeStats( &stats );
  _freeWLabelList( &label0 );
  _freeWLabelList( &label1 );
  _freeWCoupleList( &couples );
  if ( ptrTrsf != (bal_transformation*)NULL ) BAL_FreeTransformation( &auxTrsf );
  return( 1 );
}





/**********************************************************************
 *
 * interpolation
 *
 **********************************************************************/


typedef enum _typeInterpolation {
  _GREYLEVEL_,
  _LABELS_
} _typeInterpolation ;


/* transformation is supposed to go from #1 to #0
 * ie T_{0 <- 1}
 * resImage0 is theImage0 resampled at t
 * resImage1 is theImage1 resampled at t
 * resImage = (1-t) * resImage0 + t * resImage1;
 */
static int _interpolateImages( bal_image *theImage0, bal_image *theImage1,
                               bal_image *resImage0, bal_image *resImage1,
                               bal_image *resImage,
                               bal_transformation *theTrsf,
                               typeCellCorrespondence *cc,
                               float t,
                               enumTransformationInterpolation theInterpolation,
                               _typeInterpolation interpolationMode,
                               float rmax )
{
  char *proc = "_interpolateImages";
  bal_image *ptrRes0 = (bal_image*)NULL;
  bal_image *ptrRes1 = (bal_image*)NULL;
  bal_image *ptrTmpRes0 = (bal_image*)NULL;
  bal_image *ptrTmpRes1 = (bal_image*)NULL;
  bal_image tmpRes0;
  bal_image tmpRes1;
  bal_image *ptrTemplate = (bal_image*)NULL;

  bal_transformation invTrsf;
  bal_transformation partialTrsf;

  enumTransformationInterpolation interpolation = theInterpolation;

  if ( t < 0.0 || t > 1.0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: t out of bounds\n", proc );
    return( -1 );
  }

  switch( interpolationMode ) {
  default :
  case _GREYLEVEL_ :
    if ( resImage0 == (bal_image*)NULL
         && resImage1 == (bal_image*)NULL && resImage == (bal_image*)NULL )
      return( 1 );
  case _LABELS_ :
    if ( resImage0 == (bal_image*)NULL && resImage1 == (bal_image*)NULL )
      return( 1 );
  }

  if ( resImage != (bal_image*)NULL )
    ptrTemplate = resImage;
  else if ( resImage0 != (bal_image*)NULL )
    ptrTemplate = resImage0;
  else if ( resImage1 != (bal_image*)NULL )
    ptrTemplate = resImage1;

  switch( interpolationMode ) {
  default :
    break;
  case _LABELS_ :
    interpolation = NEAREST;
  }



  /* auxiliary image allocation
   */
  if ( (resImage != (bal_image*)NULL && interpolationMode != _LABELS_)
       || interpolationMode == _LABELS_ ) {
    if ( resImage0 == (bal_image*)NULL ) {
      if ( BAL_AllocImageFromImage( &tmpRes0, (char*)NULL,
                                        resImage, theImage0->type ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary result image #0\n", proc );
        return( -1 );
      }
      ptrTmpRes0 = &tmpRes0;
      ptrRes0 = &tmpRes0;
    }
    else {
      ptrRes0 = resImage0;
    }
    if ( resImage1 == (bal_image*)NULL ) {
      if ( BAL_AllocImageFromImage( &tmpRes1, (char*)NULL,
                                        resImage, theImage1->type ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to allocate auxiliary result image #1\n", proc );
        return( -1 );
      }
      ptrTmpRes1 = &tmpRes1;
      ptrRes1 = &tmpRes1;
    }
    else {
      ptrRes1 = resImage1;
    }
  }
  else {
    if ( resImage0 != (bal_image*)NULL )  ptrRes0 = resImage0;
    if ( resImage1 != (bal_image*)NULL )  ptrRes1 = resImage1;
  }


  /************************************************************
   *
   ************************************************************/

  if ( t > 0.99 ) {

    /* here, we are at t = 1,
     * recall that the transformation is supposed to be from 1 to 0
     * ie T_{0<-1}
     */

    switch( interpolationMode ) {
    default :
    case _GREYLEVEL_ :
     /* we just have to copy image1, and to resample image0
     */
      if ( ptrRes1 != (bal_image*)NULL ) {
        if ( BAL_CopyImage( theImage1, ptrRes1 ) != 1 ) {
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy image #0\n", proc );
          return( -1 );
        }
      }

      if ( ptrRes0 != (bal_image*)NULL ) {
        if ( BAL_ResampleImage( theImage0, ptrRes0, theTrsf, interpolation ) != 1 ) {
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to resample image #1\n", proc );
          return( -1 );
        }
      }
      break;
    case _LABELS_ :
      /* we have T_{0<-1}
       * we exchange #0 and #1, thus it becomes T_{1<-0}
       */
      _invertCellCorrespondence( cc );
      if ( BAL_LabelInterpolation( theTrsf, theImage1, theImage0,
                                   cc, ptrRes1, ptrRes0, 1.0-t, rmax ) != 1 ) {
        if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
        if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to interpolate label images %f\n", proc, t );
        return( -1 );
      }
      _invertCellCorrespondence( cc );
      break;
    }

  }

  else if ( t < 0.01 ) {

    /* here, we are at t = 0,
     * we compute T_{1<-0}
     */
    if ( BAL_AllocTransformation( &invTrsf, theTrsf->type, theImage0 ) != 1 ) {
      if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
      if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary inverse transformation\n", proc );
      return( -1 );
    }
    if ( BAL_InverseTransformation( theTrsf, &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
      if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to invert transformation\n", proc );
      return( -1 );
    }

    switch( interpolationMode ) {
    default :
    case _GREYLEVEL_ :
      /* we just have to copy image0, and to resample image1
       */
      if ( ptrRes0 != (bal_image*)NULL ) {
        if ( BAL_CopyImage( theImage0, ptrRes0 ) != 1 ) {
          BAL_FreeTransformation( &invTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to copy image #0\n", proc );
          return( -1 );
        }
      }

      if ( ptrRes1 != (bal_image*)NULL ) {
        if ( BAL_ResampleImage( theImage1, ptrRes1, &invTrsf, interpolation ) != 1 ) {
          BAL_FreeTransformation( &invTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to resample image #1\n", proc );
          return( -1 );
        }
      }
      break;
    case _LABELS_ :
      if ( BAL_LabelInterpolation( &invTrsf, theImage0, theImage1,
                                   cc, ptrRes0, ptrRes1, t, rmax ) != 1 ) {
        BAL_FreeTransformation( &invTrsf );
        if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
        if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to interpolate label images at t=%f\n", proc, t );
        return( -1 );
      }
      break;
    }
    BAL_FreeTransformation( &invTrsf );
  }

  else {

    if ( BAL_AllocTransformation( &partialTrsf, theTrsf->type, ptrTemplate ) != 1 ) {
      if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
      if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate intermediary transformation\n", proc );
      return( -1 );
    }

    /* here, we compute T_{1<-t} from T_{0<-1}
     */
    if ( BAL_ComputeTrsfAtIntermediaryTime( theTrsf, &partialTrsf, t, _TOWARDS_1 ) != 1 ) {
      BAL_FreeTransformation( &partialTrsf );
      if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
      if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when computing intermediary transformation at %f\n", proc, t );
      return( -1 );
    }

    switch( interpolationMode ) {
    default :
    case _GREYLEVEL_ :
      if ( ptrRes1 != (bal_image*)NULL ) {
        if ( BAL_ResampleImage( theImage1, ptrRes1, &partialTrsf, interpolation ) != 1 ) {
          BAL_FreeTransformation( &partialTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to resample image #1 at %f\n", proc, t );
          return( -1 );
        }
      }
      if ( ptrRes0 != (bal_image*)NULL ) {
        if ( BAL_MultiplyTrsf( &partialTrsf, &partialTrsf, -t/(1.0-t) ) != 1 ) {
          BAL_FreeTransformation( &partialTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: error when multiplying transformation by %f\n", proc, -t/(1.0-t) );
          return( -1 );
        }
        if ( BAL_ResampleImage( theImage0, ptrRes0, &partialTrsf, interpolation ) != 1 ) {
          BAL_FreeTransformation( &partialTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to resample image #0 at %f\n", proc, t );
          return( -1 );
        }
      }
      break;
    case _LABELS_ :
        if ( BAL_LabelInterpolation( &partialTrsf, theImage0, theImage1,
                                     cc, ptrRes0, ptrRes1, t, rmax ) != 1 ) {
          BAL_FreeTransformation( &partialTrsf );
          if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
          if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to interpolate label images %f\n", proc, t );
          return( -1 );
        }
      break;
    }

    BAL_FreeTransformation( &partialTrsf );
  }

  if ( resImage != (bal_image*)NULL && interpolationMode != _LABELS_ ) {
    if ( BAL_ImageLinearCombination( 1.0-t, ptrRes0, t, ptrRes1, resImage ) != 1 ) {
      if ( ptrTmpRes1 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes1 );
      if ( ptrTmpRes0 != (bal_image*)NULL ) BAL_FreeImage( ptrTmpRes0 );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to combine interpolated images\n", proc );
      return( -1 );
    }
  }

  return( 1 );
}





int BAL_InterpolateGreyLevelImages( bal_image *theImage0, bal_image *theImage1,
                                    bal_image *resImage0, bal_image *resImage1,
                                    bal_image *resImage,
                                    bal_transformation *theTrsf,
                                    float t,
                                    enumTransformationInterpolation interpolation )
{
  char *proc = "BAL_InterpolateGreyLevelImages";

  if ( _interpolateImages( theImage0, theImage1, resImage0, resImage1, resImage,
                           theTrsf, (typeCellCorrespondence *)NULL,
                           t, interpolation, _GREYLEVEL_, -1.0 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when interpolating at %f\n", proc, t );
    return( -1 );
  }
  return( 1 );
}





int BAL_InterpolateLabelImages( bal_image *theImage0, bal_image *theImage1,
                                bal_image *resImage0, bal_image *resImage1,
                                bal_transformation *theTrsf,
                                typeCellCorrespondence *cc,
                                float t,
                                float rmax )
{
  char *proc = "BAL_InterpolateLabelImages";

  if ( _interpolateImages( theImage0, theImage1, resImage0, resImage1, (bal_image *)NULL,
                           theTrsf, cc,
                           t, NEAREST, _LABELS_, rmax ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when interpolating at %f\n", proc, t );
    return( -1 );
  }
  return( 1 );
}






