/****************************************************
 * topological-voxel-tree.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 11:10:49 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <histogram.h>
#include <vtmalloc.h>

#include <topological-voxel-tree.h>








/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalVoxelTree( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalVoxelTree(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalVoxelTree(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalVoxelTree( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalVoxelTree(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalVoxelTree(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




















/**************************************************
 *
 * intern tree
 * this is an intermediary tree structure
 * dedicated to computation
 *
 **************************************************/








static void _initVoxelTreeIntPoint( voxelTreeIntPoint *p )
{
  p->x = -1;
  p->y = -1;
  p->z = -1;
  p->anchor = 0;
  p->order_branch = 0;
  p->label_branch = 0;
}

static void _initVoxelTreeFloatPoint( voxelTreeFloatPoint *p )
{
  p->x = 0.0;
  p->y = 0.0;
  p->z = 0.0;
  p->anchor = 0;
  p->order_branch = 0;
  p->label_branch = 0;
}

static void _initInfoBranchingPoint( infoBranchingPoint *d )
{
  d->label_branch = 0;
  d->order = 0;
  d->length = 0;
  d->inqueue = 1;
}



void initVoxelTreeComponent( voxelTreeComponent *c )
{
  c->type = _INTERN_TREE_UNKNOWN_COMPONENT_;
  c->n_data = 0;
  c->n_allocated_data = 0;
  c->data = (voxelTreeIntPoint*)NULL;
  _initVoxelTreeFloatPoint( &(c->center) );
  c->minCorner[0] = c->minCorner[1] = c->minCorner[2] = -1;
  c->maxCorner[0] = c->maxCorner[1] = c->maxCorner[2] = -1;
  c->indexFirstJunction = -1;
  c->indexLastJunction = -1;
  c->length = 0.0;
  _initInfoBranchingPoint( &(c->infoFirstPoint) );
  _initInfoBranchingPoint( &(c->infoLastPoint) );
  _initInfoBranchingPoint( &(c->infoJunction) );
  c->canbepruned = 1;
  c->pointIndexOffset = 0;
}



void freeVoxelTreeComponent( voxelTreeComponent *c )
{
  if ( c->data != (voxelTreeIntPoint*)NULL )
      vtfree( c->data );
  initVoxelTreeComponent( c );
}



static int _size_to_be_allocated_ = 10;

int addPointToVoxelTreeComponent( voxelTreeComponent *l, voxelTreeIntPoint *p )
{
  char *proc = "addPointToVoxelTreeComponent";
  int s =  l->n_allocated_data;
  voxelTreeIntPoint *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (voxelTreeIntPoint*)vtmalloc( s * sizeof(voxelTreeIntPoint), "data", proc );
    if ( data == (voxelTreeIntPoint*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(voxelTreeIntPoint) );
      vtfree( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  _initVoxelTreeIntPoint( &(l->data[l->n_data]) );

  l->data[l->n_data] = *p;

  if ( l->n_data == 0 ) {
    l->minCorner[0] = l->maxCorner[0] = p->x;
    l->minCorner[1] = l->maxCorner[1] = p->y;
    l->minCorner[2] = l->maxCorner[2] = p->z;
  }
  else {
    if ( l->minCorner[0] > p->x ) l->minCorner[0] = p->x;
    if ( l->minCorner[1] > p->y ) l->minCorner[1] = p->y;
    if ( l->minCorner[2] > p->z ) l->minCorner[2] = p->z;
    if ( l->maxCorner[0] < p->x ) l->maxCorner[0] = p->x;
    if ( l->maxCorner[1] < p->y ) l->maxCorner[1] = p->y;
    if ( l->maxCorner[2] < p->z ) l->maxCorner[2] = p->z;
  }

  l->n_data ++;

  return( 1 );
}



void  initVoxelTree( typeVoxelTree *t )
{
  t->data = (voxelTreeComponent*)NULL;
  t->n_data = 0;
  t->n_allocated_data = 0;
  t->voxelSize[0] = t->voxelSize[1] = t->voxelSize[2] = 1.0;
}



void freeVoxelTree( typeVoxelTree *t )
{
  int i;
  if ( t->data != (voxelTreeComponent*)NULL ) {
    for ( i=0; i<t->n_data; i++ )
      freeVoxelTreeComponent( &(t->data[i]) );
    vtfree( t->data );
  }
   initVoxelTree( t );
}





static int _allocVoxelTree( typeVoxelTree *t, int n )
{
  char *proc = "_allocVoxelTree";
  int i;

  t->data = (voxelTreeComponent*)vtmalloc( n*sizeof(voxelTreeComponent), "t->data", proc );
  if ( t->data == (voxelTreeComponent*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  t->n_data = n;
  t->n_allocated_data = n;
  for ( i=0; i<n; i++ )
    initVoxelTreeComponent( &(t->data[i]) );

  return( 1 );
}





int getNewComponentFromVoxelTree( typeVoxelTree *t )
{
  char *proc = "getNewComponentFromVoxelTree";
  int i, s =  t->n_allocated_data;
  voxelTreeComponent *data;

  if ( t->n_data < t->n_allocated_data ) {
    t->n_data ++;
    return( t->n_data-1 );
  }

  if ( t->n_data == t->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (voxelTreeComponent*)vtmalloc( s * sizeof(voxelTreeComponent), "data", proc );
    if ( data == (voxelTreeComponent*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    }
    if ( t->n_allocated_data > 0 ) {
      (void)memcpy( data, t->data, t->n_allocated_data*sizeof(voxelTreeComponent) );
      vtfree( t->data );
    }
    t->n_allocated_data = s;
    t->data = data;
    for ( i=t->n_data; i<t->n_allocated_data; i++ )
        initVoxelTreeComponent( &(t->data[i]) );
  }

  t->n_data ++;

  return( t->n_data-1 );
}








/**************************************************
 *
 *
 *
 **************************************************/

void initVoxelTreeList( typeVoxelTreeList *l )
{
  l->data = (typeVoxelTree*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;

}



void freeVoxelTreeList( typeVoxelTreeList *l )
{
  int i;
  if ( l->data != (typeVoxelTree*)NULL ) {
    for ( i=0; i<l->n_data; i++ )
      freeVoxelTree( &(l->data[i]) );
    vtfree( l->data );
  }
   initVoxelTreeList( l );
}



int allocVoxelTreeList( typeVoxelTreeList *l, int n )
{
  char *proc = "allocVoxelTreeList";
  int i;

  l->data = (typeVoxelTree*)vtmalloc( n*sizeof(typeVoxelTree), "l->data", proc );
  if ( l->data == (typeVoxelTree*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  l->n_data = n;
  l->n_allocated_data = n;
  for ( i=0; i<n; i++ )
    initVoxelTree( &(l->data[i]) );

  return( 1 );
}







/**************************************************
 *
 *
 *
 **************************************************/

static void _fprintBBComponentVoxelTree( FILE *f, voxelTreeComponent *c )
{
    fprintf( f, "[%3d %3d %3d]x[%3d %3d %3d], ",
             c->minCorner[0], c->minCorner[1], c->minCorner[2],
             c->maxCorner[0], c->maxCorner[1], c->maxCorner[2] );
}



void fprintEdgeVoxelTree( FILE *f, typeVoxelTree *tree, int i, int writebranches )
{
    int n;

    n = tree->data[i].n_data - 1;
    fprintf( f, "    #%3d: %4d points, ",
             i, tree->data[i].n_data );
    _fprintBBComponentVoxelTree( f, &(tree->data[i]) );
    fprintf( f, " l=%f, ", tree->data[i].length );
    if ( tree->data[i].indexFirstJunction == -1 )
        fprintf( f, "NULL" );
    else
        fprintf( f, "%4d", tree->data[i].indexFirstJunction );
    fprintf( f, " <- " );
    fprintf( f, "[(%d %d %d), ..., (%d %d %d)]",
             tree->data[i].data[0].x, tree->data[i].data[0].y, tree->data[i].data[0].z,
             tree->data[i].data[n].x, tree->data[i].data[n].y, tree->data[i].data[n].z );
    fprintf( f, " -> " );
    if ( tree->data[i].indexLastJunction == -1 )
        fprintf( f, "NULL" );
    else
        fprintf( f, "%4d", tree->data[i].indexLastJunction );
    fprintf( f, "\n" );
    if ( writebranches ) {
        fprintf( f, "\t [b=%d, order=%d,l=%f] <- ... -> [b=%d, order=%d,l=%f]\n",
                 tree->data[i].infoFirstPoint.label_branch,
                 tree->data[i].infoFirstPoint.order,
                 tree->data[i].infoFirstPoint.length,
                 tree->data[i].infoLastPoint.label_branch,
                 tree->data[i].infoLastPoint.order,
                 tree->data[i].infoLastPoint.length );
    }
}



static void _fprintfInfoBranchingPoint( FILE *f, infoBranchingPoint *i, char *desc )
{
  fprintf( f, "   - ");
  if ( desc != (char*)NULL ) fprintf( f, "'%s' ", desc );
  fprintf( f, "[label=%d, order=%d, length=%f, q=%d]\n",
           i->label_branch, i->order, i->length, i->inqueue );

}



void fprintfVoxelTreeComponent( FILE *f, voxelTreeComponent *c )
{
  fprintf( f, "   - type=" );
  switch( c->type ) {
  default : fprintf( f, "default\n" ); break;
  case _INTERN_TREE_UNKNOWN_COMPONENT_ : fprintf( f, "_INTERN_TREE_UNKNOWN_COMPONENT_\n" ); break;
  case _INTERN_TREE_REMOVED_COMPONENT_ : fprintf( f, "_INTERN_TREE_REMOVED_COMPONENT_\n" ); break;
  case _INTERN_TREE_MERGED_COMPONENT_ : fprintf( f, "_INTERN_TREE_MERGED_COMPONENT_\n" ); break;
  case _INTERN_TREE_INNER_EDGE_ : fprintf( f, "_INTERN_TREE_INNER_EDGE_\n" ); break;
  case _INTERN_TREE_END_EDGE_ : fprintf( f, "_INTERN_TREE_END_EDGE_\n" ); break;
  case _INTERN_TREE_EDGE_ : fprintf( f, "_INTERN_TREE_EDGE_\n" ); break;
  case _INTERN_TREE_JUNCTION_ : fprintf( f, "_INTERN_TREE_JUNCTION_\n" ); break;
  }
  fprintf( f, "   - indexFirstJunction=%d\n", c->indexFirstJunction );
  fprintf( f, "   - indexLastJunction =%d\n", c->indexLastJunction );
  fprintf( f, "   - length =%f\n", c->length );
  _fprintfInfoBranchingPoint( f, &(c->infoFirstPoint), "infoFrt" );
  _fprintfInfoBranchingPoint( f, &(c->infoLastPoint), "infoLst" );
  _fprintfInfoBranchingPoint( f, &(c->infoJunction), "infoJct" );
}


void fprintVoxelTree( FILE *f, typeVoxelTree *tree, int writebranches )
{
  int nends;
  int nedges;
  int njunctions;
  int i;

  for ( i=1, nedges=0, njunctions=0, nends=0; i<tree->n_data; i++ ) {
    if ( tree->data[i].type == _INTERN_TREE_END_EDGE_) nends++;
    if ( tree->data[i].type == _INTERN_TREE_INNER_EDGE_) nedges++;
    if ( tree->data[i].type == _INTERN_TREE_JUNCTION_) njunctions++;
  }

  fprintf( f, "===   voxel tree   ===\n" );
  fprintf( f, "=== %d used components, %d inner edges, %d end edges, %d junctions\n",
           tree->n_data-1, nedges, nends, njunctions );
  fprintf( f, "--- inner edges ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_INNER_EDGE_  ) continue;
    fprintEdgeVoxelTree( f, tree, i, writebranches );
  }
  fprintf( f, "--- end edges ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_END_EDGE_  ) continue;
    fprintEdgeVoxelTree( f, tree, i, writebranches );
  }
  fprintf( f, "--- junctions ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_JUNCTION_ ) continue;
    fprintf( f, "    #%3d: %d points, ", i, tree->data[i].n_data );
    _fprintBBComponentVoxelTree( f, &(tree->data[i]) );
    fprintf( f, "(%.1f %.1f %.1f)\n",
             tree->data[i].center.x, tree->data[i].center.y, tree->data[i].center.z );
  }
  fprintf( f, "--- others ---\n" );
  fprintf( f, "    removed components ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_REMOVED_COMPONENT_ ) continue;
    fprintf( f, "    #%3d: %d points\n",
             i, tree->data[i].n_data );
  }
  fprintf( f, "    merged components ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_MERGED_COMPONENT_ ) continue;
    fprintf( f, "    #%3d: %d points\n",
             i, tree->data[i].n_data );
  }
  fprintf( f, "    unknown components ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_UNKNOWN_COMPONENT_ ) continue;
    fprintf( f, "    #%3d: %d points\n",
             i, tree->data[i].n_data );
  }
  fprintf( f, "    other components ---\n" );
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type == _INTERN_TREE_UNKNOWN_COMPONENT_
         || tree->data[i].type == _INTERN_TREE_REMOVED_COMPONENT_
         || tree->data[i].type == _INTERN_TREE_MERGED_COMPONENT_
         || tree->data[i].type == _INTERN_TREE_INNER_EDGE_
         || tree->data[i].type == _INTERN_TREE_END_EDGE_
         || tree->data[i].type == _INTERN_TREE_EDGE_
         || tree->data[i].type == _INTERN_TREE_JUNCTION_ ) continue;
    fprintf( f, "    #%3d: %d points\n",
             i, tree->data[i].n_data );
  }

  if ( writebranches ) {
      fprintf( f, "=== branches\n" );
      for ( i=1; i<tree->n_data; i++ ) {
        if ( tree->data[i].type != _INTERN_TREE_END_EDGE_ )
          continue;
        fprintf( f, "   component %3d: ", i );
        if ( tree->data[i].indexFirstJunction == - 1 ) {
          fprintf( f, " branch %3d, order %2d, length %f",
                   tree->data[i].infoFirstPoint.label_branch,
                   tree->data[i].infoFirstPoint.order,
                   tree->data[i].infoFirstPoint.length );
        }
        if ( tree->data[i].indexLastJunction == - 1 ) {
          fprintf( f, " branch %3d, order %2d, length %f",
                   tree->data[i].infoLastPoint.label_branch,
                   tree->data[i].infoLastPoint.order,
                   tree->data[i].infoLastPoint.length );
        }
        fprintf( f, "\n" );
      }
  }

  fprintf( f, "======================\n" );
}



void fprintfSummaryVoxelTree( FILE *f, typeVoxelTree *c, char *desc )
{
  char *proc = "fprintfSummaryVoxelTree";
  int i;
  int nunknown, funknown, lunknown;
  int nremove, fremove, lremove;
  int nmerge, fmerge, lmerge;
  int ninner, finner, linner;
  int nend, fend, lend;
  int nedge, fedge, ledge;
  int njunction, fjunction, ljunction;

  nunknown = 0;
  funknown =lunknown = -1;
  nremove = 0;
  fremove = lremove = -1;
  nmerge = 0;
  fmerge = lmerge = -1;
  ninner = 0;
  finner = linner = -1;
  nend = 0;
  fend = lend = -1;
  nedge = 0;
  fedge = ledge = -1;
  njunction = 0;
  fjunction = ljunction = -1;

  for ( i=1; i<c->n_data; i++ ) {
    switch( c->data[i].type ) {
    default :
    case _INTERN_TREE_UNKNOWN_COMPONENT_ :
      if ( nunknown == 0 ) {
        funknown = lunknown = i;
      }
      else {
        if ( funknown > i ) funknown = i;
        else if ( lunknown < i ) lunknown = i;
      }
      nunknown ++;
      break;
    case _INTERN_TREE_REMOVED_COMPONENT_ :
      if ( nremove == 0 ) {
        fremove = lremove = i;
      }
      else {
        if ( fremove > i ) fremove = i;
        else if ( lremove < i ) lremove = i;
      }
      nremove ++;
      break;
    case _INTERN_TREE_MERGED_COMPONENT_ :
      if ( nmerge == 0 ) {
        fmerge = lmerge = i;
      }
      else {
        if ( fmerge > i ) fmerge = i;
        else if ( lmerge < i ) lmerge = i;
      }
      nmerge ++;
      break;
    case _INTERN_TREE_INNER_EDGE_ :
      if ( ninner == 0 ) {
        finner = linner = i;
      }
      else {
        if ( finner > i ) finner = i;
        else if ( linner < i ) linner = i;
      }
      ninner ++;
      break;
    case _INTERN_TREE_END_EDGE_ :
      if ( nend == 0 ) {
        fend = lend = i;
      }
      else {
        if ( fend > i ) fend = i;
        else if ( lend < i ) lend = i;
      }
      nend ++;
      break;
    case _INTERN_TREE_EDGE_ :
      if ( nedge == 0 ) {
        fedge = ledge = i;
      }
      else {
        if ( fedge > i ) fedge = i;
        else if ( ledge < i ) ledge = i;
      }
      nedge ++;
      break;
    case _INTERN_TREE_JUNCTION_ :
      if ( njunction == 0 ) {
        fjunction = ljunction = i;
      }
      else {
        if ( fjunction > i ) fjunction = i;
        else if ( ljunction < i ) ljunction = i;
      }
      njunction ++;
      break;
    }
  }

  if ( desc != (char*)NULL )
    fprintf( f, "%s: \n", desc );
  else
    fprintf( f, "%s: \n", proc );

  if ( nunknown > 0 ) {
    fprintf( f, "\t %d unknown labels in [%d %d]\n",
             nunknown, funknown, lunknown );
  }
  if ( nremove > 0 ) {
    fprintf( f, "\t %d removed component labels in [%d %d]\n",
             nremove, fremove, lremove );
  }
  if ( nmerge > 0 ) {
    fprintf( f, "\t %d merged component labels in [%d %d]\n",
             nmerge, fmerge, lmerge );
  }
  if ( ninner > 0 ) {
    fprintf( f, "\t %d inner edge labels in [%d %d]\n",
             ninner, finner, linner );
  }
  if ( nend > 0 ) {
    fprintf( f, "\t %d end edge labels in [%d %d]\n",
             nend, fend, lend );
  }
  if ( nedge > 0 ) {
    fprintf( f, "\t %d edge labels in [%d %d]\n",
             nedge, fedge, ledge );
  }
  if ( njunction > 0 ) {
    fprintf( f, "\t %d junctions labels in [%d %d]\n",
             njunction, fjunction, ljunction );
  }


}






/**************************************************
 *
 *
 *
 **************************************************/



/* count the number of points that have the same label
 */
static int _countSameNeighors( typeComponentImage *treeImage,
                               voxelTreeIntPoint *p )
{
  char *proc = "_countSameNeighors";
  int n = 0;
  int i, j, k;
  int a, b;
  int x = p->x;
  int y = p->y;
  int z = p->z;

  a = (z * treeImage->theDim[1] + y) * treeImage->theDim[0] + x;

#define _COUNTSAMENEIGHBORS( TYPE ) {                         \
  TYPE *theCC = (TYPE*)treeImage->componentBuf;               \
  if ( treeImage->theDim[2] == 1 ) {                          \
    for ( n=0, j=-1; j<=1; j++ ) {                            \
      if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;     \
      for ( i=-1; i<=1; i++ ) {                               \
        if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue;   \
        b = a + j * treeImage->theDim[0] + i;                 \
        if ( theCC[b] == theCC[a] ) n++;                      \
      }                                                       \
    }                                                         \
  }                                                           \
  else {                                                      \
    for ( n=0, k=-1; k<=1; k++ ) {                            \
      if ( z+k<0 || z+k>=treeImage->theDim[2] ) continue;     \
      for ( j=-1; j<=1; j++ ) {                               \
        if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;   \
        for ( i=-1; i<=1; i++ ) {                             \
          if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue; \
          b = a + (k * treeImage->theDim[1] + j) * treeImage->theDim[0] + i; \
          if ( theCC[b] == theCC[a] ) n++;                    \
        }                                                     \
      }                                                       \
    }                                                         \
  }                                                           \
}

  switch( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
    _COUNTSAMENEIGHBORS( u8 );
    return( n-1 );
    break;

  case SSHORT :
    _COUNTSAMENEIGHBORS( s16 );
    return( n-1 );
    break;

  case USHORT :
    _COUNTSAMENEIGHBORS( u16 );
    return( n-1 );
    break;
  }

  return( -1 );
}





/* extrait le ou les voisins d'un point 'center'
 * qui ont le meme label que 'center'
 */
static int _extractEdgeNeighbors( typeComponentImage *treeImage,
                              voxelTreeIntPoint *center,
                              voxelTreeIntPoint *p,
                              voxelTreeIntPoint *l )
{
  char *proc = "_extractEdgeNeighbors";
  int n = 0;
  int i, j, k;
  int a, b, label;
  int x = center->x;
  int y = center->y;
  int z = center->z;

  a = (z * treeImage->theDim[1] + y) * treeImage->theDim[0] + x;

#define _EXTRACTEDGENEIGHBORS( TYPE ) {                       \
  TYPE *theCC = (TYPE*)treeImage->componentBuf;               \
  label = theCC[a];                                           \
  if ( treeImage->theDim[2] == 1 ) {                          \
    for ( n=0, j=-1; j<=1; j++ ) {                            \
      if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;     \
      for ( i=-1; i<=1; i++ ) {                               \
        if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue;   \
        if ( j == 0 && i == 0 ) continue;                     \
        b = a + j * treeImage->theDim[0] + i;                 \
        if ( theCC[b] == label ) {                            \
          if ( n == 0 ) {                                     \
            p->x = x+i;   p->y = y+j;   p->z = z;             \
          }                                                   \
          else if ( n == 1 ) {                                \
            l->x = x+i;   l->y = y+j;   l->z = z;             \
          }                                                   \
          n++;                                                \
        }                                                     \
      }                                                       \
    }                                                         \
  }                                                           \
  else {                                                      \
    for ( n=0, k=-1; k<=1; k++ ) {                            \
      if ( z+k<0 || z+k>=treeImage->theDim[2] ) continue;     \
      for ( j=-1; j<=1; j++ ) {                               \
        if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;   \
        for ( i=-1; i<=1; i++ ) {                             \
          if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue; \
          if ( k == 0 && j == 0 && i == 0 ) continue;         \
          b = a + (k * treeImage->theDim[1] + j) * treeImage->theDim[0] + i; \
          if ( theCC[b] == label ) {                          \
            if ( n == 0 ) {                                   \
              p->x = x+i;   p->y = y+j;   p->z = z + k;       \
            }                                                 \
            else if ( n == 1 ) {                              \
              l->x = x+i;   l->y = y+j;   l->z = z + k;       \
            }                                                 \
            n++;                                              \
          }                                                   \
        }                                                     \
      }                                                       \
    }                                                         \
  }                                                           \
}

  switch( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
    _EXTRACTEDGENEIGHBORS( u8 );
    return( n );
    break;

  case SSHORT :
    _EXTRACTEDGENEIGHBORS( s16 );
    return( n );
    break;

  case USHORT :
    _EXTRACTEDGENEIGHBORS( u16 );
    return( n );
    break;
  }

  return( -1 );
}



/* extrait le ou les voisins d'un point 'center'
 * qui ont un label different de 'center'
 */
static int _extractComponentNeighbors( typeComponentImage *treeImage,
                                       voxelTreeIntPoint *center,
                                       int *l1, int *l2 )
{
  char *proc = "_extractComponentNeighbors";
  int n = 0;
  int i, j, k;
  int a, b, label;
  int x = center->x;
  int y = center->y;
  int z = center->z;

  *l1 = *l2 = -1;
  a = (z * treeImage->theDim[1] + y) * treeImage->theDim[0] + x;

#define _EXTRACTCOMPONENTNEIGHBORS( TYPE ) {                  \
  TYPE *theCC = (TYPE*)treeImage->componentBuf;               \
  label = theCC[a];                                           \
  if ( treeImage->theDim[2] == 1 ) {                          \
    for ( n=0, j=-1; j<=1; j++ ) {                            \
      if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;     \
      for ( i=-1; i<=1; i++ ) {                               \
        if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue;   \
        if ( j == 0 && i == 0 ) continue;                     \
        b = a + j * treeImage->theDim[0] + i;                 \
        if ( theCC[b] > 0 && theCC[b] != label ) {            \
          if ( n == 0 ) {                                     \
            *l1 = theCC[b];                                   \
          }                                                   \
          else if ( n == 1 ) {                                \
            *l2 = theCC[b];                                   \
          }                                                   \
          n++;                                                \
        }                                                     \
      }                                                       \
    }                                                         \
  }                                                           \
  else {                                                      \
    for ( n=0, k=-1; k<=1; k++ ) {                            \
      if ( z+k<0 || z+k>=treeImage->theDim[2] ) continue;     \
      for ( j=-1; j<=1; j++ ) {                               \
        if ( y+j<0 || y+j>=treeImage->theDim[1] ) continue;   \
        for ( i=-1; i<=1; i++ ) {                             \
          if ( x+i<0 || x+i>=treeImage->theDim[0] ) continue; \
          if ( k == 0 && j == 0 && i == 0 ) continue;         \
          b = a + (k * treeImage->theDim[1] + j) * treeImage->theDim[0] + i; \
          if ( theCC[b] > 0 && theCC[b] != label ) {          \
            if ( n == 0 ) {                                   \
              *l1 = theCC[b];                                 \
            }                                                 \
            else if ( n == 1 ) {                              \
              *l2 = theCC[b];                                 \
            }                                                 \
            n++;                                              \
          }                                                   \
        }                                                     \
      }                                                       \
    }                                                         \
  }                                                           \
}

  switch( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
    _EXTRACTCOMPONENTNEIGHBORS( u8 );
    return( n );
    break;

  case SSHORT :
    _EXTRACTCOMPONENTNEIGHBORS( s16 );
    return( n );
    break;

  case USHORT :
    _EXTRACTCOMPONENTNEIGHBORS( u16 );
    return( n );
    break;
  }

  return( -1 );
}




int extractJunctionEdgeVoxelTree( typeComponentImage *treeImage,
                                  voxelTreeComponent *c )
{
  char *proc = "extractJunctionEdgeVoxelTree";
  int label1, label2;
  int ncomponents;


  if ( c->type != _INTERN_TREE_INNER_EDGE_ &&  c->type != _INTERN_TREE_END_EDGE_ && c->type != _INTERN_TREE_EDGE_ ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: component seems to be not an edge\n", proc );
    return( -1 );
  }

  if ( c->n_data <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: component does not have data\n", proc );
    return( -1 );
  }

  /* neighbors of the first point
   */
  ncomponents = _extractComponentNeighbors( treeImage, &(c->data[0]), &label1, &label2 );
  if ( ncomponents == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to extract component neighbor for first point\n", proc );
    return( -1 );
  }



  /* case of the one-point edge
   * get second connecting junction, if any
   */
  if ( c->n_data == 1 ) {
    c->indexFirstJunction = c->indexLastJunction = -1;
    if ( ncomponents == 0 ) {
      ;
    }
    else if ( ncomponents == 1 ) {
      c->indexFirstJunction = label1;
    }
    else if ( ncomponents == 2 ) {
      c->indexFirstJunction = label1;
      c->indexLastJunction = label2;
    }
    else {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird situation for one-point edge\n", proc );
      return( -1 );
    }
    c->type = ( c->indexFirstJunction == -1 || c->indexLastJunction == -1 ) ? _INTERN_TREE_END_EDGE_ : _INTERN_TREE_INNER_EDGE_;
    return( 1 );
  }



  /* general case
   */
  if ( ncomponents == 0 ) {
    c->indexFirstJunction = -1;
  }
  else if ( ncomponents == 1 ) {
    c->indexFirstJunction = label1;
  }
  if ( ncomponents >= 2 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: weird situation for multiple-point edge (first point)\n", proc );
      fprintf( stderr, "\t #data=%d: ncomponents=%d, label1=%d label2=%d\n",
               c->n_data, ncomponents, label1, label2 );
    }
    return( -1 );
  }



  /* neighbors of the last point
   */
  ncomponents = _extractComponentNeighbors( treeImage, &(c->data[c->n_data-1]), &label1, &label2 );
  if ( ncomponents == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to extract component neighbor for last point\n", proc );
    return( -1 );
  }
  if ( ncomponents == 0 ) {
    c->indexLastJunction = -1;
  }
  else if ( ncomponents == 1 ) {
    c->indexLastJunction = label1;
  }
  else {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: weird situation for multiple-point edge (last point)\n", proc );
      fprintf( stderr, "\t #data=%d: ncomponents=%d, label1=%d label2=%d\n",
               c->n_data, ncomponents, label1, label2 );
    }
    return( -1 );
  }

  c->type = ( c->indexFirstJunction == -1 || c->indexLastJunction == -1 ) ? _INTERN_TREE_END_EDGE_ : _INTERN_TREE_INNER_EDGE_;
  return( 1 );
}





/* from a first point, extract the whole edge
 * and extract the connecting junctions
 */
int extractEdgeVoxelTree( typeComponentImage *treeImage,
                         voxelTreeComponent *c, voxelTreeIntPoint *firstPoint )
{
  char *proc = "extractEdgeVoxelTree";
  voxelTreeIntPoint current;
  voxelTreeIntPoint prev;
  voxelTreeIntPoint pt1;
  voxelTreeIntPoint pt2;
  int nneighbors;


  _initVoxelTreeIntPoint( &pt1 );
  _initVoxelTreeIntPoint( &pt2 );

  current = *firstPoint;

  /* add first point
   */
  if ( addPointToVoxelTreeComponent( c, &current ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add first point to component\n", proc );
    return( -1 );
  }

  /* identify edge neighbor for first point
   * at most one neighbor (r=0,1)
   * nneighbors=0 means this is a one-point edge
   */
  nneighbors = _extractEdgeNeighbors( treeImage, &current, &pt1, &pt2 );
  if ( nneighbors == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to extract edge neighbor(s) for first point\n", proc );
    return( -1 );
  }
  else if ( nneighbors > 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: first point has multiple neighbors\n", proc );
    return( -1 );
  }



  /* case of one-point edge
   */
  if ( nneighbors == 0 ) {
    if ( extractJunctionEdgeVoxelTree( treeImage, c ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to extract connecting junctions (one-point edge)\n", proc );
      return( -1 );
    }
    return( 1 );
  }



  prev = current;
  current = pt1;

  do {
    if ( addPointToVoxelTreeComponent( c, &current ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add point to component\n", proc );
      return( -1 );
    }

    nneighbors = _extractEdgeNeighbors( treeImage, &current, &pt1, &pt2 );
    if ( nneighbors == -1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to extract edge neighbor(s) for point\n", proc );
      return( -1 );
    }
    else if ( nneighbors > 2 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: point (%d %d %d) has multiple (%d) neighbors\n",
                 proc, current.x, current.y, current.z, nneighbors );
      return( -1 );
    }

    if ( nneighbors == 2 ) {
      if ( pt1.x == prev.x && pt1.y == prev.y && pt1.z == prev.z ) {
        prev = current;
        current = pt2;
      }
      else if ( pt2.x == prev.x && pt2.y == prev.y && pt2.z == prev.z ) {
        prev = current;
        current = pt1;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: previous point is not found among neighbors\n", proc );
        return( -1 );
      }
    }

  } while ( nneighbors == 2 );


  if ( extractJunctionEdgeVoxelTree( treeImage, c ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to extract connecting junctions (general case)\n", proc );
    return( -1 );
  }

  return( 1 );
}





int voxelTreeJunctionCenter( voxelTreeComponent *c )
{
  char *proc = "voxelTreeJunctionCenter";
  int j, isanchor;

  if ( c->type != _INTERN_TREE_JUNCTION_ ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: component seems to be not a junction\n", proc );
    return( -1 );
  }

  if ( c->n_data <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: component does not have data\n", proc );
    return( -1 );
  }

  _initVoxelTreeFloatPoint( &(c->center) );
  c->center.x = 0.0;
  c->center.y = 0.0;
  c->center.z = 0.0;
  isanchor = c->data[0].anchor;
  for ( j=0; j<c->n_data; j++ ) {
    c->center.x += c->data[j].x;
    c->center.y += c->data[j].y;
    c->center.z += c->data[j].z;
    if ( isanchor < c->data[j].anchor )
      isanchor = c->data[j].anchor;
  }
  c->center.x /= (float)c->n_data;
  c->center.y /= (float)c->n_data;
  c->center.z /= (float)c->n_data;
  c->center.anchor = isanchor;

  return( 1 );
}





void lengthEdgeVoxelTree( typeVoxelTree *tree,
                         int index )
{
  int i, j;
  float px, py, pz, nx, ny, nz;
  float d;
  float l = 0.0;
  voxelTreeComponent *c = &(tree->data[index]);

  c->length = 0.0;

  px = (float)c->data[0].x * tree->voxelSize[0];
  py = (float)c->data[0].y * tree->voxelSize[1];
  pz = (float)c->data[0].z * tree->voxelSize[2];

  if ( c->n_data > 1 ) {
    for ( i=1; i<c->n_data; i++ ) {
      nx = (float)c->data[i].x * tree->voxelSize[0];
      ny = (float)c->data[i].y * tree->voxelSize[1];
      nz = (float)c->data[i].z * tree->voxelSize[2];
      l += sqrt( (nx-px)*(nx-px) + (ny-py)*(ny-py) + (nz-pz)*(nz-pz) );
      px = nx;
      py = ny;
      pz = nz;
    }
    c->length = l;
  }


  /* connect to junction center at edge ends
   */
  j = c->indexFirstJunction;
  if ( j > 0 && tree->data[j].type == _INTERN_TREE_JUNCTION_ ) {
    px = (float)c->data[0].x;
    py = (float)c->data[0].y;
    pz = (float)c->data[0].z;
    d = 0.0;
    d += ( px - tree->data[j].center.x ) * ( px - tree->data[j].center.x ) * tree->voxelSize[0] * tree->voxelSize[0];
    d += ( py - tree->data[j].center.y ) * ( py - tree->data[j].center.y ) * tree->voxelSize[1] * tree->voxelSize[1];
    d += ( pz - tree->data[j].center.z ) * ( pz - tree->data[j].center.z ) * tree->voxelSize[2] * tree->voxelSize[2];
    c->length += sqrt( d );
  }

  j = c->indexLastJunction;
  if ( j > 0 && tree->data[j].type == _INTERN_TREE_JUNCTION_ ) {
    px = (float)c->data[c->n_data-1].x;
    py = (float)c->data[c->n_data-1].y;
    pz = (float)c->data[c->n_data-1].z;
    d = 0.0;
    d += ( px - tree->data[j].center.x ) * ( px - tree->data[j].center.x ) * tree->voxelSize[0] * tree->voxelSize[0];
    d += ( py - tree->data[j].center.y ) * ( py - tree->data[j].center.y ) * tree->voxelSize[1] * tree->voxelSize[1];
    d += ( pz - tree->data[j].center.z ) * ( pz - tree->data[j].center.z ) * tree->voxelSize[2] * tree->voxelSize[2];
    c->length += sqrt( d );
  }

}





/**************************************************
 *
 *
 *
 **************************************************/



/* un 'treeImage' est une image ou les composantes
 * ont ete etiquetees (voir imageToTreeImage())
 * il y a les composantes 'branches' qui vont de
 * treeImage->firstComponentLabel = 1
 * a treeImage->lastComponentLabel = r
 * et les composantes 'junctions' qui vont de
 * tree->firstJunctionLabel = 1+tree->lastComponentLabel;
 * a tree->lastJunctionLabel;
 *
 * le 'typeVoxelTree' est donc compose de tree->lastJunctionLabel+1 composantes
 * (0 n'est pas utilise), chaque composante contient la liste des points
 */
int componentImageToVoxelTree( typeComponentImage *treeImage,
                                   void *theMark,
                                   bufferType typeMark,
                                   typeVoxelTree *tree )
{
  char *proc = "componentImageToVoxelTree";
  int x, y, z, i, j;
  int r;
  voxelTreeIntPoint pt;
  voxelTreeIntPoint *ptr;

  /* la numerotation commence avec les composantes branches
   * puis viennent les composantes jonctions
   * tree->firstJunctionLabel and tree->lastJunctionLabel are set to -1
   * if there is no junctions
   */

  initVoxelTree( tree );
  if ( _allocVoxelTree( tree, treeImage->n_component ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  tree->voxelSize[0] = treeImage->voxelSize[0];
  tree->voxelSize[1] = treeImage->voxelSize[1];
  tree->voxelSize[2] = treeImage->voxelSize[2];

  /* pour les jonctions:
   * - on ajoute le point a la composante jonction
   *   on met a jour la bounding box
   *
   * pour les aretes
   * - on attend le cas ou on a 0 ou 1 voisin du meme type
   *   on extrait la composante
   */
#define _TREEIMAGETOVOXELTREE( TYPE ) {                                      \
  TYPE *theCC = (TYPE*)treeImage->componentBuf;                              \
  for ( i=0, z=0; z<treeImage->theDim[2]; z++ )                              \
  for ( y=0; y<treeImage->theDim[1]; y++ )                                   \
  for ( x=0; x<treeImage->theDim[0]; x++, i++ ) {                            \
    if ( theCC[i] == 0 ) continue;                                           \
    pt.x = x;   pt.y = y;   pt.z = z;                                        \
    if ( treeImage->component[ theCC[i] ].type == _JUNCTION_COMPONENT_ ) {   \
      tree->data[theCC[i]].type = _INTERN_TREE_JUNCTION_;                    \
      if ( addPointToVoxelTreeComponent( &(tree->data[theCC[i]]), &pt ) != 1 ) { \
        freeVoxelTree( tree );                                               \
        return( -1 );                                                        \
      }                                                                      \
    }                                                                        \
    else if ( treeImage->component[ theCC[i] ].type == _EDGE_COMPONENT_ ) {  \
      if ( tree->data[theCC[i]].n_data > 0 ) continue;                       \
      r = _countSameNeighors( treeImage, &pt );                              \
      if ( r == -1 ) {                                                       \
        freeVoxelTree( tree );                                               \
        if ( _verbose_ )                                                     \
          fprintf( stderr, "%s: error when counting neigbors of (%d,%d,%d)\n", proc, x, y, z ); \
        return( -1 );                                                        \
      }                                                                      \
      if ( r == 1 || r == 0 ) {                                              \
        tree->data[theCC[i]].type = _INTERN_TREE_EDGE_;                      \
        if ( extractEdgeVoxelTree( treeImage, &(tree->data[theCC[i]]), &pt ) != 1 ) { \
            freeVoxelTree( tree );                                           \
            if ( _verbose_ )                                                 \
              fprintf( stderr, "%s: error when extracting edge #%d from (%d,%d,%d)\n", proc, theCC[i], x, y, z ); \
            return( -1 );                                                    \
        }                                                                    \
      }                                                                      \
    }                                                                        \
    else {                                                                   \
      freeVoxelTree( tree );                                                 \
      if ( _verbose_ )                                                       \
        fprintf( stderr, "%s: weird case, %d is not a valid label\n", proc, theCC[i] ); \
      return( -1 );                                                          \
    }                                                                        \
  }                                                                          \
}

  /* build tree
   */

  _initVoxelTreeIntPoint( &pt );

  switch( treeImage->componentBufferType ) {
  default :
    freeVoxelTree( tree );
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
      _TREEIMAGETOVOXELTREE( u8 );
      break;

  case SSHORT :
      _TREEIMAGETOVOXELTREE( s16 );
      break;

  case USHORT :
      _TREEIMAGETOVOXELTREE( u16 );
      break;
  }

  /* label point
   * this is useful to differentiate the main axon
   * from the built outgoing branches
   * default label is 0 (from _initVoxelTreeIntPoint())
   * label of main axon will be 1
   */
#define _MARKVoxelTree( TYPE ) {               \
  TYPE *bufMark = (TYPE*)theMark;              \
  for ( j=0; j<tree->n_data; j++ ) {           \
    for ( i=0; i<tree->data[j].n_data; i++ ) { \
      ptr = &(tree->data[j].data[i]);          \
      r = (ptr->z * treeImage->theDim[1] + ptr->y) * treeImage->theDim[0] + ptr->x; \
      ptr->anchor = ( bufMark[r] > 0 ) ? 1 : 0; \
    }                                          \
  }                                            \
}

  if ( theMark != (void*)NULL && typeMark != TYPE_UNKNOWN ) {
    switch( typeMark ) {
    default :
      freeVoxelTree( tree );
      if ( _verbose_ )
        fprintf( stderr, "%s: such mark image type not handled yet\n", proc );
      return( -1 );

    case UCHAR :
      _MARKVoxelTree( u8 );
      break;
    }
  }

  /* compute barycenter for junctions
   */
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_JUNCTION_ ) continue;
    if ( voxelTreeJunctionCenter( &(tree->data[i]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to compute barycenter for junction #%d\n", proc, i );
      return( -1 );
    }
  }


  /* test if branch is an end-branch
   * compute edge length
   */
  for ( i=1; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && tree->data[i].type != _INTERN_TREE_END_EDGE_
         && tree->data[i].type != _INTERN_TREE_EDGE_ ) continue;
    if ( tree->data[i].indexFirstJunction == -1 || tree->data[i].indexLastJunction == -1 )
      tree->data[i].type = _INTERN_TREE_END_EDGE_;
    else
      tree->data[i].type = _INTERN_TREE_INNER_EDGE_;
    lengthEdgeVoxelTree( tree, i );
  }

  return( 1 );
}





int labelBranchImageToBranches( typeComponentImage *treeImage,
                                      typeVoxelTree *tree )
{
  char *proc = "labelBranchImageToBranches";
  float max;
  int nbranches;
  int x, y, z, i;
  int r;
  voxelTreeIntPoint pt;

  int n;
  voxelTreeIntPoint *data;
  int tmpJunction;
  infoBranchingPoint tmpInfo;



  if ( maxValue( treeImage->componentBuf, treeImage->componentBufferType,
                 treeImage->theDim, &max ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute number of branches\n", proc );
    return( -1 );
  }

  nbranches = (int)(max + 0.5);

  if ( _debug_ )
    fprintf( stderr, "%s: found %d labels\n", proc, nbranches );

  initVoxelTree( tree );
  if ( _allocVoxelTree( tree, nbranches+1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  tree->voxelSize[0] = treeImage->voxelSize[0];
  tree->voxelSize[1] = treeImage->voxelSize[1];
  tree->voxelSize[2] = treeImage->voxelSize[2];

#define _LABELBRANCHIMAGETOVOXELTREE( TYPE ) {                               \
  TYPE *theCC = (TYPE*)treeImage->componentBuf;                              \
  for ( i=0, z=0; z<treeImage->theDim[2]; z++ )                              \
  for ( y=0; y<treeImage->theDim[1]; y++ )                                   \
  for ( x=0; x<treeImage->theDim[0]; x++, i++ ) {                            \
    if ( theCC[i] == 0 || theCC[i] == 1 ) continue;                          \
    pt.x = x;   pt.y = y;   pt.z = z;                                        \
    r = _countSameNeighors( treeImage, &pt );                                \
    if ( r == -1 ) {                                                         \
      freeVoxelTree( tree );                                                 \
      if ( _verbose_ )                                                       \
        fprintf( stderr, "%s: error when counting neigbors of (%d,%d,%d)\n", proc, x, y, z ); \
      return( -1 );                                                          \
    }                                                                        \
    if ( r == 2 ) continue;                                                  \
    if ( r == 1 || r == 0 ) {                                                \
      if ( tree->data[theCC[i]].n_data > 0 ) {                               \
        if ( 0 && _debug_ )                                                  \
          fprintf( stderr, "%s: branch #%d already extracted\n", proc, theCC[i] ); \
        continue;                                                            \
      }                                                                      \
      tree->data[theCC[i]].type = _INTERN_TREE_EDGE_;                        \
      if ( extractEdgeVoxelTree( treeImage, &(tree->data[theCC[i]]), &pt ) != 1 ) { \
        freeVoxelTree( tree );                                               \
        if ( _verbose_ )                                                     \
          fprintf( stderr, "%s: error when extracting edge #%d from (%d,%d,%d)\n", proc, theCC[i], x, y, z ); \
        return( -1 );                                                        \
      }                                                                      \
    }                                                                        \
    else {                                                                   \
      freeVoxelTree( tree );                                                 \
      if ( _verbose_ )                                                       \
        fprintf( stderr, "%s: weird number of neighbors (%d)\n", proc, r );  \
      return( -1 );                                                          \
    }                                                                        \
  }                                                                          \
}


  /* build tree
   */

  _initVoxelTreeIntPoint( &pt );

  switch( treeImage->componentBufferType ) {
  default :
    freeVoxelTree( tree );
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case UCHAR :
      _LABELBRANCHIMAGETOVOXELTREE( u8 );
      break;

  case SSHORT :
      _LABELBRANCHIMAGETOVOXELTREE( s16 );
      break;

  case USHORT :
      _LABELBRANCHIMAGETOVOXELTREE( u16 );
      break;
  }

  /* each component is a branch connected to the '1' subtree
   * for further easier computation, begin each component at '1'
   */
  for ( i=2; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && tree->data[i].type != _INTERN_TREE_END_EDGE_
         && tree->data[i].type != _INTERN_TREE_EDGE_ ) continue;
    if ( tree->data[i].indexFirstJunction == 1 ) continue;

    data = (voxelTreeIntPoint *)vtmalloc( tree->data[i].n_data * sizeof(voxelTreeIntPoint),
                                          "data", proc );
    if ( data == (voxelTreeIntPoint *)NULL ) {
      freeVoxelTree( tree );
      if ( _verbose_ )
        fprintf( stderr, "%s: allocation errorn", proc );
      return( -1 );
    }

    for ( n=0; n<tree->data[i].n_data; n++ )
      data[n] = tree->data[i].data[ tree->data[i].n_data - 1 - n ];
    vtfree( tree->data[i].data );
    tree->data[i].data = data;
    tree->data[i].n_allocated_data = tree->data[i].n_data;

    tmpJunction = tree->data[i].indexFirstJunction;
    tree->data[i].indexFirstJunction = tree->data[i].indexLastJunction;
    tree->data[i].indexLastJunction = tmpJunction;

    tmpInfo = tree->data[i].infoFirstPoint;
    tree->data[i].infoFirstPoint = tree->data[i].infoLastPoint;
    tree->data[i].infoLastPoint = tmpInfo;
  }

  /* compute length estimation
   * the length to the connected junction ismising
   */
  for ( i=2; i<tree->n_data; i++ ) {
    if ( tree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && tree->data[i].type != _INTERN_TREE_END_EDGE_
         && tree->data[i].type != _INTERN_TREE_EDGE_ ) continue;
    lengthEdgeVoxelTree( tree, i );
  }

  if ( _debug_ ) {
    fprintf( stderr, "%s: after branch extraction\n", proc );
    fprintVoxelTree( stderr, tree, 1 );
  }

  return( 1 );
}




/**************************************************
 *
 *
 *
 **************************************************/

int countEndEdgeInImage( void* theBuf, bufferType theType, int *theDim )
{
  char *proc = "countEndEdgeInImage";
  typeComponentImage treeImage;
  typeVoxelTree iTree;
  int i, e;

  /* build tree image
   */

  initComponentImage( &treeImage );
  if ( allocComponentImage( &treeImage, USHORT, theDim ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error for tree image\n", proc );
    return( -1 );
  }

  if ( imageToComponentImage( theBuf, theType, &treeImage ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when labeling components in tree\n", proc );
    return( -1 );
  }

  /* build internal tree
   */

  initVoxelTree( &iTree );
  if ( componentImageToVoxelTree( &treeImage, (void*)NULL, TYPE_UNKNOWN, &iTree ) != 1 ) {
    freeVoxelTree( &iTree );
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building internal tree\n", proc );
    return( -1 );
  }

  for ( e=0, i=0; i<iTree.n_data; i++ ) {
    if ( iTree.data[i].type == _INTERN_TREE_END_EDGE_ ) {
        e++;
    }
  }

  freeVoxelTree( &iTree );
  freeComponentImage( &treeImage );

  return( e );
}










/**************************************************
 *
 *
 *
 **************************************************/

int voxelTreeToImage( typeVoxelTree *iTree,
                      void *resBuf, bufferType resType, int *resDim,
                      enumVoxelTreeAttribute attribute )
{
  char *proc = "";
  int i, n;
  size_t ind;
  size_t dimx = resDim[0];
  size_t dimy = resDim[1];
  size_t v = dimx*dimy*(size_t)resDim[2];
  voxelTreeComponent *c;



#define _ZEROFILL_( TYPE ) {          \
  TYPE *buf = (TYPE*)resBuf;         \
  for ( ind=0; ind<v; ind++, buf++ ) \
    *buf = 0;                        \
}

  switch( resType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _ZEROFILL_( u8 );
    break;
  case USHORT :
    _ZEROFILL_( u16 );
    break;
  }



#define _FILL_NO_ATTRIBUTE_( TYPE, VAL ) { \
  TYPE *buf = (TYPE*)resBuf;               \
  for ( i=1; i<iTree->n_data; i++ ) {      \
    c = &(iTree->data[i]);                 \
    switch( c->type ) {                    \
    default :                              \
      break;                               \
    case _INTERN_TREE_INNER_EDGE_ :        \
    case _INTERN_TREE_END_EDGE_ :          \
    case _INTERN_TREE_EDGE_ :              \
    case _INTERN_TREE_JUNCTION_ :          \
      for ( n=0; n<c->n_data; n++ ) {      \
        ind = ((size_t)c->data[n].z * dimy+ (size_t)c->data[n].y) * dimx + (size_t)c->data[n].x; \
        buf[ind] = VAL;                    \
      }                                    \
    }                                      \
  }                                        \
}



#define _FILL_COMPONENT_LABEL_( TYPE ) {   \
  TYPE *buf = (TYPE*)resBuf;               \
  for ( i=1; i<iTree->n_data; i++ ) {      \
    c = &(iTree->data[i]);                 \
    switch( c->type ) {                    \
    default :                              \
      break;                               \
    case _INTERN_TREE_INNER_EDGE_ :        \
    case _INTERN_TREE_END_EDGE_ :          \
    case _INTERN_TREE_EDGE_ :              \
    case _INTERN_TREE_JUNCTION_ :          \
      for ( n=0; n<c->n_data; n++ ) {      \
        ind = ((size_t)c->data[n].z * dimy+ (size_t)c->data[n].y) * dimx + (size_t)c->data[n].x; \
        buf[ind] = i;                      \
      }                                    \
    }                                      \
  }                                        \
}



#define _FILL_BRANCH_LABEL_( TYPE ) {      \
  TYPE *buf = (TYPE*)resBuf;               \
  for ( i=1; i<iTree->n_data; i++ ) {      \
    c = &(iTree->data[i]);                 \
    switch( c->type ) {                    \
    default :                              \
      break;                               \
    case _INTERN_TREE_INNER_EDGE_ :        \
    case _INTERN_TREE_END_EDGE_ :          \
    case _INTERN_TREE_EDGE_ :              \
    case _INTERN_TREE_JUNCTION_ :          \
      for ( n=0; n<c->n_data; n++ ) {      \
        ind = ((size_t)c->data[n].z * dimy+ (size_t)c->data[n].y) * dimx + (size_t)c->data[n].x; \
        buf[ind] = c->data[n].label_branch; \
      }                                    \
    }                                      \
  }                                        \
}



#define _FILL_BRANCH_ORDER_( TYPE ) {      \
  TYPE *buf = (TYPE*)resBuf;               \
  for ( i=1; i<iTree->n_data; i++ ) {      \
    c = &(iTree->data[i]);                 \
    switch( c->type ) {                    \
    default :                              \
      break;                               \
    case _INTERN_TREE_INNER_EDGE_ :        \
    case _INTERN_TREE_END_EDGE_ :          \
    case _INTERN_TREE_EDGE_ :              \
    case _INTERN_TREE_JUNCTION_ :          \
      for ( n=0; n<c->n_data; n++ ) {      \
        ind = ((size_t)c->data[n].z * dimy+ (size_t)c->data[n].y) * dimx + (size_t)c->data[n].x; \
        buf[ind] = c->data[n].order_branch; \
      }                                    \
    }                                      \
  }                                        \
}


  switch( attribute ) {
  default :
  case _NO_ATTRIBUTE_ :
  case _POINT_NEIGHBORS_ :
    switch( resType ) {
    default :
      if ( _verbose_ )
       fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _FILL_NO_ATTRIBUTE_( u8, 255 );
      break;
    case USHORT :
      _FILL_NO_ATTRIBUTE_( u16, 65535 );
      break;
    }
    break;
  case _COMPONENT_LABEL_ :
    switch( resType ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _FILL_COMPONENT_LABEL_( u8 );
      break;
    case USHORT :
      _FILL_COMPONENT_LABEL_( u16 );
      break;
    }
    break;
  case _BRANCH_LABEL_ :
    switch( resType ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _FILL_BRANCH_LABEL_( u8 );
      break;
    case USHORT :
      _FILL_BRANCH_LABEL_( u16 );
      break;
    }
    break;
  case _BRANCH_ORDER_ :
    switch( resType ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such image type not handled yet\n", proc );
      return( -1 );
    case UCHAR :
      _FILL_BRANCH_ORDER_( u8 );
      break;
    case USHORT :
      _FILL_BRANCH_ORDER_( u16 );
      break;
    }
    break;



  }


  if ( attribute == _POINT_NEIGHBORS_ ) {
    if ( countNeighborsInImage( resBuf, resType, resBuf, resType, resDim ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when counting neighbors in image tree\n", proc );
      return( -1 );
    }
  }

  return( 1 );
}











/**************************************************
 *
 *
 *
 **************************************************/

int writeComponentLengthFromVoxelTree( char *name, typeVoxelTree *iTree, int fullwrite )
{
  char *proc = "writeComponentLengthFromVoxelTree";
  int i;
  FILE *f;
  voxelTreeComponent *c;

  if ( name == (char*)NULL || name[0] == '\0' || name[0] == '>' )
    return( 1 );

  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
    return( -1 );
  }

  if ( fullwrite ) {
    fprintf( f, "# label length\n" );
  }

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch( c->type ) {
    default :
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      if ( fullwrite )
        fprintf( f, "%3d ", i );
      fprintf( f, "%f\n", c->length );
    }
  }

  fclose( f );

  return( 1 );
}



int writeBranchLengthFromVoxelTree( char *name, typeVoxelTree *iTree, int fullwrite )
{
  char *proc = "writeBranchLengthFromVoxelTree";
  int i;
  FILE *f;
  voxelTreeComponent *c;

  if ( name == (char*)NULL || name[0] == '\0' || name[0] == '>' )
    return( 1 );

  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
    return( -1 );
  }

  if ( fullwrite ) {
    fprintf( f, "# label order length\n" );
  }

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch( c->type ) {
    default :
      break;
    case _INTERN_TREE_END_EDGE_ :
      if ( c->infoFirstPoint.label_branch == 1 && c->infoLastPoint.label_branch == 1 )
        continue;
      if ( c->infoFirstPoint.length > c->infoLastPoint.length ) {
        if ( fullwrite )
          fprintf( f, "%3d %2d ", c->infoFirstPoint.label_branch, c->infoFirstPoint.order );
        fprintf( f, "%f\n", c->infoFirstPoint.length );
      }
      if ( c->infoLastPoint.length > c->infoFirstPoint.length ) {
        if ( fullwrite )
          fprintf( f, "%3d %2d ", c->infoLastPoint.label_branch, c->infoLastPoint.order );
        fprintf( f, "%f\n", c->infoLastPoint.length );
      }
    }
  }

  fclose( f );

  return( 1 );
}










