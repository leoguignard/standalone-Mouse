/****************************************************
 * topological-generic-tree.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2015, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 10:28:22 CEST
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

#include <vtmalloc.h>


#include <topological-generic-tree.h>








/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalGenericTree( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalGenericTree(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalGenericTree(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalGenericTree( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalGenericTree(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalGenericTree(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}








/**************************************************
 *
 *
 *
 **************************************************/



static void initTreeFloatPoint( typeTreeFloatPoint *p )
{
  p->type = _TREE_UNKNOWN_POINT_;
  p->x = -1.0;
  p->y = -1.0;
  p->z = -1.0;
  p->anchor = 0;
  p->label_branch = 0;
  p->label_component = 0;
  p->color = 0;
}










/**************************************************
 *
 *
 *
 **************************************************/



static void _initTreeEdge( typeTreeEdge *e )
{
  e->type = _UNDEFINED_EDGE_;
  e->label = 0;
  e->n_data = 0;
  e->n_allocated_data = 0;
  e->data = (int*)NULL;
}



static void _freeTreeEdge( typeTreeEdge *e )
{
  if ( e->data != (int*)NULL )
    vtfree( e->data );
  _initTreeEdge( e );
}



void initTree( typeTree *t )
{
  t->n_points = 0;
  t->n_allocated_points = 0;
  t->points = (typeTreeFloatPoint*)NULL;
  t->n_edges = 0;
  t->n_allocated_edges = 0;
  t->edges = (typeTreeEdge*)NULL;
}



void freeTree( typeTree *t )
{
  int i;
  if ( t->points != (typeTreeFloatPoint*)NULL )
    vtfree( t->points );
  if ( t->n_edges > 0 && t->edges != (typeTreeEdge*)NULL ) {
    for ( i=0; i<t->n_edges; i++ )
      _freeTreeEdge( &(t->edges[i]) );
    vtfree( t->edges );
  }
  initTree( t );
}





/**************************************************
 *
 *
 *
 **************************************************/

int voxelTreeToTree( typeVoxelTree *iTree, typeTree *tree )
{
  char *proc = "voxelTreeToTree";
  int npointsedges;
  int npoints;
  int nends;
  int nedges;
  int njunctions;
  int nelements;
  int i, j, k, v, e;

  /* at this point, points of the intern Tree can have a mark
   * (label field) that is
   * - either 0 (default value)
   * - or a value extracted from the buffer theMark (if not NULL)
   * for junction points, it is the maximal value of marks (since a
   * junction may be made of several poins)
   *
   * this mark can be used for coloring the tree
   */

  /* extract information
   */

  for ( i=1, nends=0, nedges=0, njunctions=0; i<iTree->n_data; i++ ) {
    if ( iTree->data[i].type == _INTERN_TREE_END_EDGE_ ) nends++;
    if ( iTree->data[i].type == _INTERN_TREE_INNER_EDGE_ ) nedges++;
    if ( iTree->data[i].type == _INTERN_TREE_JUNCTION_ ) njunctions++;
  }

  if ( _debug_  ) {
    fprintf( stderr, "%s: found %d end edges, %d inner edges, %d junctions\n",
             proc, nends, nedges, njunctions );
    fprintVoxelTree( stderr, iTree, 0 );
  }


  /* recall that iTree->data[0] remains unknown
   */
  /* this stands if there was no pruning
   *
  if ( nedges+nends+njunctions != iTree->n_data-1 ) {
     if ( _verbose_ )
       fprintf( stderr, "%s: unknown components in internal tree\n", proc );
     return( -1 );
  }
  */
  nedges += nends;

  /* count the number of points
   * recall that junction components will be represented by a single point
   */
  for ( i=1, npointsedges=0; i<iTree->n_data; i++ ) {
    if ( iTree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && iTree->data[i].type != _INTERN_TREE_END_EDGE_ ) continue;
    npointsedges += iTree->data[i].n_data;
  }
  npoints = npointsedges + njunctions;

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: found %d points, %d edges, %d junctions\n",
             proc, npoints, nedges, njunctions );
  }


  /* allocate tree
   */
  tree->points = (typeTreeFloatPoint*)vtmalloc( npoints*sizeof(typeTreeFloatPoint),
                                                "tree->points", proc );
  if ( tree->points == (typeTreeFloatPoint*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: tree points allocation failed\n", proc );
    return( -1 );
  }
  tree->n_points = npoints;
  tree->n_allocated_points = npoints;
  for ( i=0; i<tree->n_points; i++ )
    initTreeFloatPoint( &(tree->points[i]) );

  tree->edges = (typeTreeEdge*)vtmalloc( nedges*sizeof(typeTreeEdge),
                                         "tree->edges", proc );
  if ( tree->edges == (typeTreeEdge*)NULL ) {
    freeTree( tree );
    if ( _verbose_ )
      fprintf( stderr, "%s: tree points allocation failed\n", proc );
    return( -1 );
  }

  /* branches
   */
  tree->n_edges = nedges;
  tree->n_allocated_edges = nedges;
  for ( i=0; i<tree->n_edges; i++ ) {
    _initTreeEdge( &(tree->edges[i]) );
  }

  /* fill tree from internal tree
   * 1. fill all the points
   *    label is either 0 (default value)
   *    or is the value of the buffer theMark (if not NULL)
   */
  for ( i=1, e=1, v=0; i<iTree->n_data; i++ ) {
    switch( iTree->data[i].type ) {
    default :
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
      iTree->data[i].pointIndexOffset = v;
      for ( j=0; j<iTree->data[i].n_data; j++, v++ ) {
        tree->points[v].type = _TREE_EDGE_POINT_;
        tree->points[v].x = iTree->data[i].data[j].x * iTree->voxelSize[0];
        tree->points[v].y = iTree->data[i].data[j].y * iTree->voxelSize[1];
        tree->points[v].z = iTree->data[i].data[j].z * iTree->voxelSize[2];
        tree->points[v].anchor = iTree->data[i].data[j].anchor;
        tree->points[v].order_branch = iTree->data[i].data[j].order_branch;
        tree->points[v].label_branch = iTree->data[i].data[j].label_branch;
        tree->points[v].label_component = e;
        e ++;
      }
      break;
    case _INTERN_TREE_JUNCTION_ :
      iTree->data[i].pointIndexOffset = v;
      tree->points[v].type = _TREE_JUNCTION_POINT_;
      tree->points[v].x = iTree->data[i].center.x * iTree->voxelSize[0];
      tree->points[v].y = iTree->data[i].center.y * iTree->voxelSize[1];
      tree->points[v].z = iTree->data[i].center.z * iTree->voxelSize[2];
      tree->points[v].anchor = iTree->data[i].center.anchor;
      tree->points[v].order_branch = iTree->data[i].center.order_branch;
      tree->points[v].label_branch = iTree->data[i].center.label_branch;
      tree->points[v].label_component = 0;
      v++;
      break;
    }
  }

  /* fill tree from internal tree
   * 2. build edges
   *    junctions are added at the end if required
   *
   */
  for ( i=1, e=0; i<iTree->n_data; i++ ) {
    /* skip junction component
     */
    if ( iTree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && iTree->data[i].type != _INTERN_TREE_END_EDGE_ ) continue;
    /* go for a regular branch
     */
    nelements = iTree->data[i].n_data;
    if ( iTree->data[i].indexFirstJunction > 0 ) nelements++;
    if ( iTree->data[i].indexLastJunction > 0 ) nelements++;
    tree->edges[e].data = (int*)vtmalloc( nelements*sizeof(int), "tree->edges[e].data", proc );
    if ( tree->edges[e].data == (int*)NULL ) {
      freeTree( tree );
      if ( _verbose_ )
        fprintf( stderr, "%s: tree edge allocation failed for edge #%d\n", proc, e );
      return( -1 );
    }
    tree->edges[e].n_data = nelements;
    tree->edges[e].n_allocated_data = nelements;
    tree->edges[e].label = i;
    /* fill the data with the point indexes
     * first point is a junction if any
     * middle points are the edge points
     * last point is a junction if any
     */
    k = 0;
    if ( iTree->data[i].indexFirstJunction > 0 ) {
      tree->edges[e].data[k] = iTree->data[ iTree->data[i].indexFirstJunction ].pointIndexOffset;
      k++;
    }
    /* iTree->data[i].pointIndexOffset est l'index du premier point
     * de la composante dans les points tree->points
     * pour les jonctions, il s'agit d'un seul point
     */
    for ( j=0; j<iTree->data[i].n_data; j++, k++ ) {
      tree->edges[e].data[k] = iTree->data[i].pointIndexOffset + j;
      if ( (j == 0 && iTree->data[i].indexFirstJunction <= 0)
           || (j == iTree->data[i].n_data-1 && iTree->data[i].indexLastJunction <= 0) ) {
          tree->points[ tree->edges[e].data[k] ].type = _TREE_END_POINT_;
      }
    }
    if ( iTree->data[i].indexLastJunction > 0 ) {
      tree->edges[e].data[k] = iTree->data[ iTree->data[i].indexLastJunction ].pointIndexOffset;
      k++;
    }
    /* keep trace of the branch type
     */
    if ( iTree->data[i].indexFirstJunction > 0 && iTree->data[i].indexLastJunction > 0 ) {
      tree->edges[e].type = _INNER_EDGE_;
    }
    else {
      tree->edges[e].type = _END_EDGE_;
    }

    e++;
  }

  return( 1 );
}





int treeImageToTree( typeComponentImage *treeImage,
                     void *theMark, bufferType typeMark,
                     typeTree *tree )
{
  char *proc = "treeImageToTree";
  typeVoxelTree iTree;


  /* build internal tree
   */
  initVoxelTree( &iTree );
  if ( componentImageToVoxelTree( treeImage, theMark, typeMark, &iTree ) != 1 ) {
    freeVoxelTree( &iTree );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building internal tree\n", proc );
    return( -1 );
  }

  if ( voxelTreeToTree( &iTree, tree ) != 1 ) {
    freeTree( tree );
    freeVoxelTree( &iTree );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building tree\n", proc );
    return( -1 );
  }


  freeVoxelTree( &iTree );
  return( 1 );
}








/**************************************************
 *
 *
 *
 **************************************************/



void initTreeOutputParam( typeTreeOutputParam *p )
{
  p->color = _NO_COLOR_;
  p->filetype = _SINGLE_FILE_;
  p->colortable = _USER_COLORTABLE_;
  p->table[0][0] = 1.0;
  p->table[0][1] = 1.0;
  p->table[0][2] = 1.0;
  p->table[1][0] = 1.0;
  p->table[1][1] = 0.0;
  p->table[1][2] = 0.0;
}





void fprintfTreeOutputParam( FILE *f, typeTreeOutputParam *p, char *str )
{
  int i;

  fprintf( f, "- " );
  if ( str != (char*)NULL ) fprintf( f, "%s->", str );

  fprintf( f, "- " );
  if ( str != (char*)NULL ) fprintf( f, "%s->", str );
  fprintf( f, "color = " );
  switch( p->color ) {
  default :
    fprintf( f, "not defined, how embarassing\n" ); break;
  case _NO_COLOR_ :
    fprintf( f, "_NO_COLOR_\n" ); break;
  case _COLOR_ANCHOR_ :
    fprintf( f, "_COLOR_ANCHOR_\n" ); break;
  case _COLOR_BRANCH_ORDER_ :
    fprintf( f, "_COLOR_BRANCH_ORDER_\n" ); break;
  case _COLOR_BRANCH_LABEL_ :
    fprintf( f, "_COLOR_BRANCH_LABEL_\n" ); break;
  case _COLOR_COMPONENT_LABEL_ :
    fprintf( f, "_COLOR_COMPONENT_LABEL_\n" ); break;
  }

  fprintf( f, "- " );
  if ( str != (char*)NULL ) fprintf( f, "%s->", str );
  fprintf( f, "filetype = " );
  switch( p->filetype ) {
  default :
    fprintf( f, "not defined, how embarassing\n" ); break;
  case _MULTIPLE_FILES_ :
    fprintf( f, "_MULTIPLE_FILES_\n" ); break;
  case _SINGLE_FILE_ :
    fprintf( f, "_SINGLE_FILE_\n" ); break;
  }

  fprintf( f, "- " );
  if ( str != (char*)NULL ) fprintf( f, "%s->", str );
  fprintf( f, "colortable = " );
  switch( p->colortable ) {
  default :
    fprintf( f, "not defined, how embarassing\n" ); break;
  case _DEFAULT_COLORTABLE_ :
    fprintf( f, "_DEFAULT_COLORTABLE_\n" ); break;
  case _USER_COLORTABLE_ :
    fprintf( f, "_USER_COLORTABLE_\n" ); break;
  }

  for ( i=0; i<2; i++ ) {
    fprintf( f, "- " );
    if ( str != (char*)NULL ) fprintf( f, "%s->", str );
    fprintf( f, "table[%d] = ", i );
    fprintf( f, "[%f %f %f]\n", p->table[i][0], p->table[i][1], p->table[i][2] );
  }
}










/**************************************************
 *
 *
 *
 **************************************************/


int _colorTree( typeTree *tree, enumTreePointColor color )
{
  char *proc = "_colorTree";
  int i;
  typeTreeFloatPoint *pt;

  switch( color ) {
  default :
  case _NO_COLOR_ :
    for ( i=0; i<tree->n_points; i++ ) {
      tree->points[i].color = 0;
    }
    return( 1 );

  case _COLOR_ANCHOR_ :
    for ( i=0; i<tree->n_points; i++ ) {
      tree->points[i].color = ( tree->points[i].anchor ) ? 0 : 1;
    }
    return( 1 );

  case _COLOR_BRANCH_ORDER_ :
    for ( i=0; i<tree->n_points; i++ ) {
        tree->points[i].color = tree->points[i].order_branch;
    }
    return( 1 );

  case _COLOR_BRANCH_LABEL_ :
    for ( i=0; i<tree->n_points; i++ ) {
        tree->points[i].color = tree->points[i].label_branch;
    }
    return( 1 );

  case _COLOR_COMPONENT_LABEL_ :
    for ( i=0; i<tree->n_points; i++ ) {
      pt = &(tree->points[ i ]);
      if ( pt->type == _TREE_JUNCTION_POINT_ ) continue;
      pt->color = pt->label_component;
    }
    return( 1 );
  }

  if ( _verbose_ ) {
    fprintf( stderr, "%s: weird, this should not be reached\n", proc );
  }

  return( -1 );
}







/**************************************************
 *
 *
 *
 **************************************************/




/* write one file per edge
 */
int _writeVTK_singleEdges( char *format, typeTree *tree, typeTreeOutputParam *par, char *desc )
{
  char *proc = "_writeVTK_singleEdges";
  int i, j, l;
  typeTreeFloatPoint *pt;
  FILE *f;
  char name[STRINGLENGTH];

  for ( i=0; i<tree->n_edges; i++ ) {

    sprintf( name, format, i+1 );
    f = fopen( name, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
      return( -1 );
    }

    fprintf( f, "# vtk DataFile Version 2.0\n" );

    if ( desc != (char*)NULL )
        fprintf( f, "%s", desc );
    fprintf( f, "\n" );

    fprintf( f, "ASCII\n" );
    fprintf( f, "DATASET POLYDATA\n" );

    fprintf( f, "POINTS %d float\n", tree->edges[i].n_data );
    for ( j=0; j<tree->edges[i].n_data; j++ ) {
      l = tree->edges[i].data[j];
      fprintf( f, "%f %f %f\n", tree->points[l].x, tree->points[l].y, tree->points[l].z );
    }
    fprintf( f, "\n" );

    fprintf( f, "LINES 1 %d\n", tree->edges[i].n_data + 1 );
    fprintf( f, "%d", tree->edges[i].n_data );
    for ( j=0; j<tree->edges[i].n_data; j++ )
      fprintf( f, " %d", j );
    fprintf( f, "\n" );


    /* write point's color
     */
    switch( par->color ) {
    default :
    case _NO_COLOR_ :
      break;

    case _COLOR_ANCHOR_ :
    case _COLOR_BRANCH_ORDER_ :
    case _COLOR_BRANCH_LABEL_ :
    case _COLOR_COMPONENT_LABEL_ :

      fprintf( f, "POINT_DATA %d\n", tree->edges[i].n_data );
      fprintf( f, "SCALARS volume float\n" );
      switch( par->colortable ) {
      default :
      case _DEFAULT_COLORTABLE_ :
        fprintf( f, "LOOKUP_TABLE default\n" );
        break;
      case _USER_COLORTABLE_ :
        if ( par->color == _COLOR_ANCHOR_ )
          fprintf( f, "LOOKUP_TABLE userdefined\n" );
        else
          fprintf( f, "LOOKUP_TABLE default\n" );
      }

      for ( j=0; j<tree->edges[i].n_data; j++ ) {
          pt = &(tree->points[ tree->edges[i].data[j] ]);
          /* for first and last points, if they are junction points
           * pick the color of the neighboring point
           * recall
           * - that junction points are replicated, thus multiple colors
           *   are allowed
           * - that colors are interpolated between points, so a single color
           *   at junction won't make sense
           */
          if ( j == 0 && pt->type == _TREE_JUNCTION_POINT_ ) {
            if ( tree->edges[i].n_data <= 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, can not specify color for point #%d of edge #%d\n", proc, j, i );
            }
            else
              pt = &(tree->points[ tree->edges[i].data[1] ]);
          }
          else if ( j == tree->edges[i].n_data-1 && pt->type == _TREE_JUNCTION_POINT_ ) {
            if ( tree->edges[i].n_data <= 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, can not specify color for point #%d of edge #%d\n", proc, j, i );
            }
            else
              pt = &(tree->points[ tree->edges[i].data[tree->edges[i].n_data-2] ]);
          }
          fprintf( f, "%d\n", pt->color );
      }
      fprintf( f, "\n" );

      switch( par->colortable ) {
      default :
      case _DEFAULT_COLORTABLE_ :
        break;
      case _USER_COLORTABLE_ :
        if ( par->color == _COLOR_ANCHOR_ ) {
          fprintf( f, "LOOKUP_TABLE userdefined 2\n" );
          fprintf( f, "%f %f %f 1.0\n", par->table[0][0], par->table[0][1], par->table[0][2] );
          fprintf( f, "%f %f %f 1.0\n", par->table[1][0], par->table[1][1], par->table[1][2] );
          fprintf( f, "\n" );
        }
      }

      break;
    }

    fclose( f );
  }

  return( 1 );
}





int _writeVTK_multipleEdges( char *name, typeTree *tree, typeTreeOutputParam *par, char *desc )
{
  char *proc = "_writeVTK_multipleEdges";
  int nendedgepts, ninneredgepts;
  int nendedges, ninneredges;
  int npts=0, nedges=0;
  int i, j, l;
  typeTreeFloatPoint *pt;
  FILE *f;



  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening '%s'\n", proc, name );
    return( -1 );
  }

  fprintf( f, "# vtk DataFile Version 2.0\n" );

  if ( desc != (char*)NULL )
      fprintf( f, "%s", desc );
  fprintf( f, "\n" );

  fprintf( f, "ASCII\n" );
  fprintf( f, "DATASET POLYDATA\n" );



  /* on compte les points et les aretes devant etre
   * ecrites
   * les jonctions sont dupliquees autant que necessaire
   * donc le nombre de points est celui des points des
   * aretes
   */
  for ( i=0, nendedgepts=0, ninneredgepts=0, nendedges=0, ninneredges=0; i<tree->n_edges; i++ ) {
    if ( tree->edges[i].type == _END_EDGE_ ) {
      nendedgepts += tree->edges[i].n_data;
      nendedges ++;
    }
    else if ( tree->edges[i].type == _INNER_EDGE_ ) {
      ninneredgepts += tree->edges[i].n_data;
      ninneredges ++;
    }
    else {
      fclose( f );
      if ( _verbose_ )
        fprintf( stderr, "%s: unknown type for edge #%d\n", proc, i );
      return( -1 );
    }
  }

  /* on ecrit les points et les aretes
   */
  nedges = ninneredges+nendedges;
  npts = ninneredgepts+nendedgepts;
  fprintf( f, "POINTS %d float\n", npts );
  for ( i=0; i<tree->n_edges; i++ ) {
    for ( j=0; j<tree->edges[i].n_data; j++ ) {
      l = tree->edges[i].data[j];
      fprintf( f, "%f %f %f\n", tree->points[l].x, tree->points[l].y, tree->points[l].z );
    }
  }
  fprintf( f, "\n" );
  fprintf( f, "LINES %d %d\n", nedges, nedges+npts );
  for ( i=0, l=0; i<tree->n_edges; i++ ) {
    fprintf( f, "%d", tree->edges[i].n_data );
    for ( j=0; j<tree->edges[i].n_data; j++, l++ )
      fprintf( f, " %d", l );
    fprintf( f, "\n" );
  }
  fprintf( f, "\n" );



  /* couleur
   */
  switch( par->color ) {
  default :
  case _NO_COLOR_ :
    break;

  case _COLOR_ANCHOR_ :
  case _COLOR_BRANCH_ORDER_ :
  case _COLOR_BRANCH_LABEL_ :
  case _COLOR_COMPONENT_LABEL_ :

    fprintf( f, "POINT_DATA %d\n", npts );
    fprintf( f, "SCALARS volume float\n" );
    switch( par->colortable ) {
    default :
    case _DEFAULT_COLORTABLE_ :
      fprintf( f, "LOOKUP_TABLE default\n" );
      break;
    case _USER_COLORTABLE_ :
      if ( par->color == _COLOR_ANCHOR_ )
        fprintf( f, "LOOKUP_TABLE userdefined\n" );
      else
        fprintf( f, "LOOKUP_TABLE default\n" );
    }

    for ( i=0; i<tree->n_edges; i++ ) {
      for ( j=0; j<tree->edges[i].n_data; j++ ) {
          pt = &(tree->points[ tree->edges[i].data[j] ]);
          /* for first and last points, if they are junction points
           * pick the color of the neighboring point
           * recall
           * - that junction points are replicated, thus multiple colors
           *   are allowed
           * - that colors are interpolated between points, so a single color
           *   at junction won't make sense
           */
          if ( j == 0 && pt->type == _TREE_JUNCTION_POINT_ ) {
            if ( tree->edges[i].n_data <= 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, can not specify color for point #%d of edge #%d\n", proc, j, i );
            }
            else
              pt = &(tree->points[ tree->edges[i].data[1] ]);
          }
          else if ( j == tree->edges[i].n_data-1 && pt->type == _TREE_JUNCTION_POINT_ ) {
            if ( tree->edges[i].n_data <= 1 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, can not specify color for point #%d of edge #%d\n", proc, j, i );
            }
            else
              pt = &(tree->points[ tree->edges[i].data[tree->edges[i].n_data-2] ]);
          }
          fprintf( f, "%d\n", pt->color );

      }
    }
    fprintf( f, "\n" );


    switch( par->colortable ) {
    default :
    case _DEFAULT_COLORTABLE_ :
      break;
    case _USER_COLORTABLE_ :
      if ( par->color == _COLOR_ANCHOR_ ) {
        fprintf( f, "LOOKUP_TABLE userdefined 2\n" );
        fprintf( f, "%f %f %f 1.0\n", par->table[0][0], par->table[0][1], par->table[0][2] );
        fprintf( f, "%f %f %f 1.0\n", par->table[1][0], par->table[1][1], par->table[1][2] );
        fprintf( f, "\n" );
      }
    }

    break;
  }

  fclose( f );
  return( 1 );
}





int writeVTKLegacyFile( char *name, typeTree *tree, typeTreeOutputParam *par, char *desc )
{
  char *proc = "writeVTKLegacyFile";

  if ( name == (char*)NULL || name[0] == '\0' || name[0] == '>' )
    return( 1 );

  /* at this point, the color (marked field of typeTreeFloatPoint)
   * is either 0 or a value given by a main axon image
   *
   *
   */
  _colorTree( tree, par->color );

  switch( par->filetype ) {

  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such type not handled yet\n", proc );
      fprintf(stderr, "\t '%s' will not be written\n", name );
    }
    return( -1 );
  case _MULTIPLE_FILES_ :
      return( _writeVTK_singleEdges( name, tree, par, desc ) );
  case _SINGLE_FILE_ :
    return( _writeVTK_multipleEdges( name, tree, par, desc ) );
  }

  return( 1 );
}












