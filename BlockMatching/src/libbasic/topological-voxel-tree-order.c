/****************************************************
 * topological-voxel-tree-order.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 15 jui 2017 09:29:15 CEST
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

#include <threshold.h>

#include <topological-voxel-tree-order.h>







/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalVoxelTreeOrder( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalVoxelTreeOrder(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalVoxelTreeOrder(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalVoxelTreeOrder( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalVoxelTreeOrder(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalVoxelTreeOrder(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}










/**************************************************
 *
 *
 *
 **************************************************/

static void _fprintfVoxelTreeOrder( FILE *f, typeVoxelTree *iTree )
{
  int i;
  voxelTreeComponent *c;

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch ( c->type ) {
    default :
      fprintf( f, "#%3d - unknown  -", i );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      fprintf( f, "#%3d - edge     -", i );
      if ( c->infoFirstPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " (%3d, %d)", c->indexFirstJunction, c->infoFirstPoint.order );
      fprintf( f, " <- ->" );
      fprintf( f, " (%3d, %d)", c->indexLastJunction, c->infoLastPoint.order );
      if ( c->infoLastPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_JUNCTION_ :
      fprintf( f, "#%3d - junction -", i );
      if ( c->infoJunction.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " %d", c->infoJunction.order );
      fprintf( f, "\n" );
      break;
    }
  }
}





/* initialization
 * - length is 0 for edges with both ends labeled
 * - length != 0 for edges for one edge labeled
 * thoses edges are no more in queue
 */

static int _initOrderInVoxelTree( typeComponentImage *treeImage, typeVoxelTree *iTree, int *lastLabel )
{
  char *proc = "_initOrderInVoxelTree";
  int i, k;
  float dx, dy, dz;
  voxelTreeComponent *c, *j;
  float maxlength;


  /* compute the sum of all length
   * it is a maximal value for length computation
   */
  for ( maxlength=0.0, i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    if ( c->type != _INTERN_TREE_INNER_EDGE_
         && c->type != _INTERN_TREE_END_EDGE_
         && c->type != _INTERN_TREE_EDGE_ )
      continue;
    maxlength += iTree->data[i].length;
  }

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    if ( c->type == _INTERN_TREE_INNER_EDGE_
         || c->type == _INTERN_TREE_END_EDGE_
         || c->type == _INTERN_TREE_EDGE_ ) {
      c->infoFirstPoint.length = c->infoLastPoint.length = maxlength;
    }
    else if ( c->type == _INTERN_TREE_JUNCTION_ ) {
      c->infoJunction.length = maxlength;
    }
  }

  (*lastLabel) = 1;

  /* process edges
   * labels attributed so far are assumed to be
   * - either 1: typically belongs to the main axon
   * - or 0 : not the main axon
   */

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    if ( c->type != _INTERN_TREE_INNER_EDGE_
         && c->type != _INTERN_TREE_END_EDGE_
         && c->type != _INTERN_TREE_EDGE_ )
      continue;
    if ( c->n_data <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: empty component #%d ?!\n", proc, i );
      return( -1 );
    }

    /* one-point edge
     */
    if ( c->n_data == 1 ) {
      if ( c->data[0].anchor == 1 ) {
        c->infoFirstPoint.order   = c->infoLastPoint.order = 1;
        c->infoFirstPoint.length  = c->infoLastPoint.length = 0.0;
        c->infoFirstPoint.inqueue = c->infoLastPoint.inqueue = 0;
      }
      continue;
    }


    if ( c->data[0].anchor ) {

      c->infoFirstPoint.label_branch = 1;
      c->infoFirstPoint.order = 1;
      c->infoFirstPoint.length = 0.0;
      c->infoFirstPoint.inqueue = 0;

      if ( c->data[c->n_data-1].anchor ) {
        /* both ends are labeled
         */
        c->infoLastPoint.label_branch = 1;
        c->infoLastPoint.order = 1;
        c->infoLastPoint.length = 0.0;
        c->infoLastPoint.inqueue = 0;
      }
      else  {
        /* new component
         */
        (*lastLabel) ++;
        c->infoLastPoint.label_branch = (*lastLabel);
        c->infoLastPoint.order = 2;
        c->infoLastPoint.length = 0.0;
        for ( k=c->n_data-1; k>1 && c->data[k].anchor == 0; k-- ) {
          c->data[k].label_branch = (*lastLabel);
          dx = (float)(c->data[k].x - c->data[k-1].x) * treeImage->voxelSize[0];
          dy = (float)(c->data[k].y - c->data[k-1].y) * treeImage->voxelSize[1];
          dz = (float)(c->data[k].z - c->data[k-1].z) * treeImage->voxelSize[2];
          c->infoLastPoint.length += sqrt( dx*dx + dy*dy + dz*dz );
        }
        c->infoLastPoint.inqueue = 0;
      }

    }
    else {

      if ( c->data[c->n_data-1].anchor ) {

        c->infoLastPoint.label_branch = 1;
        c->infoLastPoint.order = 1;
        c->infoLastPoint.length = 0.0;
        c->infoLastPoint.inqueue = 0;

        /* new component
         */
        (*lastLabel) ++;
        c->infoFirstPoint.label_branch = (*lastLabel);
        c->infoFirstPoint.order = 2;
        c->infoFirstPoint.length = 0.0;
        for ( k=0; k<c->n_data-1 && c->data[k].anchor == 0; k++ ) {
          c->data[k].label_branch = (*lastLabel);
          dx = (float)(c->data[k].x - c->data[k+1].x) * treeImage->voxelSize[0];
          dy = (float)(c->data[k].y - c->data[k+1].y) * treeImage->voxelSize[1];
          dz = (float)(c->data[k].z - c->data[k+1].z) * treeImage->voxelSize[2];
          c->infoFirstPoint.length += sqrt( dx*dx + dy*dy + dz*dz );
        }
        c->infoFirstPoint.inqueue = 0;
      }
      else  {
        /* do nothing
         * both ends are not labeled
         */
        ;
      }

    }

  }

  /* process junctions
   */
  for ( i=1; i<iTree->n_data; i++ ) {

    c = &(iTree->data[i]);
    if ( c->type != _INTERN_TREE_INNER_EDGE_
         && c->type != _INTERN_TREE_END_EDGE_
         && c->type != _INTERN_TREE_EDGE_ )
      continue;
    if ( c->n_data <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: empty component #%d ?!\n", proc, i );
      return( -1 );
    }

    if ( c->infoFirstPoint.inqueue != c->infoLastPoint.inqueue ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: weird, component #%d has been uncompletely processed\n", proc, i );
      }
      return( -1 );
    }

    if ( c->infoFirstPoint.inqueue == 1 && c->infoLastPoint.inqueue == 1 )
      continue;

    if ( c->indexFirstJunction != -1 ) {
      j = &(iTree->data[c->indexFirstJunction]);
      if ( c->infoFirstPoint.length > 0.0 ) {
        j->infoJunction.label_branch = c->infoFirstPoint.label_branch;
        j->infoJunction.order = c->infoFirstPoint.order;
        dx = ((float)(c->data[0].x) - j->center.x) * treeImage->voxelSize[0];
        dy = ((float)(c->data[0].y) - j->center.y) * treeImage->voxelSize[1];
        dz = ((float)(c->data[0].z) - j->center.z) * treeImage->voxelSize[2];
        j->infoJunction.length = c->infoFirstPoint.length + sqrt( dx*dx + dy*dy + dz*dz );
        j->infoJunction.inqueue = 0;
      }
      else {
        j->infoJunction.label_branch = 1;
        j->infoJunction.order = 1;
        j->infoJunction.length = 0.0;
        j->infoJunction.inqueue = 0;
      }
    }

    if ( c->indexLastJunction != -1 ) {
      j = &(iTree->data[c->indexLastJunction]);
      if ( c->infoLastPoint.length > 0.0 ) {
        j->infoJunction.label_branch = c->infoLastPoint.label_branch;
        j->infoJunction.order = c->infoLastPoint.order;
        dx = ((float)(c->data[c->n_data-1].x) - j->center.x) * treeImage->voxelSize[0];
        dy = ((float)(c->data[c->n_data-1].y) - j->center.y) * treeImage->voxelSize[1];
        dz = ((float)(c->data[c->n_data-1].z) - j->center.z) * treeImage->voxelSize[2];
        j->infoJunction.length = c->infoLastPoint.length + sqrt( dx*dx + dy*dy + dz*dz );
        j->infoJunction.inqueue = 0;
      }
      else {
        j->infoJunction.label_branch = 1;
        j->infoJunction.order = 1;
        j->infoJunction.length = 0.0;
        j->infoJunction.inqueue = 0;
      }
    }

  }

  return( 1 );
}










/**************************************************
 *
 *
 *
 **************************************************/

static void _fprintfVoxelTreeLength( FILE *f, typeVoxelTree *iTree )
{
  int i;
  voxelTreeComponent *c;

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch ( c->type ) {
    default :
      fprintf( f, "#%3d - unknown  -", i );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      fprintf( f, "#%3d - edge     -", i );
      if ( c->infoFirstPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " (%3d, %f)", c->indexFirstJunction, c->infoFirstPoint.length );
      fprintf( f, " <- %f ->", c->length );
      fprintf( f, " (%3d, %f)", c->indexLastJunction, c->infoLastPoint.length );
      if ( c->infoLastPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_JUNCTION_ :
      fprintf( f, "#%3d - junction -", i );
      if ( c->infoJunction.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " %f", c->infoJunction.length );
      fprintf( f, "\n" );
      break;
    }
  }
}





#define NEDGES 10

static int _getJunctionNeighbors( typeVoxelTree *iTree, int junction,
                                  int *edges )
{
  char *proc = "_getJunctionNeighbors";
  int i, nedges = 0;
  voxelTreeComponent *c;

  for ( i=0; i<NEDGES; i++ )
    edges[i] = 0;

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    if ( c->type != _INTERN_TREE_INNER_EDGE_
         && c->type != _INTERN_TREE_END_EDGE_
         && c->type != _INTERN_TREE_EDGE_ )
      continue;
    if ( c->indexFirstJunction == junction || c->indexLastJunction == junction ) {
      if ( nedges >= NEDGES ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: too many neighbors\n", proc );
        return( -1 );
      }
      edges[ nedges ] = i;
      nedges ++;
    }
  }

  return( nedges );
}





static int _propagateLength( typeComponentImage *treeImage, typeVoxelTree *iTree )
{
  char *proc = "_propagateLength";
  float dl;
  int i, j, changes;
  int iteration = 0;
  voxelTreeComponent *c, *fj, *lj, *e;
  int nedges, edges[NEDGES];
  float minlength;


  dl = ( treeImage->voxelSize[0] < treeImage->voxelSize[1] ) ? treeImage->voxelSize[0]  : treeImage->voxelSize[1];
  if ( dl > treeImage->voxelSize[2] ) dl = treeImage->voxelSize[2];
  dl /= 2.0;

  do {

    for ( changes=0, i=1; i<iTree->n_data; i++ ) {

      c = &(iTree->data[i]);

      if ( c->type == _INTERN_TREE_REMOVED_COMPONENT_ ||
           c->type == _INTERN_TREE_MERGED_COMPONENT_ )
          continue;

      /* edge case
       * - update length at both ends with junctions,
       * - update length at both ends
       */
      if ( c->type == _INTERN_TREE_INNER_EDGE_ || c->type == _INTERN_TREE_END_EDGE_ ) {

        if ( c->infoFirstPoint.inqueue ==  0 ||  c->infoLastPoint.inqueue == 0 ) {
          if ( c->infoFirstPoint.inqueue !=  0 || c->infoLastPoint.inqueue != 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, edge #%d seems not to be completely processed\n", proc, i );
            return( -1 );
          }
          continue;
        }

        fj = (voxelTreeComponent*)NULL;
        lj = (voxelTreeComponent*)NULL;
        if ( c->indexFirstJunction > 0 ) {
          fj = &(iTree->data[c->indexFirstJunction]);
          c->infoFirstPoint.length =  fj->infoJunction.length;
        }
        if ( c->indexLastJunction > 0 ) {
          lj = &(iTree->data[c->indexLastJunction]);
          c->infoLastPoint.length =  lj->infoJunction.length;
        }

        if ( c->infoLastPoint.length > c->infoFirstPoint.length + c->length + dl ) {
          c->infoLastPoint.length = c->infoFirstPoint.length + c->length;
          if ( c->indexLastJunction > 0 )
            lj->infoJunction.length = c->infoLastPoint.length;
          changes ++;
        }
        else if ( c->infoFirstPoint.length > c->infoLastPoint.length + c->length + dl ) {
          c->infoFirstPoint.length = c->infoLastPoint.length + c->length;
          if ( c->indexFirstJunction > 0 )
            fj->infoJunction.length = c->infoFirstPoint.length;
          changes ++;
        }
        else {
          ;
        }

      }
      else if ( c->type == _INTERN_TREE_JUNCTION_ ) {

        if ( c->infoJunction.inqueue == 0 )
          continue;

        nedges = _getJunctionNeighbors( iTree, i, edges );
        if ( nedges == -1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when looking for neighbors of junction #%d\n", proc, i );
          return( -1 );
        }
        if ( nedges <= 2 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: weird number of neighbors (%d) for junction #%d\n", proc, nedges, i );
          return( -1 );
        }
        minlength = c->infoJunction.length;
        for ( j=0; j<nedges; j++ ) {
          e = &(iTree->data[ edges[j] ]);
          if ( e->indexFirstJunction == i ) {
            if ( minlength > e->infoFirstPoint.length )
              minlength = e->infoFirstPoint.length;
          }
          else if ( e->indexLastJunction == i ) {
            if ( minlength > e->infoLastPoint.length )
              minlength = e->infoLastPoint.length;
          }
          else {
            if ( _verbose_ )
              fprintf( stderr, "%s: edge #%d is not connected to junction #%d?!\n", proc, edges[j], i );
            return( -1 );
          }
        }
        if ( c->infoJunction.length > minlength + dl ) {
          c->infoJunction.length = minlength;
          changes ++;
        }


      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unexpected component type\n", proc );
        return( -1 );
      }


    }

    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: iteration %2d, changes = %3d\n", proc, iteration, changes );
    iteration ++;

  } while ( changes > 0 );


  return( 1 );
}









/**************************************************
 *
 *
 *
 **************************************************/

static void _fprintfVoxelTreeLabel( FILE *f, typeVoxelTree *iTree )
{
  int i;
  voxelTreeComponent *c;

  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch ( c->type ) {
    default :
      fprintf( f, "#%3d - unknown  -", i );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      fprintf( f, "#%3d - edge     -", i );
      if ( c->infoFirstPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " (%3d, %d)", c->indexFirstJunction, c->infoFirstPoint.label_branch );
      fprintf( f, " <- ->" );
      fprintf( f, " (%3d, %d)", c->indexLastJunction, c->infoLastPoint.label_branch );
      if ( c->infoLastPoint.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, "\n" );
      break;
    case _INTERN_TREE_JUNCTION_ :
      fprintf( f, "#%3d - junction -", i );
      if ( c->infoJunction.label_branch == 1 )
        fprintf( f, " M" );
      else
        fprintf( f, "  " );
      fprintf( f, " %d", c->infoJunction.label_branch );
      fprintf( f, "\n" );
      break;
    }
  }
}





static int _propagateLabels( typeVoxelTree *iTree, int order, int *lastLabel )
{
  char *proc = "_propagateLabel";
  int i, j, changes;
  int iteration = 0;
  voxelTreeComponent *c, *fj, *lj, *e;
  int nedges, edges[NEDGES];

  /* les composantes avec inqueue == 0 ont deja ete traitees
   * le champ order peut etre a
   * - 0 : initialisation
   * - order-1 (junction)
   * - order (edges non termines) et apres
   * le champ label peut etre a
   * - 0
   * - >0
   */
  do {

    for ( changes=0, i=1; i<iTree->n_data; i++ ) {

      c = &(iTree->data[i]);

      if ( c->type == _INTERN_TREE_REMOVED_COMPONENT_ ||
           c->type == _INTERN_TREE_MERGED_COMPONENT_ )
          continue;

      /* edge case
       */
      if ( c->type == _INTERN_TREE_INNER_EDGE_ || c->type == _INTERN_TREE_END_EDGE_ ) {

        if ( c->infoFirstPoint.inqueue ==  0 || c->infoLastPoint.inqueue == 0 ) {
          if ( c->infoFirstPoint.inqueue !=  0 || c->infoLastPoint.inqueue != 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, edge #%d seems not to be completely processed\n", proc, i );
            return( -1 );
          }
          continue;
        }

        if ( c->infoFirstPoint.label_branch > 0 && c->infoLastPoint.label_branch > 0 )
          continue;

        if ( (c->infoFirstPoint.label_branch == 0 && c->infoLastPoint.label_branch > 0)
             || (c->infoFirstPoint.label_branch > 0 && c->infoLastPoint.label_branch == 0) ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: weird, edge #%d seems not to be completely processed\n", proc, i );
          return( -1 );
        }

        if ( (c->infoFirstPoint.order == 0 && c->infoLastPoint.order > 0)
             || (c->infoFirstPoint.order > 0 && c->infoLastPoint.order == 0)
             || (c->infoFirstPoint.order > 0 && c->infoLastPoint.order > 0) ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: weird, edge #%d seems not to be completely processed\n", proc, i );
          return( -1 );
        }

        fj = (voxelTreeComponent*)NULL;
        lj = (voxelTreeComponent*)NULL;

        if ( c->indexFirstJunction > 0 ) {
          fj = &(iTree->data[c->indexFirstJunction]);
          if ( fj->infoJunction.order == order-1 ) {
            (*lastLabel) ++;
            c->infoFirstPoint.label_branch = c->infoLastPoint.label_branch = (*lastLabel);
            c->infoFirstPoint.order = c->infoLastPoint.order = order;
            changes ++;
          }
          else if ( fj->infoJunction.order == order ) {
            c->infoFirstPoint.label_branch = c->infoLastPoint.label_branch = fj->infoJunction.label_branch;
            c->infoFirstPoint.order = c->infoLastPoint.order = order;
            changes ++;
          }
          else if ( fj->infoJunction.order > 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, edge #%d is connected to junction #%d of order %d\n",
                       proc, i, c->indexFirstJunction, fj->infoJunction.order );
            return( -1 );
          }
        }

        if ( c->indexLastJunction > 0 ) {
          lj = &(iTree->data[c->indexLastJunction]);
          if ( lj->infoJunction.order > 0 && c->infoFirstPoint.order > 0 && c->infoLastPoint.order > 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, edge #%d seems to get order from junction #%d\n", proc, i, c->indexFirstJunction );
            return( -1 );
          }
          if ( lj->infoJunction.order == order-1 ) {
            (*lastLabel) ++;
            c->infoFirstPoint.label_branch = c->infoLastPoint.label_branch = (*lastLabel);
            c->infoFirstPoint.order = c->infoLastPoint.order = order;
            changes ++;
          }
          else if ( lj->infoJunction.order == order ) {
            c->infoFirstPoint.label_branch = c->infoLastPoint.label_branch = lj->infoJunction.label_branch;
            c->infoFirstPoint.order = c->infoLastPoint.order = order;
            changes ++;
          }
          else if ( lj->infoJunction.order > 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: weird, edge #%d is connected to junction #%d of order %d\n",
                       proc, i, c->indexLastJunction, lj->infoJunction.order );
            return( -1 );
          }
        }

      }
      else if ( c->type == _INTERN_TREE_JUNCTION_ ) {

        if ( c->infoJunction.inqueue == 0 )
          continue;

        if ( c->infoJunction.label_branch > 0 )
          continue;

        if ( c->infoJunction.order > 0 )
          continue;

        nedges = _getJunctionNeighbors( iTree, i, edges );
        if ( nedges == -1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: error when looking for neighbors of junction #%d\n", proc, i );
          return( -1 );
        }
        if ( nedges <= 2 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: weird number of neighbors (%d) for junction #%d\n", proc, nedges, i );
          return( -1 );
        }

        for ( j=0; j<nedges; j++ ) {
          e = &(iTree->data[ edges[j] ]);

          if ( e->indexFirstJunction == i ) {
            if ( e->infoFirstPoint.order == order ) {
              c->infoJunction.order = e->infoFirstPoint.order;
              c->infoJunction.label_branch = e->infoFirstPoint.label_branch;
              changes++;
            }
            else if ( e->infoFirstPoint.order > 0 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, junction #%d is connected to edge #%d of order %d\n",
                         proc, i, edges[j], e->infoFirstPoint.order );
              return( -1 );
            }
          }
          else if ( e->indexLastJunction == i ) {
            if ( e->infoLastPoint.order == order ) {
              c->infoJunction.order = e->infoLastPoint.order;
              c->infoJunction.label_branch = e->infoLastPoint.label_branch;
              changes++;
            }
            else if ( e->infoLastPoint.order > 0 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: weird, junction #%d is connected to edge #%d of order %d\n",
                         proc, i, edges[j], e->infoLastPoint.order );
              return( -1 );
            }
          }
          else {
            if ( _verbose_ )
              fprintf( stderr, "%s: edge #%d is not connected to junction #%d?!\n", proc, edges[j], i );
            return( -1 );
          }
        }

      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unexpected component type\n", proc );
        return( -1 );
      }

    }

    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: iteration %2d, changes = %3d\n", proc, iteration, changes );
    iteration ++;

  } while ( changes > 0 );


  return( 1 );
}





static int _extractMaxBranches( typeVoxelTree *iTree, int order )
{
  char *proc= "_extractMaxBranches";
  int i, j, label;
  voxelTreeComponent *c, *e;
  int imax, next, current;
  float lengthmax;
  int nedges, edges[NEDGES];

  for ( i=1; i<iTree->n_data; i++ ) {

    c = &(iTree->data[i]);

    if ( c->type != _INTERN_TREE_END_EDGE_ ) continue;

    if ( c->infoFirstPoint.inqueue ==  0 || c->infoLastPoint.inqueue == 0 ) {
      if ( c->infoFirstPoint.inqueue !=  0 || c->infoLastPoint.inqueue != 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: weird, edge #%d seems not to be completely processed\n", proc, i );
        return( -1 );
      }
      continue;
    }

    if ( c->infoFirstPoint.label_branch == 0 || c->infoLastPoint.label_branch == 0 ) continue;
    if ( c->infoFirstPoint.order == 0 || c->infoLastPoint.order == 0 ) continue;

    /* found an end branch for label
     */
    if ( c->infoFirstPoint.order != order || c->infoLastPoint.order != order ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird, edge #%d has orders (%d,%d) instead of %d\n", proc, i,
                 c->infoFirstPoint.order, c->infoLastPoint.order, order );
      return( -1 );
    }
    label = c->infoFirstPoint.label_branch;
    imax = i;
    if ( c->infoFirstPoint.length > c->infoLastPoint.length ) {
      lengthmax = c->infoFirstPoint.length;
    }
    else {
      lengthmax = c->infoLastPoint.length;
    }

    /* look for a better maximum
     */
    for ( j=1; j<iTree->n_data; j++ ) {
      c = &(iTree->data[j]);
      if ( c->type != _INTERN_TREE_END_EDGE_ ) continue;
      if ( c->infoFirstPoint.label_branch != label && c->infoLastPoint.label_branch != label )
        continue;
      if ( (c->infoFirstPoint.label_branch != label && c->infoLastPoint.label_branch == label)
           || (c->infoFirstPoint.label_branch == label && c->infoLastPoint.label_branch != label) ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: weird, edge #%d has labels (%d,%d) instead of %d\n", proc, i,
                   c->infoFirstPoint.label_branch, c->infoLastPoint.label_branch, label );
        return( -1 );
      }
      if ( c->infoFirstPoint.length > lengthmax ) {
        lengthmax = c->infoFirstPoint.length;
        imax = j;
      }
      if ( c->infoLastPoint.length > lengthmax ) {
        lengthmax = c->infoLastPoint.length;
        imax = j;
      }
    }



    /* maximum length of subtree label is reached for component imax
     * back track
     */

    current = imax;
    do {
      c = &(iTree->data[current]);
      if ( c->type == _INTERN_TREE_INNER_EDGE_ || c->type == _INTERN_TREE_END_EDGE_ ) {

        c->infoFirstPoint.inqueue = c->infoLastPoint.inqueue = 0;
        if ( c->indexFirstJunction > 0 ) {
          if ( c->indexLastJunction > 0 ) {
            if ( iTree->data[c->indexFirstJunction].infoJunction.length < iTree->data[c->indexLastJunction].infoJunction.length )
              next = c->indexFirstJunction;
            else
              next = c->indexLastJunction;
          }
          else {
            next = c->indexFirstJunction;
          }
        }
        else if ( c->indexLastJunction > 0 ) {
          next = c->indexLastJunction;
        }
        /* end condition
         */
        if ( iTree->data[ next ].infoJunction.inqueue == 0 )
          break;

      }
      else if ( c->type == _INTERN_TREE_JUNCTION_ ) {

        c->infoJunction.inqueue = 0;
        nedges = _getJunctionNeighbors( iTree, current, edges );
        for ( j=0; j<nedges; j++ ) {
          e = &(iTree->data[ edges[j] ]);
          if ( e->indexFirstJunction == current ) {
            if ( e->infoLastPoint.length < c->infoJunction.length )
              next = edges[j];
          }
          else if ( e->indexLastJunction == current ) {
            if ( e->infoFirstPoint.length < c->infoJunction.length )
              next = edges[j];
          }
          else {
            if ( _verbose_ )
              fprintf( stderr, "%s: edge #%d is not connected to junction #%d?!\n", proc, edges[j], current );
            return( -1 );
          }
        }
        /* end condition
         */
        if ( iTree->data[ next ].infoFirstPoint.inqueue == 0 || iTree->data[ next ].infoLastPoint.inqueue == 0 )
          break;

      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unexpected component type\n", proc );
        return( -1 );
      }
      current = next;
    } while ( 1 );



    /* set label and order to 0 for other parts of the subtree
     */
    for ( j=1; j<iTree->n_data; j++ ) {
      c = &(iTree->data[j]);

      if ( c->type == _INTERN_TREE_REMOVED_COMPONENT_ ||
           c->type == _INTERN_TREE_MERGED_COMPONENT_ )
          continue;

      if ( c->type == _INTERN_TREE_INNER_EDGE_ || c->type == _INTERN_TREE_END_EDGE_ ) {
        if ( c->infoFirstPoint.label_branch != label && c->infoLastPoint.label_branch != label )
          continue;
        if ( c->infoFirstPoint.inqueue == 0 && c->infoLastPoint.inqueue == 0 )
          continue;
        c->infoFirstPoint.order = c->infoLastPoint.order = 0;
        c->infoFirstPoint.label_branch = c->infoLastPoint.label_branch = 0;
      }
      else if ( c->type == _INTERN_TREE_JUNCTION_ ) {
        if ( c->infoJunction.label_branch != label )
          continue;
        if ( c->infoJunction.inqueue == 0 )
          continue;
        c->infoJunction.order = 0;
        c->infoJunction.label_branch = 0;
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: unexpected component type\n", proc );
        return( -1 );
      }
    }

  }

  return( 1 );
}






/**************************************************
 *
 *
 *
 **************************************************/

int branchesVoxelTree( typeComponentImage *treeImage, typeVoxelTree *iTree )
{
  char *proc = "branchesVoxelTree";
  int lastLabel = 1;
  int previousLastLabel;
  int order = 2;
  int i, n;
  voxelTreeComponent *c;




  if ( _initOrderInVoxelTree( treeImage, iTree, &lastLabel ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when initializing propagation\n", proc );
    return( -1 );
  }



  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: after initialization, last label = %d\n", proc, lastLabel );
  }

  if ( _propagateLength( treeImage, iTree ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when propagate distance\n", proc );
    return( -1 );
  }

  if ( 0 && _verbose_ >= 2 ) {
    fprintf( stderr, "LENGTH \n" );
    _fprintfVoxelTreeLength( stderr, iTree );
    fprintf( stderr, "\n" );
  }

  if ( 0 && _verbose_ >= 2 ) {
    fprintf( stderr, "ORDER \n" );
    _fprintfVoxelTreeOrder( stderr, iTree );
    fprintf( stderr, "\n" );
  }

  if ( 0 && _verbose_ >= 2 ) {
    fprintf( stderr, "LABEL \n" );
    _fprintfVoxelTreeLabel( stderr, iTree );
    fprintf( stderr, "\n" );
  }


  do {
    /* labeling subtrees
     */
    previousLastLabel = lastLabel;
    if ( _propagateLabels( iTree, order, &lastLabel ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when propagating labels for order %d\n", proc, order );
      return( -1 );
    }
    /* end condition
     * no more new branches
     */
    if ( previousLastLabel == lastLabel )
      break;

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: after propagating labels for order %d, last label = %d\n", proc, order, lastLabel );
    }

    if ( _extractMaxBranches( iTree, order ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when extracting maximal branches for order %d\n", proc, order );
      return( -1 );
    }


    order ++;
  } while( 1 );


  /* propagate branch label and order to points
   */
  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch( c->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such component type not handled yet\n", proc );
      return( -1 );
    case _INTERN_TREE_REMOVED_COMPONENT_ :
    case _INTERN_TREE_MERGED_COMPONENT_ :
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      for ( n=0; n<c->n_data; n++ ) {
        if ( c->infoFirstPoint.label_branch == c->infoLastPoint.label_branch ) {
          c->data[n].label_branch = c->infoFirstPoint.label_branch;
          c->data[n].order_branch = c->infoFirstPoint.order;
        }
        else {
          if ( c->data[n].anchor ) {
            if ( c->data[0].anchor ) {
              c->data[n].label_branch = c->infoFirstPoint.label_branch;
              c->data[n].order_branch = c->infoFirstPoint.order;
            }
            else if ( c->data[c->n_data-1].anchor ) {
              c->data[n].label_branch = c->infoLastPoint.label_branch;
              c->data[n].order_branch = c->infoLastPoint.order;
            }
            else {
              if ( _verbose_ )
                fprintf( stderr, "%s: anchor point in an edge ended by non-anchor points ?!\n", proc );
              return( -1 );
            }
          }
          else {
            if ( c->data[0].anchor == 0 ) {
              c->data[n].label_branch = c->infoFirstPoint.label_branch;
              c->data[n].order_branch = c->infoFirstPoint.order;
            }
            else if ( c->data[c->n_data-1].anchor == 0 ) {
              c->data[n].label_branch = c->infoLastPoint.label_branch;
              c->data[n].order_branch = c->infoLastPoint.order;
            }
            else {
              if ( _verbose_ )
                fprintf( stderr, "%s: non-anchor point in an edge ended by anchor points ?!\n", proc );
              return( -1 );
            }
          }
        }
      }
      break;
    case _INTERN_TREE_JUNCTION_ :
      for ( n=0; n<c->n_data; n++ ) {
          c->data[n].label_branch = c->infoJunction.label_branch;
          c->data[n].order_branch = c->infoJunction.order;
      }
      break;
    }
  }




  return( 1 );
}






/**************************************************
 *
 *
 *
 **************************************************/

int maxOrderVoxelTree( typeVoxelTree *iTree )
{
  char *proc = "maxOrderVoxelTree";
  int maxorder = 0;
  int i;
  voxelTreeComponent *c;


  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    switch( c->type ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: such component type not handled yet\n", proc );
      return( -1 );
    case _INTERN_TREE_REMOVED_COMPONENT_ :
    case _INTERN_TREE_MERGED_COMPONENT_ :
      break;
    case _INTERN_TREE_INNER_EDGE_ :
    case _INTERN_TREE_END_EDGE_ :
    case _INTERN_TREE_EDGE_ :
      if ( maxorder < c->infoFirstPoint.order )
          maxorder = c->infoFirstPoint.order;
      if ( maxorder < c->infoLastPoint.order )
          maxorder = c->infoLastPoint.order;
      break;
    case _INTERN_TREE_JUNCTION_ :
        if ( maxorder < c->infoJunction.order )
            maxorder = c->infoJunction.order;
      break;
    }
  }

  return( maxorder );
}





/**************************************************
 *
 *
 *
 **************************************************/

int branchesTreeImage( void* theBuf, bufferType theType,
                       void* ancBuf, bufferType ancType,
                       void* resBuf, bufferType resType,
                       int *theDim, float *theSize )
{
  char *proc = "branchesTreeImage";
  typeComponentImage treeImage;
  typeVoxelTree iTree;

  /* build tree image
   */

  initComponentImage( &treeImage );
  if ( allocComponentImage( &treeImage, USHORT, theDim ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error for tree image\n", proc );
    return( -1 );
  }
  treeImage.voxelSize[0] = theSize[0];
  treeImage.voxelSize[1] = theSize[1];
  treeImage.voxelSize[2] = theSize[2];

  if ( imageToComponentImage( theBuf, theType, &treeImage ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when labeling components in tree\n", proc );
    return( -1 );
  }

  /* build internal tree
   * all edges can be pruned after voxel tree building
   */

  initVoxelTree( &iTree );
  if ( componentImageToVoxelTree( &treeImage, ancBuf, ancType, &iTree ) != 1 ) {
    freeVoxelTree( &iTree );
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building internal tree\n", proc );
    return( -1 );
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "\n" );
    fprintfSummaryVoxelTree( stderr, &iTree, "before pruning" );
  }

  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\n" );
    fprintVoxelTree( stderr, &iTree, 0 );
  }


  /* compute branches by order
   */


  if ( branchesVoxelTree(  &treeImage, &iTree ) != 1 ) {
    freeVoxelTree( &iTree );
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when processing internal tree\n", proc );
    return( -1 );
  }


  freeVoxelTree( &iTree );

  if ( thresholdBuffer( treeImage.componentBuf, treeImage.componentBufferType,
                      resBuf, resType, treeImage.theDim, 1.0 ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when converting component image\n", proc );
    return( -1 );
  }

  freeComponentImage( &treeImage );

  return( 1 );
}


















