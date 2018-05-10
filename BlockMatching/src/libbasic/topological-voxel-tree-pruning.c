/****************************************************
 * topological-voxel-tree-pruning.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Lun 12 jui 2017 11:30:52 CEST
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

#include <connexe.h>
#include <convert.h>
#include <t06t26.h>
#include <threshold.h>
#include <topological-operations-common.h>
#include <vtmalloc.h>


#include <topological-component-image.h>
#include <topological-voxel-tree.h>
#include <topological-voxel-tree-pruning.h>







/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalVoxelTreePruning( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalVoxelTreePruning(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalVoxelTreePruning(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalVoxelTreePruning( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalVoxelTreePruning(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalVoxelTreePruning(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}











/**************************************************
 *
 *
 *
 **************************************************/

typedef enum enumSmallestEdgeTest {
  _CHECK_REAL_LENGTH_FIRST_,
  _CHECK_VOXEL_LENGTH_FIRST_
} enumSmallestEdgeTest;



static int _getSmallestEndEdge( typeVoxelTree *iTree,
                             enumSmallestEdgeTest test )
{
  char *proc = "_getSmallestEndEdge";
  int i;
  int label = -1;

  if ( _debug_ )
    fprintf( stderr, "%s:\n", proc );

  for ( i=0; i<iTree->n_data; i++ ) {
    if ( iTree->data[i].type != _INTERN_TREE_END_EDGE_ ) {
        continue;
    }
    if ( iTree->data[i].canbepruned == 0 ) {
        if ( _debug_ )
          fprintf( stderr, "    skip component #%d: can not be pruned\n", i );
        continue;
    }

    if ( label == -1 ) {
      label = i;
    }
    else {
        switch( test ) {
        default :
        case _CHECK_REAL_LENGTH_FIRST_ :
            if ( iTree->data[label].length > iTree->data[i].length )
              label = i;
            break;
        case _CHECK_VOXEL_LENGTH_FIRST_ :
            if ( iTree->data[label].n_data > iTree->data[i].n_data )
              label = i;
            break;
        }
    }

  }

  return( label );
}





/* delete the end edge
 * returns:
 * -1 : error
 *  0 : uncomplete deletion
 *  1 : complete deletion
 */
static int _eraseEndEdge( typeComponentImage *treeImage,
                          typeVoxelTree *iTree, int label )
{
  char *proc = "_eraseEndEdge";
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int p, ndeletion = 0;
  int stoppruning;
  int i, j, k;
  size_t indp;


  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *resBuf = (u16*)(treeImage->componentBuf);

      if ( _verbose_ >= 3 ) {
        fprintf( stderr, "\t - remove component %d in [%d %d %d]x[%d %d %d]\n",
                 label, iTree->data[label].minCorner[0], iTree->data[label].minCorner[1],
                 iTree->data[label].minCorner[2], iTree->data[label].maxCorner[0],
                 iTree->data[label].maxCorner[1], iTree->data[label].maxCorner[2]);
      }

      /* the component is removed from the end point
       *
       */
      if ( iTree->data[label].indexFirstJunction <= 0 ) {
        for ( p=0, stoppruning=0, ndeletion=0; p<iTree->data[label].n_data && stoppruning==0; p++ ) {
          if ( iTree->data[label].data[p].anchor ) {
            stoppruning = 1;
            if ( _verbose_ >= 3 ) {
              fprintf( stderr, "\t - stop pruning component #%d at mark (%d,%d,%d)", label,
                       iTree->data[label].data[p].x,
                       iTree->data[label].data[p].y,
                       iTree->data[label].data[p].z );
              fprintf( stderr, ", remove %d points out of %d", p, iTree->data[label].n_data-1 );
              fprintf( stderr, "\n" );
            }
            continue;
          }
          i = iTree->data[label].data[p].x;
          j = iTree->data[label].data[p].y;
          k = iTree->data[label].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          resBuf[indp] = 0;
          ndeletion ++;
        }
        if ( stoppruning == 1 ) {
          for ( p=0; p<iTree->data[label].n_data-ndeletion; p++ ) {
            iTree->data[label].data[p] = iTree->data[label].data[p+ndeletion];
          }
          iTree->data[label].n_data -= ndeletion;
        }
      }
      else if ( iTree->data[label].indexLastJunction <= 0 ) {
        for ( p=iTree->data[label].n_data-1, stoppruning=0, ndeletion=0; p>=0 && stoppruning==0; p-- ) {
          if ( iTree->data[label].data[p].anchor ) {
            stoppruning = 1;
            if ( _verbose_ >= 3 ) {
              fprintf( stderr, "\t - stop pruning component #%d at mark (%d,%d,%d)", label,
                       iTree->data[label].data[p].x,
                       iTree->data[label].data[p].y,
                       iTree->data[label].data[p].z );
              fprintf( stderr, ", remove %d points out of %d", iTree->data[label].n_data-1-p, iTree->data[label].n_data-1 );
              fprintf( stderr, "\n" );
            }
            continue;
          }
          i = iTree->data[label].data[p].x;
          j = iTree->data[label].data[p].y;
          k = iTree->data[label].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          resBuf[indp] = 0;
          ndeletion ++;
        }
        if ( stoppruning == 1 ) {
          iTree->data[label].n_data -= ndeletion;
        }
      }
      else {
        if ( _verbose_ )
          fprintf( stderr, "%s: weird neighbor configuration for component to be removed\n", proc );
        if ( _verbose_ >= 2 ) {
          fprintf( stderr, "\t   component %d is connected to %d and %d\n",
                   label, iTree->data[label].indexFirstJunction,
                   iTree->data[label].indexLastJunction );
        }
        return( -1 );
      }
    }
    break;
  }


  /* the component has not been completely deleted
   * update its length
   */
  if ( stoppruning == 1 ) {
    iTree->data[label].canbepruned = 0;
    lengthEdgeVoxelTree( iTree, label );
    return( 0 );
  }

  /* the component has been completely deleted
   */
  if ( iTree->data[label].n_data == 0 )
    iTree->data[label].type = _INTERN_TREE_REMOVED_COMPONENT_;
  return( 1 );
}





/* delete the simple points from the junction
 * returns:
 * -1 : error
 *  1 : succees
 */
static int _eraseJunction( typeComponentImage *treeImage,
                           typeVoxelTree *iTree, int junction )
{
  char *proc = "_eraseJunction";
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  int totalremoved = 0;
  int p, remove;
  int i, j, k;
  int x, y, z;
  int n, neighbors[27], anchorNeighbors[27];
  size_t indn, indp;
  int t06, t26;
  voxelTreeIntPoint pt;
  voxelTreeComponent *c;
  int ptc;

  if ( _debug_ ) {
    fprintf( stderr, "%s: erase junction #%d\n", proc, junction );
  }

  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *resBuf = (u16*)(treeImage->componentBuf);

      /* remove junction points
       */
      do {
        remove = 0;
        for ( p=0; p<iTree->data[junction].n_data; p++ ) {
          i = iTree->data[junction].data[p].x;
          j = iTree->data[junction].data[p].y;
          k = iTree->data[junction].data[p].z;



          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* this point has already been deleted
           */
          if ( resBuf[indp] == 0 ) continue;

          for ( n=0, z=k-1; z<=k+1; z++ )
          for ( y=j-1; y<=j+1; y++ )
          for ( x=i-1; x<=i+1; x++, n++ ) {
           if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
             neighbors[n] = 0;
             anchorNeighbors[n] = 0;
             continue;
           }
           indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
           neighbors[n] = resBuf[indn];
           if ( neighbors[n] == 0 ) {
             anchorNeighbors[n] = 0;
           }
           else {
             c = &(iTree->data[ resBuf[indn] ]);
             for ( ptc=0; ptc<c->n_data; ptc++ ) {
               if ( c->data[ptc].x != x || c->data[ptc].y != y || c->data[ptc].z != z )
                 continue;
               anchorNeighbors[n] = ( c->data[ptc].anchor ) ? resBuf[indn] : 0;
             }
           }
          }

          if ( 0 && _debug_ ) {
              fprintf( stderr, "     neighborhood of (%d %d %d]\n", i, j, k );
              fprintfNeighborhood( stderr, neighbors, 27 );
              fprintf( stderr, "     anchor neighborhood of (%d %d %d]\n", i, j, k );
              fprintfNeighborhood( stderr, anchorNeighbors, 27 );
          }

          Compute_T06_and_T26( neighbors, &t06, &t26 );
          if ( t06 != 1 || t26 != 1 ) continue;

          if ( iTree->data[junction].data[p].anchor ) {
            Compute_T06_and_T26( anchorNeighbors, &t06, &t26 );
            if ( t06 != 1 || t26 != 1 ) continue;
            if ( _verbose_ >= 3 )
              fprintf( stderr, "\t   remove *anchor* point [%d,%d,%d]\n", i, j, k );
          }
          else {
            if ( _verbose_ >= 3 )
              fprintf( stderr, "\t   remove point [%d,%d,%d]\n", i, j, k );
          }
          resBuf[indp] = 0;
          remove ++;

          /* if not the last point
           * switch with the last point, update the number of points
           */
          if ( p < iTree->data[junction].n_data-1 ) {
            pt = iTree->data[junction].data[ iTree->data[junction].n_data-1 ];
            iTree->data[junction].data[ iTree->data[junction].n_data-1 ] = iTree->data[junction].data[ p ];
            iTree->data[junction].data[ p ] = pt;
            p--;
          }
          iTree->data[junction].n_data --;
        }
        totalremoved += remove;
      } while ( remove > 0 );
    }
    break;
  }



  return( totalremoved );
}





static void _updateIncidentEdges( typeVoxelTree *iTree,
                                 int junction )
{
  int i;

  for ( i=1; i<iTree->n_data; i++ ) {
    if ( iTree->data[i].type != _INTERN_TREE_INNER_EDGE_
         && iTree->data[i].type != _INTERN_TREE_END_EDGE_
         && iTree->data[i].type != _INTERN_TREE_EDGE_ )
      continue;
    if ( iTree->data[i].indexFirstJunction != junction
         && iTree->data[i].indexLastJunction != junction )
      continue;
    if ( _verbose_ >= 3 )
      fprintf( stderr, "\t   update length of edge %d \n", i );
    lengthEdgeVoxelTree( iTree, i );
  }
}










/**************************************************
 *
 *
 *
 **************************************************/

typedef enum _enumJunctionFuture {
  _FUTURE_UNKNOWN_,
  _FUTURE_IS_JUNCTION_,
  _FUTURE_IS_EDGE_,
  _FUTURE_IS_MIXED_
} _enumJunctionFuture;



/* count the neighbors for each point of a junction
 */

static _enumJunctionFuture _getJunctionFuture( typeComponentImage *treeImage,
                                            typeVoxelTree *iTree,
                                            int junction )
{
  char *proc = "_getJunctionFuture";
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  int njunctionpoints = 0;
  int nedgespoints = 0;
  int p;
  int i, j, k;
  int x, y, z;
  int n;
  size_t indn, indp;

  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *resBuf = (u16*)(treeImage->componentBuf);

      for ( p=0; p<iTree->data[junction].n_data; p++ ) {
        i = iTree->data[junction].data[p].x;
        j = iTree->data[junction].data[p].y;
        k = iTree->data[junction].data[p].z;
        indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
        if ( resBuf[indp] == 0 ) continue;

        for ( n=0, z=k-1; z<=k+1; z++ )
        for ( y=j-1; y<=j+1; y++ )
        for ( x=i-1; x<=i+1; x++ ) {
          if ( x == i && y == j && z == k )
            continue;
          if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
            continue;
          }
          indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
          if ( resBuf[indn] > 0 ) n++;
        }

        if ( n <= 2 ) nedgespoints ++;
        else njunctionpoints ++;
      }
    }
    break;
  }

  if ( njunctionpoints > 0 && nedgespoints == 0 ) {
    return( _FUTURE_IS_JUNCTION_ );
  }

  if ( njunctionpoints == 0 && nedgespoints > 0 ) {
    return( _FUTURE_IS_EDGE_ );
  }

  if ( njunctionpoints > 0 && nedgespoints > 0 ) {
    return( _FUTURE_IS_MIXED_ );
  }

  return( _FUTURE_UNKNOWN_ );
}











/**************************************************
 *
 * merge 2 edges and the connecting junction
 *
 **************************************************/

#define NEDGES 10





/* counts the neighboring components of the junction
 * returns:
 * -1 : error
 *  2 : no more a junction => merge
 * >2 : still a junction
 */
static int _getJunctionNeighbors( typeComponentImage *treeImage,
                               typeVoxelTree *iTree, int junction,
                               int *edges )
{
  char *proc = "_getJunctionNeighbors";
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  int p;
  int i, j, k;
  int x, y, z;
  int n, neighbors[27];
  size_t indn, indp;
  int t06, t26;
  int isajunction;
  int nedges;


  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *resBuf = (u16*)(treeImage->componentBuf);

      /* is the junction still a junction ?
       */
      isajunction = -1;
      for ( p=0; p<iTree->data[junction].n_data; p++ ) {
        i = iTree->data[junction].data[p].x;
        j = iTree->data[junction].data[p].y;
        k = iTree->data[junction].data[p].z;
        indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
        if ( resBuf[indp] == 0 ) continue;

        for ( n=0, z=k-1; z<=k+1; z++ )
        for ( y=j-1; y<=j+1; y++ )
        for ( x=i-1; x<=i+1; x++, n++ ) {
         if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
           neighbors[n] = 0;
           continue;
         }
         indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
         neighbors[n] = resBuf[indn];
        }

        Compute_T06_and_T26( neighbors, &t06, &t26 );
        if ( t26 == 1 ) {
          /* t26 == 1 may occur for "large" junctions
           */
          if ( 0 && _verbose_ ) {
            fprintf( stderr, "%s: weird neighborhood for junction #%d ( t26 == 1 )\n", proc, junction );
            fprintf(stderr, "     at point (%d,%d,%d)\n", i, j, k );
            fprintfNeighborhood( stderr, neighbors, 27 );
          }
        } else if ( t26 == 2 ) {
          if ( isajunction == -1 ) isajunction = 0;
        } else if ( t26 > 2 ) {
          isajunction = 1;
        } else {
          if ( _verbose_ ) {
            fprintf( stderr, "%s: weird neighborhood for junction #%d ( t26 < 1 )\n", proc, junction );
            fprintfNeighborhood( stderr, neighbors, 27 );
          }
        }
      }

      /* if the junction is not a junction,
       * it should be a joint between two branchs
       */
      if ( isajunction == -1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: weird, junction #%d not recognized as junction or edge\n", proc, junction );
        }
        return( -1 );
      }
      else if ( isajunction == 1 ) {
        /* it is still a junction,
         * there are at least 3 neighbors
         */
        return( 3 );
      }
      else if ( isajunction == 0 ) {
        /* identify neighboring components
         * be careful: all junction points may have two neighbors
         * but it still may be a junction
         * check it with the number of neighboring components
         */
        for ( n=0; n<NEDGES; n++ ) edges[n] = 0;

        nedges = 0;
        for ( p=0; p<iTree->data[junction].n_data; p++ ) {
          i = iTree->data[junction].data[p].x;
          j = iTree->data[junction].data[p].y;
          k = iTree->data[junction].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          if ( resBuf[indp] == 0 ) continue;

          for ( n=0, z=k-1; z<=k+1; z++ )
          for ( y=j-1; y<=j+1; y++ )
          for ( x=i-1; x<=i+1; x++, n++ ) {
           if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
             neighbors[n] = 0;
             continue;
           }
           indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
           neighbors[n] = resBuf[indn];
          }

          for ( n=0; n<27; n++ ) {
            if ( neighbors[n] == 0 ) continue;
            if ( neighbors[n] == junction ) continue;
            if ( nedges == NEDGES ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: too many neighboring components when extracting them\n", proc );
              return( -1 );
            }
            edges[nedges] = neighbors[n];
            nedges ++;
          }
        }

        if ( nedges < 2 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: not enough neighboring components\n", proc );
          return( -1 );
        }

        return( nedges );

      }
    }
    break;
  }

  if ( _verbose_ )
    fprintf( stderr, "%s: this point should not be reached\n", proc );
  return( -1 );
}










static int _mergeAtJunction( typeComponentImage *treeImage,
                             typeVoxelTree *iTree, int junction )
{
  char *proc = "_mergeAtJunction";
  voxelTreeIntPoint *newpts = (voxelTreeIntPoint *)NULL;
  int n_data, n_allocated_data = 0;
  int p, i, j, k;
  size_t indp, indn;
  int junctionhasbeenadded;
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  int neighbors[27];
  int x, y, z, n;
  int nedges, edges[NEDGES];
  int minlabel, maxlabel;

  /* get junction neighbors, ie the components to be merged
   */
  for ( n=0; n<NEDGES; n++ ) edges[n] = 0;
  nedges = _getJunctionNeighbors( treeImage, iTree, junction, edges );

  if ( nedges != 2 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird case when for junction #%d, nedges should be 2 but is %d\n",
               proc, junction, nedges );
    return( -1 );
  }

  if ( _debug_ ) {
    fprintf( stderr, "     component %d is connected to %d and %d\n",
             edges[0],
             iTree->data[edges[0]].indexFirstJunction,
             iTree->data[edges[0]].indexLastJunction );
    fprintf( stderr, "     component %d is connected to %d and %d\n",
             edges[1],
             iTree->data[edges[1]].indexFirstJunction,
             iTree->data[edges[1]].indexLastJunction );
  }

  /* merge the junction with the two neighboring components
   */
  if ( edges[0] < edges[1] ) {
    minlabel = edges[0];
    maxlabel = edges[1];
  }
  else {
    minlabel = edges[1];
    maxlabel = edges[0];
  }

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "\t   will merge edges " );
    fprintf( stderr, "%d and %d ", minlabel, maxlabel );
    fprintf( stderr, "with junction %d\n", junction );
  }


  /* allocate array of points
   */
  n_allocated_data = iTree->data[minlabel].n_data + iTree->data[junction].n_data + iTree->data[maxlabel].n_data;
  newpts = (voxelTreeIntPoint*)vtmalloc( n_allocated_data * sizeof(voxelTreeIntPoint), "newpts", proc );
  if ( newpts == (voxelTreeIntPoint*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }


  switch ( treeImage->componentBufferType ) {
  default :
    vtfree( newpts );
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *resBuf = (u16*)(treeImage->componentBuf);

      /* copy points
       */
      n_data = 0;

      /* copy points of iTree->data[minlabel]
       */
      if ( iTree->data[minlabel].indexFirstJunction == junction ) {
        for ( p=iTree->data[minlabel].n_data-1; p>=0; p-- ) {
          i = iTree->data[minlabel].data[p].x;
          j = iTree->data[minlabel].data[p].y;
          k = iTree->data[minlabel].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* the point has already been deleted
           */
          if ( resBuf[indp] == 0 ) continue;
          newpts[n_data] = iTree->data[minlabel].data[p];
          n_data++;
        }
        iTree->data[minlabel].indexFirstJunction = iTree->data[minlabel].indexLastJunction;
      }
      else if ( iTree->data[minlabel].indexLastJunction == junction ) {
        for ( p=0; p<iTree->data[minlabel].n_data; p++ ) {
          i = iTree->data[minlabel].data[p].x;
          j = iTree->data[minlabel].data[p].y;
          k = iTree->data[minlabel].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* the point has already been deleted
           */
          if ( resBuf[indp] == 0 ) continue;
          newpts[n_data] = iTree->data[minlabel].data[p];
          n_data++;
        }
      }
      else {
        vtfree( newpts );
        if ( _verbose_ )
          fprintf( stderr, "%s: this should not occur\n", proc );
        return( -1 );
      }

      /* copy points of junction
       * here each point has 2 neighbors
       * find the point
       */
      do {
        for ( p=0, junctionhasbeenadded=0; p<iTree->data[junction].n_data; p++ ) {
          i = iTree->data[junction].data[p].x;
          j = iTree->data[junction].data[p].y;
          k = iTree->data[junction].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* the point has already been deleted
           * or changed into minlabel
           */
          if ( resBuf[indp] != junction ) continue;

          for ( n=0, z=k-1; z<=k+1; z++ )
          for ( y=j-1; y<=j+1; y++ )
          for ( x=i-1; x<=i+1; x++, n++ ) {
           if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
             neighbors[n] = 0;
             continue;
           }
           indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
           neighbors[n] = resBuf[indn];
          }
          for ( n=0; n<27 && junctionhasbeenadded==0;n++ ) {
            if ( neighbors[n] == minlabel ) {
              newpts[n_data] = iTree->data[junction].data[p];
              n_data++;
              resBuf[indp] = minlabel;
              junctionhasbeenadded = 1;
            }
          }
        }
      } while( junctionhasbeenadded == 1 );

      iTree->data[junction].type = _INTERN_TREE_MERGED_COMPONENT_;

      /* copy points of iTree->data[maxlabel]
       */
      if ( iTree->data[maxlabel].indexFirstJunction == junction ) {
        for ( p=0; p<iTree->data[maxlabel].n_data; p++ ) {
          i = iTree->data[maxlabel].data[p].x;
          j = iTree->data[maxlabel].data[p].y;
          k = iTree->data[maxlabel].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* the point has already been deleted
           */
          if ( resBuf[indp] == 0 ) continue;
          newpts[n_data] = iTree->data[maxlabel].data[p];
          n_data++;
          resBuf[indp] = minlabel;
        }
        iTree->data[minlabel].indexLastJunction = iTree->data[maxlabel].indexLastJunction;
        if ( iTree->data[maxlabel].indexLastJunction == -1 && iTree->data[maxlabel].canbepruned == 0 )
          iTree->data[minlabel].canbepruned = 0;
      }
      else if ( iTree->data[maxlabel].indexLastJunction == junction ) {
        for ( p=iTree->data[maxlabel].n_data-1; p>=0; p-- ) {
          i = iTree->data[maxlabel].data[p].x;
          j = iTree->data[maxlabel].data[p].y;
          k = iTree->data[maxlabel].data[p].z;
          indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
          /* the point has already been deleted
           */
          if ( resBuf[indp] == 0 ) continue;
          newpts[n_data] = iTree->data[maxlabel].data[p];
          n_data++;
          resBuf[indp] = minlabel;
        }
        iTree->data[minlabel].indexLastJunction = iTree->data[maxlabel].indexFirstJunction;
        if ( iTree->data[maxlabel].indexFirstJunction == -1 && iTree->data[maxlabel].canbepruned == 0 )
          iTree->data[minlabel].canbepruned = 0;
      }
      else {
        vtfree( newpts );
        if ( _verbose_ )
          fprintf( stderr, "%s: this should not occur\n", proc );
        return( -1 );
      }

      iTree->data[maxlabel].type = _INTERN_TREE_MERGED_COMPONENT_;

    }
    break;
  }

  if ( _debug_ ) {
    fprintf( stderr, "     component %d is now connected to %d and %d\n",
             minlabel,
             iTree->data[minlabel].indexFirstJunction,
             iTree->data[minlabel].indexLastJunction );
  }

  if ( iTree->data[minlabel].indexFirstJunction == -1
       || iTree->data[minlabel].indexLastJunction == -1 )
    iTree->data[minlabel].type = _INTERN_TREE_END_EDGE_;

  vtfree( iTree->data[minlabel].data );
  iTree->data[minlabel].data = newpts;
  iTree->data[minlabel].n_allocated_data = n_allocated_data;
  iTree->data[minlabel].n_data = n_data;

  lengthEdgeVoxelTree( iTree, minlabel );

  return( 1 );
}








/**************************************************
 *
 *
 *
 **************************************************/







static int _countNeighborsMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree, int junction,
                                    u8 *resBuf, int *jctDim )
{
  char *proc = "_countNeighborsMixedJunction";
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  int p;
  int i, j, k;
  int x, y, z;
  int n;
  size_t indn, indp, indr;

  for ( n=0, z=0; z<jctDim[2]; z++ )
  for ( y=0; y<jctDim[1]; y++ )
  for ( x=0; x<jctDim[0]; x++, n++ ) {
    resBuf[n] = 0;
  }

  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *theBuf = (u16*)(treeImage->componentBuf);

      for ( p=0; p<iTree->data[junction].n_data; p++ ) {
        i = iTree->data[junction].data[p].x;
        j = iTree->data[junction].data[p].y;
        k = iTree->data[junction].data[p].z;
        indp = ((size_t)k * (size_t)dimy + (size_t)j) * (size_t)dimx + (size_t)i;
        if ( theBuf[indp] == 0 ) continue;

        for ( n=0, z=k-1; z<=k+1; z++ )
        for ( y=j-1; y<=j+1; y++ )
        for ( x=i-1; x<=i+1; x++ ) {
          if ( z < 0 || z >= dimz || y < 0 || y >=dimy || x < 0 || x >= dimx ) {
            continue;
          }
          indn = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;
          if ( theBuf[indn] > 0 ) n++;
        }

        indr = ((size_t)(k-iTree->data[junction].minCorner[2])
                * (size_t)jctDim[1] + (size_t)(j-iTree->data[junction].minCorner[1])) * (size_t)jctDim[0]
                + (size_t)(i-iTree->data[junction].minCorner[0]);
        resBuf[indr] = n;
      }
    }
    break;
  }

  return( 1 );
}







static int _getNeighborsEdgeMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree, int junction,
                                          u8 *ccBuf, int* ccDim, int label,
                                          int *neighbors )
{
  char *proc ="_getNeighborsEdgeMixedJunction";
  int i, j, k, l, npts;
  int x, y, z;
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  int dimz = treeImage->theDim[2];
  size_t ind, indn;
  int a, b, c;
  int f, n, nneighbors = 0;

  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *theBuf = (u16*)(treeImage->componentBuf);

      for ( npts=0, l=0, k=0; k<ccDim[2]; k++ )
      for ( j=0; j<ccDim[1]; j++ )
      for ( i=0; i<ccDim[0]; i++, l++ ) {
        if ( ccBuf[l] != label ) continue;
        npts ++;

        x = i + iTree->data[junction].minCorner[0];
        y = j + iTree->data[junction].minCorner[1];
        z = k + iTree->data[junction].minCorner[2];
        ind = ((size_t)z * (size_t)dimy + (size_t)y) * (size_t)dimx + (size_t)x;

        if ( theBuf[ind] != junction ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: weird label in tree image at (%d,%d,%d), should be %d, is %d\n", proc, x, y, z, junction, theBuf[ind] );
          return( -1 );
        }

        /* for neighbors
         * - check whether it is outside the component
         * - check whether it has been erased
         * - check whether it is background
         * - check whether the component has been found
         */
        for ( c=-1; c<=1; c++ )
        for ( b=-1; b<=1; b++ )
        for ( a=-1; a<=1; a++ ) {
          if ( a == 0 && b == 0 && c == 0 ) continue;
          if ( 0 <= i+a && i+a < ccDim[0] && 0 <= j+b && j+b < ccDim[1] && 0 <= k+c && k+c < ccDim[2] ) {
            indn = ( (k+c) * ccDim[1] + (j+b) ) * ccDim[0] + (i+a);
            if ( ccBuf[indn] == label ) continue;
          }
          if ( 0 <= x+a && x+a < dimx && 0 <= y+b && y+b < dimy && 0 <= z+c && z+c < dimz ) {
            ind = ((size_t)(z+c) * (size_t)dimy + (size_t)(y+b)) * (size_t)dimx + (size_t)(x+a);
            if ( theBuf[ind] == 0 ) continue;
            for ( f=0, n=0; n<nneighbors; n++ ) {
              if ( neighbors[n] == theBuf[ind] ) f = 1;
            }
            if ( f == 1 ) continue;
            neighbors[ nneighbors ] = theBuf[ind];
            nneighbors ++;
          }
        }

      }
    }
    break;
  }

  if ( npts == 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: found no points for component %d\n", proc, label );
    return( -1 );
  }

  return( nneighbors );
}





static int _relabelLocalComponentInVoxelTree( typeComponentImage *treeImage, int newLabel,
                                              typeVoxelTree *iTree, int junctionLabel,
                                              u8* labelBuf, int* labelDim, int label )
{
  char *proc = "_relabelLocalComponentInVoxelTree";
  int i, j, k, p;
  voxelTreeIntPoint pt;
  int indl;
  int dimx = treeImage->theDim[0];
  int dimy = treeImage->theDim[1];
  size_t ind;
  int n;

  switch ( treeImage->componentBufferType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such tree image type not handled yet\n", proc );
    return( -1 );

  case USHORT :
    {
      u16 *theBuf = (u16*)(treeImage->componentBuf);

      for ( n=0, p=0; p<iTree->data[junctionLabel].n_data; p++ ) {
        /* junction point
         */
        pt = iTree->data[junctionLabel].data[p];
        /* point in sub-image
         */
        i = pt.x - iTree->data[junctionLabel].minCorner[0];
        j = pt.y - iTree->data[junctionLabel].minCorner[1];
        k = pt.z - iTree->data[junctionLabel].minCorner[2];
        indl = (k * labelDim[1] + j) * labelDim[0] + i;
        /* get a point to be processed
         * - add point to new component
         * - update value in component image
         * - remove point from previous junction
         */
        if ( labelBuf[indl] != label ) continue;
        ind = ((size_t)pt.z * (size_t)dimy + (size_t)pt.y) * (size_t)dimx + (size_t)pt.x;
        theBuf[ ind ] = newLabel;
        n++;
      }
    }
    break;
  }

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "\t   - has relabel %d points of local component %d into %d\n",
             n, label, newLabel );
  }

  return( 1 );
}





static int _processJunctionsInMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree,
                                             int junctionLabel,
                                             u8* jctBuf, int* jctDim, int njunctions,
                                             int *labels, int *nlabels )
{
  char *proc = "_processJunctionsInMixedJunction";
  int l;
  int newLabel;
  voxelTreeComponent *c;
  voxelTreeIntPoint pt;
  int i, j, k, p;
  int indj;

  if ( njunctions <= 1 ) {
    return( 1 );
  }


  /* process new junctions
   */
  for ( l=2; l<=njunctions; l++ ) {

    /* get a new label
     */
    newLabel = getNewComponentFromVoxelTree( iTree );
    labels[*nlabels] = newLabel;
    (*nlabels) ++;

    /* relabel the component image
     */
    if ( _relabelLocalComponentInVoxelTree( treeImage, newLabel,
                                            iTree, junctionLabel,
                                            jctBuf, jctDim, l ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to relabel local junction %d\n", proc, l );
      return( -1 );
    }

    /* fill the component
     */
    c = &(iTree->data[newLabel]);
    c->type = _INTERN_TREE_JUNCTION_;

    if ( _verbose_ >= 3 ) {
      fprintf( stderr, "\t     create junction component %d\n", newLabel );
    }

    /* process points
     */
    for ( p=0; p<iTree->data[junctionLabel].n_data; p++ ) {
      /* junction point
       */
      pt = iTree->data[junctionLabel].data[p];
      /* point in sub-image
       */
      i = pt.x - iTree->data[junctionLabel].minCorner[0];
      j = pt.y - iTree->data[junctionLabel].minCorner[1];
      k = pt.z - iTree->data[junctionLabel].minCorner[2];
      indj = (k * jctDim[1] + j) * jctDim[0] + i;
      /* get a point to be processed
       * - add point to new component
       * - remove point from previous junction
       */
      if ( jctBuf[indj] != l ) continue;
      if ( addPointToVoxelTreeComponent( c, &pt ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to add point to new component %d\n", proc, newLabel );
        return( -1 );
      }
      if ( p < iTree->data[junctionLabel].n_data-1 ) {
        pt = iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ];
        iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ] = iTree->data[junctionLabel].data[ p ];
        iTree->data[junctionLabel].data[ p ] = pt;
        p--;
      }
      iTree->data[junctionLabel].n_data --;
    }

    /* update junction
     */
    if ( voxelTreeJunctionCenter( &(iTree->data[newLabel]) ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to compute center of junction #%d\n", proc, newLabel );
        return( -1 );
      }
    }
  }

  return( 1 );
}





static int _mergeEdgeInMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree,
                                      int edgeLabel, int junctionLabel,
                                      u8* edgBuf, int* edgDim, int e )
{
  char *proc = "_mergeEdgeInMixedJunction";
  int p;
  int i, j, k;
  int indj;
  voxelTreeComponent *c;
  voxelTreeIntPoint pt;
  voxelTreeComponent newComponent;


  /* relabel the component image
   * -> relabel junctionLabel from edgBuf to edgeLabel in treeImage
   */
  if ( _relabelLocalComponentInVoxelTree( treeImage, edgeLabel,
                                          iTree, junctionLabel,
                                          edgBuf, edgDim, e ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to relabel local edge %d\n", proc, e );
    return( -1 );
  }

  /* remove points from junction
   */
  for ( p=0; p<iTree->data[junctionLabel].n_data; p++ ) {
    /* junction point
     */
    pt = iTree->data[junctionLabel].data[p];
    /* point in sub-image
     */
    i = pt.x - iTree->data[junctionLabel].minCorner[0];
    j = pt.y - iTree->data[junctionLabel].minCorner[1];
    k = pt.z - iTree->data[junctionLabel].minCorner[2];
    indj = (k * edgDim[1] + j) * edgDim[0] + i;
    /* get a point to be processed
     * - remove point from previous junction
     */
    if ( edgBuf[indj] != e ) continue;
    if ( p < iTree->data[junctionLabel].n_data-1 ) {
      pt = iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ];
      iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ] = iTree->data[junctionLabel].data[ p ];
      iTree->data[junctionLabel].data[ p ] = pt;
      p--;
    }
    iTree->data[junctionLabel].n_data --;
  }

  /* get the component from the other end
   * it gets the junction labels
   * but does not update the length
   */
  c = &(iTree->data[ edgeLabel ]);
  initVoxelTreeComponent( &newComponent );
  newComponent.type = _INTERN_TREE_EDGE_;

  if ( c->indexFirstJunction == junctionLabel ) {
    if ( extractEdgeVoxelTree( treeImage, &newComponent, &(c->data[c->n_data-1]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to extract new edge\n", proc );
      return( -1 );
    }
    newComponent.infoFirstPoint = c->infoLastPoint;
    newComponent.infoLastPoint = c->infoFirstPoint;
  }
  else if ( c->indexLastJunction == junctionLabel ) {
    if ( extractEdgeVoxelTree( treeImage, &newComponent, &(c->data[0]) ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to extract new edge\n", proc );
      return( -1 );
    }
    newComponent.infoFirstPoint = c->infoFirstPoint;
    newComponent.infoLastPoint = c->infoLastPoint;
  }

  freeVoxelTreeComponent( &(iTree->data[ edgeLabel ]) );
  iTree->data[ edgeLabel ] = newComponent;

  /* update edge length
   */
  lengthEdgeVoxelTree( iTree, edgeLabel );

  return( 1 );
}





static int _processEdgesInMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree,
                                         int junctionLabel,
                                         u8* edgBuf, int* edgDim, int nedges )
{
  char *proc = "_processEdgesInMixedJunction";
  int e;
  int nneighbors, neighbors[NEDGES];
  int newLabel;
  voxelTreeComponent *c;
  voxelTreeIntPoint pt;
  int i, j, k, p;
  int n, ind;


  if ( nedges <= 0 ) return( 1 );


  for ( e=1; e<=nedges; e++ ) {
    /* get neighbor label for edges
     * there should be two labels per edge
     * since new junctions have been relabeled
     */
    for ( i=0; i<NEDGES; i++ ) neighbors[i] = 0;
    nneighbors = _getNeighborsEdgeMixedJunction( treeImage, iTree, junctionLabel,
                                                 edgBuf, edgDim, e, neighbors);
    if ( _verbose_ >= 3 ) {
      fprintf( stderr, "\t   - local edge %2d, found neighbors: ", e );
      for ( i=0; i<nneighbors; i++ )
        fprintf( stderr, "%d ", neighbors[i] );
      fprintf( stderr, "\n" );
    }

    if ( nneighbors != 2 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird number of neighbors (%d) for local edge %d\n", proc, nneighbors, e );
      return( -1 );
    }



    if ( ( iTree->data[ neighbors[0] ].type == _INTERN_TREE_INNER_EDGE_
           || iTree->data[ neighbors[0] ].type == _INTERN_TREE_END_EDGE_
           || iTree->data[ neighbors[0] ].type == _INTERN_TREE_EDGE_ )
         && iTree->data[ neighbors[1] ].type == _INTERN_TREE_JUNCTION_ ) {
      if ( _verbose_ >= 3 ) {
        fprintf( stderr, "\t     merge local edge %2d with edge %3d\n", e, neighbors[0] );
      }
      if ( _mergeEdgeInMixedJunction( treeImage, iTree, neighbors[0], junctionLabel,
                                      edgBuf, edgDim, e ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to merge local edge %d with edge %d\n", proc, e, neighbors[0] );
        return( -1 );
      }
    }

    else if ( ( iTree->data[ neighbors[1] ].type == _INTERN_TREE_INNER_EDGE_
                || iTree->data[ neighbors[1] ].type == _INTERN_TREE_END_EDGE_
                || iTree->data[ neighbors[1] ].type == _INTERN_TREE_EDGE_ )
              && iTree->data[ neighbors[0] ].type == _INTERN_TREE_JUNCTION_ ) {
      if ( _verbose_ >= 3 ) {
        fprintf( stderr, "\t     merge local edge %2d with edge %3d\n", e, neighbors[1] );
      }
      if ( _mergeEdgeInMixedJunction( treeImage, iTree, neighbors[1], junctionLabel,
                                      edgBuf, edgDim, e ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to merge local edge %d with edge %d\n", proc, e, neighbors[1] );
        return( -1 );
      }
    }

    else if ( iTree->data[ neighbors[0] ].type == _INTERN_TREE_JUNCTION_
              && iTree->data[ neighbors[1] ].type == _INTERN_TREE_JUNCTION_ ) {

      /* count points
       */
      for ( n=0, ind=0, k=0; k<edgDim[2]; k++ )
      for ( j=0; j<edgDim[1]; j++ )
      for ( i=0; i<edgDim[0]; i++, ind++ ) {
        if ( edgBuf[ind] == e ) n++;
      }

      /* get a new label
       */
      newLabel = getNewComponentFromVoxelTree( iTree );

      if ( _verbose_ >= 3 ) {
        fprintf( stderr, "\t     create edge component %d with %d points\n", newLabel, n );
      }

      /* relabel the component image
       */
      if ( _relabelLocalComponentInVoxelTree( treeImage, newLabel,
                                              iTree, junctionLabel,
                                              edgBuf, edgDim, e ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to relabel local edge %d\n", proc, e );
        return( -1 );
      }

      /* fill the component
       */
      c = &(iTree->data[newLabel]);
      c->type = _INTERN_TREE_EDGE_;



      if ( n >= 3 ) {
        /* here we have to recognize an end point
         * and to track the edge
         */
        if ( _verbose_ )
          fprintf( stderr, "%s: case to be developed\n", proc );
        return( -1 );
      }


      /* process points
       */
      for ( p=0; p<iTree->data[junctionLabel].n_data; p++ ) {
        /* junction point
         */
        pt = iTree->data[junctionLabel].data[p];
        /* point in sub-image
         */
        i = pt.x - iTree->data[junctionLabel].minCorner[0];
        j = pt.y - iTree->data[junctionLabel].minCorner[1];
        k = pt.z - iTree->data[junctionLabel].minCorner[2];
        ind = (k * edgDim[1] + j) * edgDim[0] + i;
        /* get a point to be processed
         * - add point to new component
         * - remove point from previous junction
         */
        if ( edgBuf[ind] != e ) continue;
        if ( addPointToVoxelTreeComponent( c, &pt ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to add point to new component %d\n", proc, newLabel );
          return( -1 );
        }
        if ( p < iTree->data[junctionLabel].n_data-1 ) {
          pt = iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ];
          iTree->data[junctionLabel].data[ iTree->data[junctionLabel].n_data-1 ] = iTree->data[junctionLabel].data[ p ];
          iTree->data[junctionLabel].data[ p ] = pt;
          p--;
        }
        iTree->data[junctionLabel].n_data --;
      }

      /* update edge
       * - connecting junction
       */
      if ( extractJunctionEdgeVoxelTree( treeImage, &(iTree->data[newLabel]) ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to compute center of junction #%d\n", proc, newLabel );
          return( -1 );
        }
      }
    }

    else {
      if ( _verbose_ )
        fprintf( stderr, "%s: weird configuration for local edge %d\n", proc, e );
      return( -1 );
    }
  }


  return( 1 );
}





static int _processMixedJunction( typeComponentImage *treeImage, typeVoxelTree *iTree, int junction )
{
  char *proc = "_processMixedJunction";
  int junctionDim[3];
  int i, v, j;
  int nedges, njunctions;
  u8 *nghBuf, *tmpBuf, *edgeBuf, *juncBuf;
  int newJunctionLabel[NEDGES], nNewJunctionLabel=0;
  voxelTreeComponent *c;
  int update;


  for ( i=0; i<3; i++ ) {
    junctionDim[i] = iTree->data[junction].maxCorner[i] - iTree->data[junction].minCorner[i] + 1;
  }
  v = junctionDim[0] * junctionDim[1] * junctionDim[2];

  tmpBuf = (u8*)vtmalloc( 4*v*sizeof(u8), "tmpBuf", proc );
  if ( tmpBuf == (u8*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: auxiliary buffer allocation failed\n", proc );
    return( -1 );
  }
  nghBuf = edgeBuf = juncBuf = tmpBuf;
  nghBuf  += v;
  edgeBuf += 2*v;
  juncBuf += 3*v;

  /* count neighbors
   * recall that center point is counted too
   */
  if ( _countNeighborsMixedJunction( treeImage, iTree, junction, nghBuf, junctionDim ) != 1 ) {
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to count neighbors\n", proc );
    return( -1 );
  }

  /* count edges
   */
  for ( i=0; i<v; i++ ) {
    tmpBuf[i] = ( nghBuf[i] == 1 || nghBuf[i] == 2 || nghBuf[i] == 3 ) ? 255 : 0;
  }
  nedges = CountConnectedComponents( tmpBuf, UCHAR, edgeBuf, UCHAR, junctionDim );

  /* count junctions
   */
  for ( i=0; i<v; i++ ) {
    tmpBuf[i] = ( nghBuf[i] > 3 ) ? 255 : 0;
  }
  njunctions = CountConnectedComponents( tmpBuf, UCHAR, juncBuf, UCHAR, junctionDim );


  /* summary
   */
  if ( _verbose_ >= 3 ) {
     fprintf( stderr, "\t   found %d tree branches and %d junctions\n", nedges, njunctions );
  }



  /* njunctions == 0
   */
  if ( njunctions == 0 || nedges == 0 ) {
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: weird number of junctions (%d) or edges (%d)\n", proc, njunctions, nedges );
    return( -1 );
  }




  /* 1. split junctions into components, keep junction labels
   *    -> add junctions
   *    junctionsLabels: array of junction labels
   *    to further update incident edges
   * Points have been removed from 'junction'
   * recall to update its center
   */
  if ( njunctions >= 2 ) {
    if ( _processJunctionsInMixedJunction( treeImage, iTree, junction,
                                           juncBuf, junctionDim, njunctions, newJunctionLabel, &nNewJunctionLabel ) != 1 ) {
      vtfree( tmpBuf );
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to relabel %d junction sub-components of junction %d\n", proc, njunctions, junction );
      return( -1 );
    }
  }

  /* 2. edges which had a neighbor into edge
   *    -> add to edge
   * 3. edges which had neighbors into junctions
   *    -> add edge
   * 4. recognize label end for edges, update length
   */
  if ( _processEdgesInMixedJunction( treeImage, iTree, junction,
                                     edgeBuf, junctionDim, nedges ) != 1 ) {
    vtfree( tmpBuf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to relabel %d edge sub-components of junction %d\n", proc, nedges, junction );
    return( -1 );
  }

  vtfree( tmpBuf );


  /* update junction center
   * recall that points have been removed/transfered to other component
   */
  if ( voxelTreeJunctionCenter( &(iTree->data[junction]) ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to compute center of junction #%d\n", proc, junction );
      return( -1 );
    }
  }


  /* update edges connected to junction or to newly created junctions
   * - connecting junctions (junction may have been splitted)
   * - length (since junctions center have been updated)
   */
  for ( i=1; i<iTree->n_data; i++ ) {
    c = &(iTree->data[i]);
    if ( c->type != _INTERN_TREE_INNER_EDGE_
         && c->type != _INTERN_TREE_END_EDGE_
         && c->type != _INTERN_TREE_EDGE_ )
        continue;
    update = 0;
    if ( c->indexFirstJunction == junction || c->indexLastJunction == junction ) update = 1;
    for ( j=0; j<nNewJunctionLabel; j++ )
      if ( c->indexFirstJunction == newJunctionLabel[j] || c->indexLastJunction == newJunctionLabel[j] ) update = 1;
    if ( update == 0 ) continue;
    if ( _verbose_ >= 4 ) {
      fprintf( stderr, "\t     update edge %3d (connected to %3d and %3d), length=%f\n", i, c->indexFirstJunction, c->indexLastJunction, c->length );
    }
    if ( extractJunctionEdgeVoxelTree( treeImage, c ) != 1 ) {
      if ( _verbose_ ) {
        fprintf( stderr, "%s: unable to update connecting junctions of component #%d\n", proc, i );
        return( -1 );
      }
    }
    lengthEdgeVoxelTree( iTree, i );
    if ( _verbose_ >= 4 ) {
      fprintf( stderr, "\t                  to (connected to %3d and %3d), length=%f\n", c->indexFirstJunction, c->indexLastJunction, c->length );
    }

  }

  return( 1 );
}










/**************************************************
 *
 *
 *
 **************************************************/



static int _removeComponent( typeComponentImage *treeImage,
                       typeVoxelTree *iTree, int label )
{
  char *proc = "_removeComponent";
  int ret;
  int junctionRemovedPoints;

  int junction;
  _enumJunctionFuture junctionFuture;



  /* check whether it is an end edge
   * and it has only one end (and thus only one neighbor)
   */

  if ( iTree->data[label].type != _INTERN_TREE_END_EDGE_ ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird type for component to be removed\n", proc );
    return( -1 );
  }

  if ( (iTree->data[label].indexFirstJunction > 0 && iTree->data[label].indexLastJunction > 0)
       || (iTree->data[label].indexFirstJunction <= 0 && iTree->data[label].indexLastJunction <= 0) ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird neighbor configuration for component to be removed\n", proc );
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "\t   component %d is connected to %d and %d\n",
               label, iTree->data[label].indexFirstJunction,
               iTree->data[label].indexLastJunction );
    }
    return( -1 );
  }



  /* erase component
   * return if the end edge has not been completely removed
   */
  ret = _eraseEndEdge( treeImage, iTree, label);
  if ( ret == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when erasing component #%d\n", proc, label );
    return( -1 );
  }
  else if ( ret == 0 ) {
    /* the component can not be entirely deleted
     * some anchor points prevent from it
     */
    iTree->data[label].canbepruned = 0;
    return( 1 );
  }



  /* here the component has been entirely deleted
   */
  iTree->data[label].type = _INTERN_TREE_REMOVED_COMPONENT_;



  /* find the connected junction
   * - junction may be -1 if there is a single branch to be removed
   */
  junction = 0;
  if ( iTree->data[label].indexFirstJunction <= 0 )
    junction = iTree->data[label].indexLastJunction;
  else if ( iTree->data[label].indexLastJunction <= 0 )
    junction = iTree->data[label].indexFirstJunction;
  else {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird neighbor configuration for component to be removed\n", proc );
    return( -1 );
  }
  if ( junction <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no junction for component to be removed ?!\n", proc );
    return( -1 );
  }

  if ( _verbose_ >= 3 ) {
    fprintf( stderr, "\t - process junction %d in [%d %d %d]x[%d %d %d]\n",
             junction, iTree->data[junction].minCorner[0], iTree->data[junction].minCorner[1],
             iTree->data[junction].minCorner[2], iTree->data[junction].maxCorner[0],
             iTree->data[junction].maxCorner[1], iTree->data[junction].maxCorner[2] );
  }


  /* remove simple points from the junction
   */
  junctionRemovedPoints = _eraseJunction( treeImage, iTree, junction );
  if ( junctionRemovedPoints == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when erasing junction %d\n", proc, junction );
    return( -1 );
  }


  /* see what must be done now
   */
  junctionFuture = _getJunctionFuture( treeImage, iTree, junction );
  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\t   junction future is " );
    switch( junctionFuture ) {
    default : fprintf( stderr, "unknown, how embarassing\n" ); break;
    case _FUTURE_UNKNOWN_ :     fprintf( stderr, "_FUTURE_UNKNOWN_\n" ); break;
    case _FUTURE_IS_JUNCTION_ : fprintf( stderr, "_FUTURE_IS_JUNCTION_\n" ); break;
    case _FUTURE_IS_EDGE_ :     fprintf( stderr, "_FUTURE_IS_EDGE_\n" ); break;
    case _FUTURE_IS_MIXED_ :    fprintf( stderr, "_FUTURE_IS_MIXED_\n" ); break;
    }
  }



  switch( junctionFuture ) {
  default :
  case _FUTURE_UNKNOWN_ :
    if ( _verbose_ )
      fprintf( stderr, "%s: unknown junction future for junction %d\n", proc, junction );
    return( -1 );
    break;
  case _FUTURE_IS_JUNCTION_ :
    if ( _verbose_ >= 3 )
      fprintf( stderr, "\t   junction %d is still a junction\n", junction );
    /* still a junction
     * update
     * - junction center
     * - edge length
     */
    if ( junctionRemovedPoints > 0 ) {
      if ( _verbose_ >= 3 )
        fprintf( stderr, "\t   update junction %d\n", junction );
      if ( voxelTreeJunctionCenter( &(iTree->data[junction]) ) != 1 ) {
        if ( _verbose_ ) {
          fprintf( stderr, "%s: unable to update center of junction #%d\n", proc, junction );
        }
        return( -1 );
      }
      _updateIncidentEdges( iTree, junction );
    }
    /* still a junction, do nothing
     */
    break;
  case _FUTURE_IS_EDGE_ :
    if ( _mergeAtJunction( treeImage, iTree, junction ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when merging components through junction %d\n", proc,
                 junction );
      return( -1 );
    }
    break;
  case _FUTURE_IS_MIXED_ :
    if ( _verbose_ >= 3 )
      fprintf( stderr, "\t   process mixed junction %d (%d points)\n", junction, iTree->data[junction].n_data );
    if ( _processMixedJunction( treeImage, iTree, junction ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when processing mixed junction #%d\n", proc, junction );
      return( -1 );
    }
    break;
  }

  return( 1 );
}





/**************************************************
 *
 *
 *
 **************************************************/





int pruneVoxelTreeWithOrder( typeComponentImage *treeImage, typeVoxelTree *iTree,
                             int maxorder )
{
  char *proc = "pruneVoxelTreeWithOrder";
  int i;
  voxelTreeComponent *c;
  int label;
  int iteration = 0;


  if ( maxorder <= 1 )
    return( 1 );

  /* before pruning
   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s\n", proc );
    fprintfSummaryVoxelTree( stderr, iTree, "before pruning" );
  }

  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\n" );
    fprintVoxelTree( stderr, iTree, 0 );
  }


  /* prune
   */
  do {
    for ( i=1, label=-1; i<iTree->n_data && label==-1; i++ ) {
      c = &(iTree->data[i]);
      if ( c->type != _INTERN_TREE_END_EDGE_ )
        continue;
      if ( c->infoFirstPoint.order > maxorder || c->infoLastPoint.order > maxorder )
        label = i;
    }

    if ( label == -1 )
        continue;

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: remove label %d at iteration %d\n", proc, label, iteration );
      fprintf( stderr, "\t real length = %f, nb of voxel = %d",
               iTree->data[label].length,
               iTree->data[label].n_data );
      fprintf( stderr, "\n" );
    }

    if ( _removeComponent( treeImage, iTree, label ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when removing component %d\n", proc, label );
      return( -1 );
    }

    iteration ++;

  } while ( label > 0 );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s\n", proc );
    fprintfSummaryVoxelTree( stderr, iTree, "after pruning" );
  }

  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\n" );
    fprintVoxelTree( stderr, iTree, 0 );
  }

  return( 1 );
}





int pruneVoxelTreeWithEdgeLength( typeComponentImage *treeImage, typeVoxelTree *iTree,
                    float realBranchLength, int voxelBranchLength,
                    int nmaxendedges )
{
  char *proc = "pruneVoxelTreeWithEdgeLength";
  enumSmallestEdgeTest test = _CHECK_REAL_LENGTH_FIRST_;
  int i;
  int nendedges;
  int label;
  int iteration = 0;


  if ( realBranchLength <= 0.0 && voxelBranchLength < 0 && nmaxendedges <= 0 ) {
    return( 1 );
  }
  if ( voxelBranchLength > 0 && realBranchLength <= 0.0 )
      test = _CHECK_VOXEL_LENGTH_FIRST_;


  /* before pruning
   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s\n", proc );
    fprintfSummaryVoxelTree( stderr, iTree, "before pruning" );
  }

  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\n" );
    fprintVoxelTree( stderr, iTree, 0 );
  }

  for ( nendedges=0, i=0; i<iTree->n_data; i++ ) {
    if ( iTree->data[i].type == _INTERN_TREE_END_EDGE_ ) {
        nendedges++;
    }
  }


  /* prune
   */
  do {

    label = _getSmallestEndEdge( iTree, test );
    if ( label == -1 )
        continue;

    /* should we delete the edge ?
     */
    if ( (nmaxendedges > 0 && nendedges > nmaxendedges)
         || (realBranchLength > 0.0 && iTree->data[label].length <= realBranchLength)
         || (voxelBranchLength > 0 && iTree->data[label].n_data <= voxelBranchLength) ) {
        ;
    }
    else {
        label = -1;
        continue;
    }

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "%s: remove label %d at iteration %d\n", proc, label, iteration );
      fprintf( stderr, "\t real length = %f, nb of voxel = %d",
               iTree->data[label].length,
               iTree->data[label].n_data );
      if ( nmaxendedges > 0 )
          fprintf( stderr, ", #end edges = %d/%d",
                   nendedges, nmaxendedges );
      fprintf( stderr, "\n" );
    }

    if ( _removeComponent( treeImage, iTree, label ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: error when removing component %d\n", proc, label );
      return( -1 );
    }

    /* update the number of end edges
     * nendedges -- should have been sufficient
     */
    for ( nendedges=0, i=0; i<iTree->n_data; i++ ) {
      if ( iTree->data[i].type == _INTERN_TREE_END_EDGE_ ) {
          nendedges++;
      }
    }


    iteration ++;

  } while ( label > 0 );

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s\n", proc );
    fprintfSummaryVoxelTree( stderr, iTree, "after pruning" );
  }

  if ( _verbose_ >= 4 ) {
    fprintf( stderr, "\n" );
    fprintVoxelTree( stderr, iTree, 0 );
  }

  return( 1 );
}





int pruneTreeImage( void* theBuf, bufferType theType,
                    void* ancBuf, bufferType ancType,
                    void* resBuf, bufferType resType,
                    int *theDim, float *theSize,
                    float realBranchLength, int voxelBranchLength,
                    int nmaxendedges )
{
  char *proc = "pruneTreeImage";
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
      fprintf( stderr, "%s: error when labeling components in tree image\n", proc );
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
      fprintf( stderr, "%s: error when building voxel tree\n", proc );
    return( -1 );
  }

  /* pruning
   */
  if ( pruneVoxelTreeWithEdgeLength( &treeImage, &iTree, realBranchLength, voxelBranchLength,
                       nmaxendedges ) != 1 ) {
    freeVoxelTree( &iTree );
    freeComponentImage( &treeImage );
    if ( _verbose_ )
        fprintf( stderr, "%s: error when pruning voxel tree\n", proc );
      return( -1 );
  }

  /* pruning is done
   */

  freeVoxelTree( &iTree );

  if ( thresholdBuffer( treeImage.componentBuf, treeImage.componentBufferType,
                        resBuf, resType, treeImage.theDim, 1.0 ) != 1 ) {
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when thresholding image tree\n", proc );
    return( -1 );
  }

  freeComponentImage( &treeImage );

  return( 1 );
}





















