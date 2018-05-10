/****************************************************
 * topological-voxel-tree-main.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Mar 20 jui 2017 11:38:23 CEST
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

#include <vtmalloc.h>

#include <topological-voxel-tree-main.h>
#include <topological-voxel-tree-order.h>
#include <topological-voxel-tree-pruning.h>







/**************************************************
 *
 *
 *
 **************************************************/


static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInTopologicalVoxelTreeMain( int v )
{
  _verbose_ = v;
}

void incrementVerboseInTopologicalVoxelTreeMain(  )
{
  _verbose_ ++;
}

void decrementVerboseInTopologicalVoxelTreeMain(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInTopologicalVoxelTreeMain( int d )
{
  _debug_ = d;
}

void incrementDebugInTopologicalVoxelTreeMain(  )
{
  _debug_ ++;
}

void decrementDebugInTopologicalVoxelTreeMain(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




















/**************************************************
 *
 *
 **************************************************/

int relabelVoxelTree( typeComponentImage *treeImage,
                             void* ancBuf, bufferType ancType,
                             typeVoxelTree *iTree )
{
  char *proc = "relabelVoxelTree";

  freeVoxelTree( iTree );
  if ( treeImage->component != (typeComponent*)NULL )
    vtfree( treeImage->component );
  treeImage->component = (typeComponent*)NULL;
  treeImage->n_component = 0;
  treeImage->n_allocated_component = 0;

  if ( imageToComponentImage( treeImage->componentBuf,
                              treeImage->componentBufferType, treeImage ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when labeling components in tree\n", proc );
    return( -1 );
  }

  initVoxelTree( iTree );
  if ( componentImageToVoxelTree( treeImage, ancBuf, ancType, iTree ) != 1 ) {
    freeVoxelTree( iTree );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building internal tree\n", proc );
    return( -1 );
  }

  if ( ancBuf != (void*)NULL ) {
    if ( branchesVoxelTree( treeImage, iTree ) != 1 ) {
      freeVoxelTree( iTree );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when extracting branches\n", proc );
        return( -1 );
    }
  }

  return( 1 );
}




/**************************************************
 *
 *
 **************************************************/


int buildVoxelTree( void* theBuf, bufferType theType,
                    void* ancBuf, bufferType ancType,
                    int *theDim, float *theSize,
                    typeVoxelTree *iTree,
                    float realBranchLength, int voxelBranchLength, int nmaxendedges,
                    int maxorder )
{
  char *proc = "buildVoxelTree";
  typeComponentImage treeImage;


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
   * internal tree is assumed to be an empty structure
   * all edges can be pruned after voxel tree building
   */

  initVoxelTree( iTree );
  if ( componentImageToVoxelTree( &treeImage, ancBuf, ancType, iTree ) != 1 ) {
    freeVoxelTree( iTree );
    freeComponentImage( &treeImage );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building internal tree\n", proc );
    return( -1 );
  }


  if ( realBranchLength > 0.0 || voxelBranchLength > 0 || nmaxendedges > 0 ) {
    /* pruning based on length
     */
    if ( pruneVoxelTreeWithEdgeLength( &treeImage, iTree,
                                       realBranchLength, voxelBranchLength,
                                       nmaxendedges ) != 1 ) {
      freeVoxelTree( iTree );
      freeComponentImage( &treeImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when pruning voxel tree\n", proc );
      return( -1 );
    }

  }



  /* branches, anchor image is required
   */
  if ( ancBuf != (void*)NULL ) {

    if ( branchesVoxelTree( &treeImage, iTree ) != 1 ) {
      freeVoxelTree( iTree );
      freeComponentImage( &treeImage );
      if ( _verbose_ )
        fprintf( stderr, "%s: error when extracting branches\n", proc );
      return( -1 );
    }

    if ( maxorder >= 2 ) {
      if ( pruneVoxelTreeWithOrder( &treeImage, iTree,
                                         maxorder ) != 1 ) {
        freeVoxelTree( iTree );
        freeComponentImage( &treeImage );
        if ( _verbose_ )
          fprintf( stderr, "%s: error when pruning voxel tree\n", proc );
        return( -1 );
      }
    }

  }

  if ( _debug_ ) {
    fprintfSummaryVoxelTree( stderr, iTree, "after pruning" );
    if ( _debug_ >= 2 )
      fprintVoxelTree( stderr, iTree, 0 );
  }

  /* if there was some pruning
   * relabel components and branches
   */

  if ( realBranchLength > 0.0 || voxelBranchLength > 0 ||
       nmaxendedges > 0 ||
       maxorder >= 2 ) {

    if ( relabelVoxelTree( &treeImage, ancBuf, ancType, iTree ) != 1 ) {
      freeVoxelTree( iTree );
      freeComponentImage( &treeImage );
      if ( _verbose_ )
          fprintf( stderr, "%s: error when relabeling tree\n", proc );
      return( -1 );
    }

  }

  if ( _debug_ ) {
    fprintf( stderr, "%s: end of procedure\n", proc );
    fprintfSummaryVoxelTree( stderr, iTree, "after relabeling" );
    if ( _debug_ >= 2 )
      fprintVoxelTree( stderr, iTree, 0 );
  }


  return( 1 );
}








