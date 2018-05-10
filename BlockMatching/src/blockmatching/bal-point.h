/*************************************************************************
 * bal-point.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 23 sep 2013 16:52:21 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_POINT_H
#define BAL_POINT_H

#ifdef __cplusplus
extern "C" {
#endif





#include <bal-stddef.h>


extern void BAL_SetVerboseInBalPoint( int v );
extern void BAL_IncrementVerboseInBalPoint(  );
extern void BAL_DecrementVerboseInBalPoint(  );

extern void BAL_SetDebugInBalPoint( int v );
extern void BAL_IncrementDebugInBalPoint(  );
extern void BAL_DecrementDebugInBalPoint(  );



typedef struct {
  bal_typeFieldPoint *data;
  int n_data;
  int n_allocated_data;

  /* is the displacement encoded in voxel or real units ?
   */
  enumUnitTransfo unit;
  typeVoxelSize vx;          /* real voxel size in X dimension */
  typeVoxelSize vy;          /* real voxel size in Y dimension */
  typeVoxelSize vz;          /* real voxel size in Z dimension */

} bal_typeFieldPointList;

extern void BAL_InitTypeFieldPointList( bal_typeFieldPointList *l );
extern void BAL_FreeTypeFieldPointList( bal_typeFieldPointList *l );
extern void BAL_PrintTypeFieldPointList( FILE *f, bal_typeFieldPointList *l );
extern int BAL_ReadTypeFieldPointList( bal_typeFieldPointList *l, char *filename );



typedef struct {
  bal_doublePoint *data;
  int n_data;
  int n_allocated_data;

  /* is the displacement encoded in voxel or real units ?
   */
  enumUnitTransfo unit;
  typeVoxelSize vx;          /* real voxel size in X dimension */
  typeVoxelSize vy;          /* real voxel size in Y dimension */
  typeVoxelSize vz;          /* real voxel size in Z dimension */

} bal_doublePointList;

extern void BAL_InitDoublePointList( bal_doublePointList *l );
extern int BAL_AllocDoublePointList( bal_doublePointList *l, int s );
extern void BAL_FreeDoublePointList( bal_doublePointList *l );
extern void BAL_PrintDoublePointList( FILE *f, bal_doublePointList *l );
extern int BAL_ReadDoublePointList( bal_doublePointList *l, char *filename );
extern int BAL_WriteDoublePointList( bal_doublePointList *l, char *filename );




#ifdef __cplusplus
}
#endif

#endif
