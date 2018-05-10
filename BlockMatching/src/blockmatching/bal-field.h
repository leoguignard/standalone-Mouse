/*************************************************************************
 * bal-field.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_FIELD_H
#define BAL_FIELD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-stddef.h>
#include <bal-matrix.h>




/* depending of the context,
   the units can be either real (ie millimetric)
   or in voxel
*/

typedef struct {
  /* origin of the displacement
     it is the center of a block,
     (ie left_corner + (size-1)/2)
  */
  bal_typeFieldPoint origin;

  /* displacement
  */
  bal_typeFieldPoint vector;

  /* similarity of the pairing
   */
  bal_flag valid;   
  typeField rho; 

  /* error wrt the computed transformation
   */
  double error;

} typeScalarWeightedDisplacement;





typedef struct {
  
  /* data related members
   */
  
  typeScalarWeightedDisplacement *data;
  typeScalarWeightedDisplacement **pointer;

  size_t n_allocated_pairs;

  /* computed pairings
   */
  size_t n_computed_pairs;

  /* retained pairings after selection
   */
  size_t n_selected_pairs;



  /* is the displacement encoded in voxel or real units ?
   */
  enumUnitTransfo unit;

  /* voxel size
     since the pairings are computed between images of the same 
     geometry, they have the same voxel size.
     It allows to go from the "real world" to the "voxel world"
   */
  typeVoxelSize vx;          /* real voxel size in X dimension */
  typeVoxelSize vy;          /* real voxel size in Y dimension */
  typeVoxelSize vz;          /* real voxel size in Z dimension */

  bal_imageGeometry geometry;
  /* matrix to transform voxel coordinates into real coordinates
   * real from voxel
   */
  _MATRIX to_real;

  /* voxel from real
   */
  _MATRIX to_voxel;

  


} FIELD;



/* FIELD management */

extern int BAL_AllocFieldGeometry( FIELD *field );
extern int BAL_SetFieldVoxelSizes( FIELD *field, typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz );

extern int BAL_AllocateField ( FIELD * field, size_t npoints );

extern void BAL_FreeField ( FIELD * field );

/* MISC */

extern void CreateFileDef(FIELD *field,
		   char *nom_image, char *nom_champ );

extern void BAL_PrintField( FILE *f, FIELD *field );
extern void BAL_PrintSelectedPairsOfField( FILE *f, FIELD *field );
extern void BAL_PrintValidPairsOfField( FILE *f, FIELD *field );

extern void BAL_ChangeFieldToVoxelUnit( FIELD *field );
extern void BAL_ChangeFieldToRealUnit( FIELD *field );


#ifdef __cplusplus
}
#endif

#endif
