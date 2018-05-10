/*************************************************************************
 * bal-interpolation.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun  4 avr 2016 22:07:46 CEST
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_INTERPOLATION_H
#define BAL_INTERPOLATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include <bal-stddef.h>
#include <bal-transformation.h>
#include <bal-transformation-tools.h>

extern void BAL_SetVerboseInBalInterpolation( int v );
extern void BAL_IncrementVerboseInBalInterpolation(  );
extern void BAL_DecrementVerboseInBalInterpolation(  );

extern void BAL_SetDebugInBalInterpolation( int d );
extern void BAL_IncrementDebugInBalInterpolation(  );
extern void BAL_DecrementDebugInBalInterpolation(  );




/* a list of labels
 * used for correspondences
 * (list of labels) - (list of labels)
 */
typedef struct typeLabelList {
  int *data;
  int n_data;
  int n_allocated_data;
} typeLabelList;



typedef struct typeCell {
  bal_integerPoint left_corner;
  bal_integerPoint right_corner;
  int npoints;
  int new_label;
  typeLabelList corresp;
} typeCell;

typedef struct typeCellList {
    typeCell *data;
    int n_data;
    int n_allocated_data;
} typeCellList;

extern int BAL_FillCellList( bal_image *image, typeCellList *l );



typedef struct typeCorrespondence {
  typeLabelList *label0;
  typeLabelList *label1;
  typeLabelList allocatedLabel0;
  typeLabelList allocatedLabel1;
} typeCorrespondence;

typedef struct typeCorrespondenceList {
  typeCorrespondence *data;
  int n_data;
  int n_allocated_data;
} typeCorrespondenceList;

extern int BAL_ReadCorrespondenceList( typeCorrespondenceList *l, char *name );
extern int BAL_WriteCorrespondenceList( typeCorrespondenceList *l, char *name );



typedef struct typeCellCorrespondence {
  typeCorrespondenceList cliques;
  typeCellList *list0;
  typeCellList *list1;
  typeCellList allocatedList0;
  typeCellList allocatedList1;
} typeCellCorrespondence;

extern void BAL_InitCellCorrespondence( typeCellCorrespondence *c );
extern void BAL_FreeCellCorrespondence( typeCellCorrespondence *c );

extern int BAL_CrossFillCellCorrespondence( typeCellCorrespondence *cc );
extern void BAL_CheckCellCorrespondence( FILE *f, typeCellCorrespondence *cc );

extern int BAL_RelabelCells( bal_image *image0, bal_image *image1,
                             typeCellCorrespondence *cc );






extern int BAL_InterpolateGreyLevelImages( bal_image *theImage0, bal_image *theImage1,
                                    bal_image *resImage0, bal_image *resImage1,
                                    bal_image *resImage,
                                    bal_transformation *theTrsf,
                                    float t,
                                    enumTransformationInterpolation interpolation );

extern int BAL_InterpolateLabelImages( bal_image *theImage0, bal_image *theImage1,
                                bal_image *resImage0, bal_image *resImage1,
                                bal_transformation *theTrsf,
                                typeCellCorrespondence *cc,
                                float t,
                                float rmax );

#ifdef __cplusplus
}
#endif

#endif
