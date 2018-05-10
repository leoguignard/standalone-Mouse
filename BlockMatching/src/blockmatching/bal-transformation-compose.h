/*************************************************************************
 * bal-transformation-compose.h -
 *
 * $Id$
 *
 * Composeright (c) INRIA 2017, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Dim 22 jan 2017 18:53:01 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_COMPOSE_H
#define BAL_TRANSFORMATION_COMPOSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-transformation.h>



extern void BAL_SetVerboseInBalTransformationCompose( int v );\
extern void BAL_IncrementVerboseInBalTransformationCompose( );
extern void BAL_DecrementVerboseInBalTransformationCompose( );
extern void BAL_SetDebugInBalTransformationCompose( int v );
extern void BAL_IncrementDebugInBalTransformationCompose( );
extern void BAL_DecrementDebugInBalTransformationCompose( );


/***************************************************
 *
 * transformation composition
 * t_res = t1 o t2
   t_{I3<-I1} = t_{I3<-I2} o t_{I2<-I1}
 *
 ***************************************************/

extern int BAL_TransformationComposition( bal_transformation *t_res,
                                          bal_transformation *t1,
                                          bal_transformation *t2 );

extern int BAL_TransformationListComposition( bal_transformation *t_res,
                                             bal_transformation **trsfsStructure,
                                             int n );

extern enumTypeTransfo BAL_TypeTransformationComposition( bal_transformation *t1,
                                                          bal_transformation *t2 );

extern int BAL_AllocTransformationComposition( bal_transformation *res,
                                               bal_transformation *t1,
                                               bal_transformation *t2,
                                               bal_image *ref );

extern enumTypeTransfo BAL_TypeTransformationListComposition( bal_transformation **array,
                                                              int n );

extern int BAL_AllocTransformationListComposition( bal_transformation *res,
                                                   bal_transformation **array,
                                                   int n,
                                                   bal_image *ref );


#ifdef __cplusplus
}
#endif

#endif
