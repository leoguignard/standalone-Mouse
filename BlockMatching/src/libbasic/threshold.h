/*************************************************************************
 * threshold.h - image thresholding
 *
 * $$
 *
 * Copyright (c) INRIA 2016, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 *
 * CREATION DATE:
 * Jeu 17 nov 2016 09:17:49 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#ifndef _threshold_h_
#define _threshold_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>

extern int thresholdBuffer( void *bufferIn,
                            bufferType typeIn,
                            void *bufferOut,
                            bufferType typeOut,
                            int *bufferDims,
                            float threshold );

#ifdef __cplusplus
}
#endif

#endif
