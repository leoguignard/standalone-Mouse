/*************************************************************************
 * bal-image.h -
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




#ifndef BAL_IMAGE_H
#define BAL_IMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-stddef.h>
#include <bal-behavior.h>
#include <bal-matrix.h>

#include <typedefs.h>
#include <linearFiltering-common.h>




typedef struct bal_image {
        size_t ncols;         /* Number of columns (X dimension) */
        size_t nrows;         /* Number of rows (Y dimension) */
        size_t nplanes;       /* Number of planes (Z dimension) */
        size_t vdim;          /* Vector size */
        bufferType type;
        void *data;        /* Generic pointer on image data buffer.
                              This pointer has to be casted in proper data type
                              depending on the type field */
        void ***array;     /* Generic 3D array pointing on each image element.
                              This pointer has to be casted in proper data type
                              depending on the type field */

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

        /* QForm and SForm codes
         */
        int qform_code;
        int sform_code;
        _MATRIX sform_to_real;

        char *name;
} bal_image;



extern int BAL_GetVerboseInBalImage(  );
extern void BAL_SetVerboseInBalImage( int v );
extern void BAL_IncrementVerboseInBalImage(  );
extern void BAL_DecrementVerboseInBalImage(  );

extern void BAL_SetDebugInBalImage( int d );
extern void BAL_IncrementDebugInBalImage(  );
extern void BAL_DecrementDebugInBalImage(  );





/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image geometry
 *
 --------------------------------------------------*/

extern int BAL_AllocImageGeometry( bal_image *image );

extern int BAL_SetImageVoxelSizes( bal_image *image,
                                   typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz );

extern int BAL_CopyImageGeometry( bal_image *theIm, bal_image *resIm );

extern int BAL_ResizeImageGeometry( bal_image *theIm, bal_image *resIm );





/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image structure
 *
 --------------------------------------------------*/

/* 1. fill the bal_image structure with dimensions
 * 2. does not set the voxel to/from real geometry
 * 3. does not allocate the buffers
 *
 * use BAL_AllocImage() to allocate the image
 * use BAL_SetImageVoxelSizes to set the voxel sizes
 *
 * do not require desallocation
 */
extern int  BAL_InitImage( bal_image *image, char *name,
                                     int dimx, int dimy, int dimz, int dimv,
                                     bufferType type );

/* 1. fill the bal_image structure with dimensions
 * 2. set the voxel to/from real geometry
 * 3. does not allocate the buffers
 *
 * use BAL_AllocImage() to allocate the image
 *
 * require desallocation with BAL_FreeImage()
 */
extern int BAL_InitFullImage( bal_image *image, char *name,
                          int dimx, int dimy, int dimz, int dimv,
                          typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz,
                          bufferType type );


/* does nothing but a call to BAL_InitImageGeometry()
 * parameters are issued from image 'from'
 */
extern int BAL_InitImageFromImage( bal_image *image, char *name,
                                 bal_image *from, bufferType type );

extern int BAL_InitScalarImageFromImage( bal_image *image, char *name,
                                  bal_image *from, bufferType type );




/*--------------------------------------------------
 *
 * IMAGE MANAGEMENT: image allocation
 *
 --------------------------------------------------*/

extern void BAL_FreeImage( bal_image *image );

/* allocate and fill the array pointer
 * to be used when the data buffer already exists
 */
extern int BAL_AllocArrayImage( bal_image *image );

/* Allocate the buffers of the bal_image structure.
 * There are two buffers:
 * - data is a 1D buffer containing the data
 * - array is a multiple pointer, that points to convenient places
 *   of 'data' so ((cast)array)[z][y][x] gives access to the voxel value
 * values are initialized at 0.
 */
extern int  BAL_AllocImage( bal_image *image );



/* does nothing but a call to BAL_InitImageGeometry()
 * followed by a call to BAL_AllocImage()
 *
 */
extern int BAL_AllocFullImage( bal_image *image, char *name,
                               int dimx, int dimy, int dimz, int dimv,
                               typeVoxelSize vx, typeVoxelSize vy, typeVoxelSize vz,
                               bufferType type );



/* does nothing but a call to BAL_InitImageFromImage()
 * followed by a call to BAL_AllocImage()
 */
extern int BAL_AllocImageFromImage( bal_image *image, char *name,
                                    bal_image *from, bufferType type );

/* same as above but set 'vdim' to 1
 * to be used for allocation from a vectorial image
 */
extern int BAL_AllocScalarImageFromImage( bal_image *image, char *name,
                                 bal_image *from, bufferType type );





extern int BAL_CopyImage( bal_image *theIm, bal_image *resIm );




/* allows to fill the buffer image with one given value
 */
extern int BAL_FillImage( bal_image *theIm, float fval );
extern int BAL_NormaImage( bal_image *theIm, bal_image *resIm );
extern size_t  BAL_ImageDataSize( bal_image *image );



/* Get (interpolate) image value at the given position
   (in voxels not in millimeters)
*/
extern double BAL_GetXYZvalue( bal_image *image, double x, double y, double z );
extern double BAL_GetXYKvalue( bal_image *image, double x, double y, int k );



/* I/O operation
 */
extern void BAL_PrintImage( FILE *f, bal_image *image, char *s );
extern void BAL_PrintParImage( FILE *f, bal_image *image, char *s );

extern int BAL_ReadImage ( bal_image *image, char *name, int normalisation );
extern int BAL_WriteImage( bal_image *image, char *name );


/* filtering
 */
extern void BAL_SetFilterType( filterType filter );
extern int BAL_SmoothImage( bal_image *theIm,
                            bal_doublePoint *theSigma );
extern int BAL_SmoothImageIntoImage( bal_image *theIm, bal_image *resIm,
                                     bal_doublePoint *theSigma );
extern int BAL_2DDerivativesOfImage( bal_image *theIm, 
                                     bal_image *theDx, bal_image *theDy,
                                     bal_doublePoint *theSigma );
extern int BAL_3DDerivativesOfImage( bal_image *theIm, 
                                     bal_image *theDx, bal_image *theDy, bal_image *theDz,
                                     bal_doublePoint *theSigma );

#ifdef __cplusplus
}
#endif

#endif
