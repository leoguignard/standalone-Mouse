


extern void Reech3DTriLinVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef  /* transformation matrix */
			     );

extern void Reech3DTriLinVectorFieldgb_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef,  /* transformation matrix */
			     float gain,
			     float bias );

extern void Reech3DNearestVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef  /* transformation matrix */
			     );

extern void Reech2DTriLinVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef  /* transformation matrix */
			     );

extern void Reech2DTriLinVectorFieldgb_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef,  /* transformation matrix */
			     float gain,
			     float bias );

extern void Reech2DNearestVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef  /* transformation matrix */
			     );
