


extern void Reech3DTriLin4x4_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech3DTriLin4x4gb_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech3DNearest4x4_TYPE ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
extern void Reech2DTriLin4x4_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat   /* transformation matrix */
			     );
extern void Reech2DTriLin4x4gb_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim, /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     double* mat,   /* transformation matrix */
			     float gain,
			     float bias );
extern void Reech2DNearest4x4_TYPE ( void* theBuf, /* buffer to be resampled */
			      int *theDim,  /* dimensions of this buffer */
			      void* resBuf, /* result buffer */
			      int *resDim,  /* dimensions of this buffer */
			      double* mat   /* transformation matrix */
			      );
