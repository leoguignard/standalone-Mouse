




/* Resampling procedure.

   Work for 3D images, not for vectorial ones.
   
   (double* mat) is the matrix which permits to get
   from resBuf into theBuf. 
   If one only have the matrix from theBuf into resBuf,
   it must be inverted first.

   Soit x le point transforme et ix=(int)x;
   nous allons distinguer les cas suivants :
    x < -0.5               => resultat = 0
    -0.5 <= x < 0.0        => ix=0, on n'interpole pas selon X
    0.0 < x && ix < dimx-1 => on interpole selon X
    x < dimx-0.5           => ix=dimx-1, on n'interpole pas selon X
    x >= dimx-0.5          => resultat = 0

*/

static void *_Cspline3D4x4_TYPE ( void *par )
{
  char *proc = "_Cspline3D4x4_TYPE";
  typeChunk *chunk = (typeChunk *)par;
  void *parameter = chunk->parameters;
  size_t first = chunk->first;
  size_t last = chunk->last;

  _CsplineResamplingParam *p = (_CsplineResamplingParam *)parameter;
  
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  double* mat = p->mat; 
  int *derivative = p->derivative;
  typeCSplineCoefficients *theCoeff = p->theCoeff;

  size_t i, j, k;
  register int ix, iy, iz;
  register double x, y, z, dx, dy, dz;
  register double res;
  size_t rdimx=resDim[0], rdimy=resDim[1];

  register TYPE *rbuf = (TYPE*)resBuf;

  register float *tbuf;
  float *bufCoeff = theCoeff->theCoeff;

  int tdimx  = theCoeff->theDim[0];
  int tdimxy = theCoeff->theDim[0] * theCoeff->theDim[1];

  int tdimx3 = theCoeff->theDim[0]-3;
  int tdimy3 = theCoeff->theDim[1]-3;
  int tdimz3 = theCoeff->theDim[2]-3;

  register double ddimx = (double)theCoeff->theDim[0]-2.0-0.5;
  register double ddimy = (double)theCoeff->theDim[1]-2.0-0.5;
  register double ddimz = (double)theCoeff->theDim[2]-2.0-0.5;

  double cx[4], cy[4], cz[4], cyz;
  size_t offset[64];

  int l, m, n, o;

  double (*bsplinx)(double);
  double (*bspliny)(double);
  double (*bsplinz)(double);

  size_t ifirst, jfirst, kfirst;
  size_t ilast, jlast, klast;
  size_t iend, jend;



  if ( derivative == (int*) NULL ) {
    bsplinx = & CubicBspline;
    bspliny = & CubicBspline;
    bsplinz = & CubicBspline;
  } else {
    switch ( derivative[0] ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: invalid derivation order in x\n", proc );
      chunk->ret = -1;
      return( (void*)NULL );
    case 0 : bsplinx = & CubicBspline;            break;
    case 1 : bsplinx = & CubicBspline_FirstDeriv; break;
    case 2 : bsplinx = & CubicBspline_SecndDeriv; break;
    }

    switch ( derivative[1] ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: invalid derivation order in y\n", proc );
      chunk->ret = -1;
      return( (void*)NULL );
    case 0 : bspliny = & CubicBspline;            break;
    case 1 : bspliny = & CubicBspline_FirstDeriv; break;
    case 2 : bspliny = & CubicBspline_SecndDeriv; break;
    }

    switch ( derivative[2] ) {
    default :
      if ( _verbose_ )
        fprintf( stderr, "%s: invalid derivation order in z\n", proc );
      chunk->ret = -1;
      return( (void*)NULL );
    case 0 : bsplinz = & CubicBspline;            break;
    case 1 : bsplinz = & CubicBspline_FirstDeriv; break;
    case 2 : bsplinz = & CubicBspline_SecndDeriv; break;
    }
  }



  for ( o = 0, l = 0; l < 4; l ++ )
  for ( m = 0; m < 4; m ++ )
  for ( n = 0; n < 4; n ++, o ++ ) {
    offset[o] = (size_t)l * (size_t)tdimxy + (size_t)m * (size_t)tdimx + (size_t)n;
  }



  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %lu\r", k );
    jend = (k==klast) ? jlast+1 : rdimy;
    for ( ; j<jend; j++, i=0 ) {
      iend = (j==jlast && k==klast) ? ilast+1 : rdimx;
      for ( ; i<iend; i++, rbuf++ ) {
        /* computation of the corresponding point after deformation

           (x,y,x) -> (ix=(int)x,iy=(int)x,iz=(int)x)
           les coefficients de la spline a prendre en compte
           sont [1+ix+{-1,0,1,2}, 1+iy+{-1,0,1,2}, 1+iz+{-1,0,1,2}]

           On ajoute 1 aux coordonnees, car on a ajoute
           une bordure aux coefficients de la spline.
         */

	x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
	if ((x <= -0.5) || ( x >= ddimx)) { *rbuf = 0; continue; }
	y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
	if ((y <= -0.5) || ( y >= ddimy)) { *rbuf = 0; continue; }
	z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
	if ((z <= -0.5) || ( z >= ddimz)) { *rbuf = 0; continue; }
	
        /* ici on a 0 <   x   < theCoeff->theDim[0]-3
                    0 <= ix   < theCoeff->theDim[0]-3
                    1 <= ix+1 < theCoeff->theDim[0]-2
                    0 <= ix+1-1 && ix+1+2 < theCoeff->theDim[0]

           maintenant on a -0.5 <         x   <        theCoeff->theDim[0]-2.5
                              0 <= ix   <= theCoeff->theDim[0]-3
                              1 <= ix+1 <= theCoeff->theDim[0]-2

           les extrema pour les coefficients verifient
                              0 <= ix+1-1 <= theCoeff->theDim[0]-3
                              3 <= ix+1+2 <= theCoeff->theDim[0]
           dans ce dernier cas, on suppose que ca va etre 0 ...
        */
        ix = (int)x;    dx = x-ix;
        iy = (int)y;    dy = y-iy;
        iz = (int)z;    dz = z-iz;

        if ( x < 0.0 ) { x = 0.0; ix = 0; dx = 0.0; }
        if ( y < 0.0 ) { y = 0.0; iy = 0; dy = 0.0; }
        if ( z < 0.0 ) { z = 0.0; iz = 0; dz = 0.0; }

        /* Pre-calcul des valeurs de bsplines pour chaque dimension:

           (*bsplinx)(dx-1), (*bsplinx)(dx), (*bsplinx)(dx+1), (*bsplinx)(dx+2)
           (*bspliny)(dy-1), (*bspliny)(dy), (*bspliny)(dy+1), (*bspliny)(dy+2)
           (*bsplinz)(dz-1), (*bsplinz)(dz), (*bsplinz)(dz+1), (*bsplinz)(dz+2)

           A. Roche

           Dans l'image des coefficients,
           (*bsplinx)(dx-1) correspond a ix,
           (*bsplinx)(dy-1) correspond a iy,
           (*bsplinx)(dz-1) correspond a iz,

        */
        for ( l = 0; l < 4; l ++ ) {
          cx[l] = (*bsplinx)( dx - (double)(l-1) );
          cy[l] = (*bspliny)( dy - (double)(l-1) );
          cz[l] = (*bsplinz)( dz - (double)(l-1) );
        }

        tbuf = bufCoeff;
        tbuf += (size_t)iz * (size_t)tdimxy + (size_t)iy * (size_t)tdimx + (size_t)ix;



        /*
         * are we on the border or not ?
         */
        if ( (ix < tdimx3) && (iy < tdimy3) && (iz < tdimz3) ) {


          for ( res = 0.0, o = 0, l = 0; l < 4; l ++ )
          for ( m = 0; m < 4; m ++ ) {
            cyz = cz[l]*cy[m];
            for ( n = 0; n < 4; n ++, o ++ ) {
              res += cx[n]*cyz*tbuf[ offset[ o ] ];
            }
          }
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }

        /*
         * here, we are sure we are on some border
         * if ix == tdimx3  => ix = theCoeff->theDim[0]-3
         */

        *rbuf = 0;
        *rbuf = (TYPE)_CONVERT_( 0 );

        if ( ix == tdimx3 ) {
          if ( iy == tdimy3 ) {
            if ( iz == tdimz3 ) {

              for ( res = 0.0, o = 0, l = 0; l < 3; l++, o+=4 )
              for ( m = 0; m < 3; m++, o++ ) {
                cyz = cz[l]*cy[m];
                for ( n = 0; n < 3; n ++, o ++ )
                  res += cx[n]*cyz*tbuf[ offset[ o ] ];
              }
              *rbuf = (TYPE)_CONVERT_( res );
              continue;
            }
            for ( res = 0.0, o = 0, l = 0; l < 4; l++, o+=4 )
            for ( m = 0; m < 3; m++, o++ ) {
              cyz = cz[l]*cy[m];
              for ( n = 0; n < 3; n ++, o ++ )
                res += cx[n]*cyz*tbuf[ offset[ o ] ];
            }
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }

          if ( iz == tdimz3 ) {
            for ( res = 0.0, o = 0, l = 0; l < 3; l++ )
            for ( m = 0; m < 4; m++, o++ ) {
              cyz = cz[l]*cy[m];
              for ( n = 0; n < 3; n ++, o ++ )
                res += cx[n]*cyz*tbuf[ offset[ o ] ];
            }
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          for ( res = 0.0, o = 0, l = 0; l < 4; l++ )
          for ( m = 0; m < 4; m++, o++ ) {
            cyz = cz[l]*cy[m];
            for ( n = 0; n < 3; n ++, o ++ )
              res += cx[n]*cyz*tbuf[ offset[ o ] ];
          }
          *rbuf = (TYPE)_CONVERT_( res );
          continue;

        }

        if ( iy == tdimy3 ) {

          if ( iz == tdimz3 ) {
            for ( res = 0.0, o = 0, l = 0; l < 3; l++, o+=4 )
            for ( m = 0; m < 3; m++ ) {
              cyz = cz[l]*cy[m];
              for ( n = 0; n < 4; n ++, o ++ )
                res += cx[n]*cyz*tbuf[ offset[ o ] ];
            }
            *rbuf = (TYPE)_CONVERT_( res );
            continue;
          }
          for ( res = 0.0, o = 0, l = 0; l < 4; l++, o+=4 )
          for ( m = 0; m < 3; m++ ) {
            cyz = cz[l]*cy[m];
            for ( n = 0; n < 4; n ++, o ++ )
              res += cx[n]*cyz*tbuf[ offset[ o ] ];
          }
          *rbuf = (TYPE)_CONVERT_( res );
          continue;

        }

        if ( iz == tdimz3 ) {
          for ( res = 0.0, o = 0, l = 0; l < 3; l++ )
          for ( m = 0; m < 4; m++ ) {
            cyz = cz[l]*cy[m];
            for ( n = 0; n < 4; n ++, o ++ )
              res += cx[n]*cyz*tbuf[ offset[ o ] ];
          }
          *rbuf = (TYPE)_CONVERT_( res );
          continue;
        }
        for ( res = 0.0, o = 0, l = 0; l < 4; l++ )
        for ( m = 0; m < 4; m++ ) {
          cyz = cz[l]*cy[m];
          for ( n = 0; n < 4; n ++, o ++ )
            res += cx[n]*cyz*tbuf[ offset[ o ] ];
        }
        *rbuf = (TYPE)_CONVERT_( res );
        continue;

      }
    }
  }

  chunk->ret = 1;
  return( (void*)NULL );
}





int Cspline3D4x4_TYPE ( typeCSplineCoefficients *theCoeff,
                         void* resBuf, /* result buffer */
                         int *resDim,  /* dimensions of this buffer */
                         double* mat,   /* transformation matrix */
                         int *derivative
                         )
{
  char *proc = "Cspline3D4x4_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _CsplineResamplingParam p;
  
  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  
  p.theCoeff = theCoeff;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.mat = mat;
  p.derivative = derivative;
  p.gain = 1.0;
  p.bias = 0.0;
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Cspline3D4x4_TYPE, &chunks, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to resample image\n", proc );
    freeChunks( &chunks );
    return( -1 );
  }
  
  freeChunks( &chunks );
  return( 1 );
}



