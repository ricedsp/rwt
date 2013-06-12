/*
File Name: MDWT.c
Last Modification Date:	06/14/95	13:15:44
Current Version: MDWT.c	2.4
File Creation Date: Wed Oct 19 10:51:58 1994
Author: Markus Lang  <lang@jazz.rice.edu>

Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
Created by Markus Lang, Department of ECE, Rice University. 

This software is distributed and licensed to you on a non-exclusive 
basis, free-of-charge. Redistribution and use in source and binary forms, 
with or without modification, are permitted provided that the following 
conditions are met:

1. Redistribution of source code must retain the above copyright notice, 
   this list of conditions and the following disclaimer.
2. Redistribution in binary form must reproduce the above copyright notice, 
   this list of conditions and the following disclaimer in the 
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software 
   must display the following acknowledgment: This product includes 
   software developed by Rice University, Houston, Texas and its contributors.
4. Neither the name of the University nor the names of its contributors 
   may be used to endorse or promote products derived from this software 
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS, 
AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY 
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For information on commercial licenses, contact Rice University's Office of 
Technology Transfer at techtran@rice.edu or (713) 348-6173

Change History: Fixed the code such that 1D vectors passed to it can be in
                either passed as a row or column vector. Also took care of 
		the code such that it will compile with both under standard
		C compilers as well as for ANSI C compilers
		Jan Erik Odegard <odegard@ece.rice.edu> Wed Jun 14 1995

%y = mdwt(x,h,L);
% 
% function computes the discrete wavelet transform y for a 1D or 2D input
% signal x.
%
%    Input:
%	x    : finite length 1D or 2D signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: midwt, mrdwt, mirdwt
*/

/*! \file dwt.c
    \brief Implementation of the discrete wavelet transform

    This file actually implements the wavelet transform. 
    It references rwt_platform.h and does not use any matlab-specific functions, 
    though it does use the mat() macro which describes the in-memory layout of 
    1D and 2D matrices in matlab
*/

#include "rwt_platform.h"

/*!
 * Perform convolution for dwt
 *
 * @param x_in input signal values
 * @param lx the length of x
 * @param h0 the low pass coefficients
 * @param h1 the high pass coefficients
 * @param lh_minus_one one less than the number of scaling coefficients
 * @param x_out_low low pass results
 * @param x_out_high high pass results
 * 
 * For the convolution we will calculate the output of the lowpass and highpass filters in parallel
 *
 */
void fpsconv(double *x_in, int lx, double *h0, double *h1, int lh_minus_one, double *x_out_low, double *x_out_high) {
  int i, j, ind;
  double x0, x1;
  for (i=lx; i<lx+lh_minus_one; i++) {
    x_in[i] = *(x_in+(i-lx));
  }
  ind = 0;
  for (i=0; i<(lx); i+=2) {
    x0 = 0;
    x1 = 0;
    for (j=0; j<=lh_minus_one; j++) {
      x0 = x0 + x_in[i+j] * h0[lh_minus_one-j];
      x1 = x1 + x_in[i+j] * h1[lh_minus_one-j];
    }
    x_out_low[ind] = x0;
    x_out_high[ind++] = x1;
  }
}


/*!
 * Allocate memory for dwt
 *
 * @param m the number of rows of the input matrix
 * @param n the number of columns of the input matrix
 * @param lh the number of scaling coefficients
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param h0
 * @param h1
 *
 * The low pass and high pass filter coefficients are the same size as the scaling coefficients
 * For the output storage area we will need as much space as the input: m*n
 * For the input storage area we will need the same plus one less than the length of the coeffiecients
 */
void dwt_allocate(int m, int n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  *xdummy       = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *h0           = (double *) rwt_calloc(lh,            sizeof(double));
  *h1           = (double *) rwt_calloc(lh,            sizeof(double));
}


/*!
 * Free memory that we allocated for dwt
 *
 * @param xdummy
 * @param y_dummy_low
 * @param y_dummy_high
 * @param h0
 * @param h1
 *
 */
void dwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*h0);
  rwt_free(*h1);
}


/*!
 * Put the scaling coeffients into a form ready for use in the convolution function
 *
 * @param lh length of h / the number of scaling coefficients
 * @param h  the wavelet scaling coefficients
 * @param h0 the high pass coefficients - reversed h
 * @param h1 the high pass coefficients - forward h, even values are sign flipped
 *
 * The coefficients of our Quadrature Mirror Filter are described by
 * \f$ g\left[lh - 1 - n \right] = (-1)^n * h\left[n\right] \f$
 *
 */
void dwt_coefficients(int lh, double *h, double **h0, double **h1) {
  int i;
  for (i=0; i<lh; i++) {
    (*h0)[i] = h[(lh-i)-1];
    (*h1)[i] = h[i];
  }
  for (i=0; i<lh; i+=2)
    (*h1)[i] = -((*h1)[i]);
}


/*!
 * Perform the discrete wavelet transform
 *
 * @param x  the input signal
 * @param m  number of rows in the input
 * @param n  number of columns in the input
 * @param h  wavelet scaling coefficients
 * @param lh length of h / the number of scaling coefficients
 * @param L  the number of levels
 * @param y  the output signal with the wavelet transform applied
 *
 * The discrete wavelet transform begins with a set of samples of a signal whose length
 * is a power of 2. This exponent we shall call 'L' as it corresponds to the number of
 * levels in the filter bank for the calculation of the wavelet transform. 
 *
 * We shall use the name 'a' to refer to the approximation coefficients.
 * However, the actual implementation will not store this information separately, 
 * but rather will place it in the available space in the output array.
 *
 */
void dwt(double *x, int m, int n, double *h, int lh, int L, double *y) {
  double  *h0, *h1, *y_dummy_low, *y_dummy_high, *xdummy;
  long i;
  int actual_L, lh_minus_one;
  int actual_m, actual_n, row_of_a, column_of_a, idx_rows, idx_columns;

  if (n==1) { /* If passed a 1d column vector, just treat it as a row vector */
    n = m;
    m = 1;
  }
  
  dwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
  dwt_coefficients(lh, h, &h0, &h1); /*! For performance, precalculate what we can outside the loops */
  lh_minus_one = lh - 1;
  actual_m = 2*m;
  actual_n = 2*n;
 
  for (actual_L=1; actual_L<=L; actual_L++) {
    if (m==1)
      actual_m = 1;
    else{
      actual_m = actual_m/2;
      row_of_a = actual_m/2;     
    }
    actual_n = actual_n/2;
    column_of_a = actual_n/2;

    for (idx_rows=0; idx_rows<actual_m; idx_rows++) {
      for (i=0; i<actual_n; i++)
	if (actual_L==1)  
	  xdummy[i] = mat(x, idx_rows, i, m);  
	else 
	  xdummy[i] = mat(y, idx_rows, i, m);  
      /* perform filtering lowpass and highpass*/
      fpsconv(xdummy, actual_n, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      idx_columns = column_of_a;
      for (i=0; i<column_of_a; i++) {    
	mat(y, idx_rows, i,             m) = y_dummy_low[i];  
	mat(y, idx_rows, idx_columns++, m) = y_dummy_high[i];  
      } 
    }  
    
    /*! For the 2d transform, we go through each of the columns after having gone through the rows */
    if (m>1) {
      for (idx_columns=0; idx_columns<actual_n; idx_columns++) { /* loop over columns */
	/*! Store in dummy variables */
	for (i=0; i<actual_m; i++)
	  xdummy[i] = mat(y, i, idx_columns, m);  
	/*! Perform filtering lowpass and highpass*/
	fpsconv(xdummy, actual_m, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
	/*! Restore dummy variables in matrix */
	idx_rows = row_of_a;
	for (i=0; i<row_of_a; i++) {
	  mat(y, i, idx_columns, m) = y_dummy_low[i];  
	  mat(y, idx_rows++, idx_columns, m) = y_dummy_high[i];  
	}
      }
    }
  }
  dwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
}

