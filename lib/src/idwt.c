/*
File Name: MIDWT.c
Last Modification Date:	06/14/95	13:01:15
Current Version: MIDWT.c	2.4
File Creation Date: Wed Oct 12 08:44:43 1994
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

		Fix minor bug to allow maximum number of levels

decription of the matlab call:
%y = midwt(x,h,L);
% 
% function computes the inverse discrete wavelet transform y for a 1D or 2D
% input signal x.
%
%    Input:
%	x    : finite length 1D or 2D input signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%
% see also: mdwt, mrdwt, mirdwt

*/

/*! \file idwt.c
    \brief Implementation of the inverse discrete wavelet transform

*/

#include "rwt_platform.h"

void bpsconv(double *x_out, int lx, double *g0, double *g1, int lh_minus_one, int lh_halved_minus_one, double *x_inl, double *x_inh) {
  int i, j, ind, tj;
  double x0, x1;

  for (i=lh_halved_minus_one-1; i > -1; i--){
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }
  ind = 0;
  for (i=0; i<(lx); i++){
    x0 = 0;
    x1 = 0;
    tj = -2;
    for (j=0; j<=lh_halved_minus_one; j++){
      tj+=2;
      x0 = x0 + x_inl[i+j]*g0[lh_minus_one-1-tj] + x_inh[i+j]*g1[lh_minus_one-1-tj] ;
      x1 = x1 + x_inl[i+j]*g0[lh_minus_one-tj] + x_inh[i+j]*g1[lh_minus_one-tj] ;
    }
    x_out[ind++] = x0;
    x_out[ind++] = x1;
  }
}


void idwt_allocate(int m, int n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  *xdummy       = (double *) rwt_calloc(max(m,n),        sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n)+lh/2-1, sizeof(double));
  *g0           = (double *) rwt_calloc(lh,              sizeof(double));
  *g1           = (double *) rwt_calloc(lh,              sizeof(double));
}


void idwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **g0, double **g1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*g0);
  rwt_free(*g1);
}


void idwt_coefficients(int lh, double *h, double **g0, double **g1) {
  int i;
  for (i=0; i<lh; i++){
    (*g0)[i] = h[i];
    (*g1)[i] = h[lh-i-1];
  }
  for (i=1; i<=lh; i+=2)
    (*g1)[i] = -((*g1)[i]);
}


void idwt(double *x, int m, int n, double *h, int lh, int L, double *y) {
  double  *g0, *g1, *y_dummy_low, *y_dummy_high, *xdummy;
  long i;
  int actual_L, actual_m, actual_n, row_of_a, column_of_a, ir, ic, lh_minus_one, lh_halved_minus_one, sample_f;

  idwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
  idwt_coefficients(lh, h, &g0, &g1);

  if (n==1){
    n = m;
    m = 1;
  }
  
  lh_minus_one = lh - 1;
  lh_halved_minus_one = lh/2 - 1;
  /* 2^L */
  sample_f = 1;
  for (i=1; i<L; i++)
    sample_f = sample_f*2;
  
  if (m>1)
    actual_m = m/sample_f;
  else 
    actual_m = 1;
  actual_n = n/sample_f;

  for (i=0; i<(m*n); i++)
    x[i] = y[i];
  
  /* main loop */
  for (actual_L=L; actual_L >= 1; actual_L--){
    row_of_a = actual_m/2;
    column_of_a = actual_n/2;
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      for (ic=0; ic<actual_n; ic++){            /* loop over column */
	/* store in dummy variables */
	ir = row_of_a;
	for (i=0; i<row_of_a; i++){    
	  y_dummy_low[i+lh_halved_minus_one] = mat(x, i, ic, m);  
	  y_dummy_high[i+lh_halved_minus_one] = mat(x, ir++, ic, m);  
	}
	/* perform filtering lowpass and highpass*/
	bpsconv(xdummy, row_of_a, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
	/* restore dummy variables in matrix */
	for (i=0; i<actual_m; i++)
	  mat(x, i, ic, m) = xdummy[i];  
      }
    }
    /* go by rows */
    for (ir=0; ir<actual_m; ir++){            /* loop over rows */
      /* store in dummy variable */
      ic = column_of_a;
      for  (i=0; i<column_of_a; i++){    
	y_dummy_low[i+lh_halved_minus_one] = mat(x, ir, i, m);  
	y_dummy_high[i+lh_halved_minus_one] = mat(x, ir, ic++, m);  
      } 
      /* perform filtering lowpass and highpass*/
      bpsconv(xdummy, column_of_a, g0, g1, lh_minus_one, lh_halved_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      for (i=0; i<actual_n; i++)
        mat(x, ir, i, m) = xdummy[i];  
    }  
    if (m==1)
      actual_m = 1;
    else
      actual_m = actual_m*2;
    actual_n = actual_n*2;
  }
  idwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &g0, &g1);
}

