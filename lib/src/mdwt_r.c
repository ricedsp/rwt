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

#include "platform.h"


void fpsconv(double *x_in, int lx, double *h0, double *h1, int lh_minus_one, double *x_outl, double *x_outh) {
  int i, j, ind;
  double x0, x1;
  for (i=lx; i<lx+lh_minus_one; i++)
    x_in[i] = *(x_in+(i-lx));
  ind = 0;
  for (i=0; i<(lx); i+=2) {
    x0 = 0;
    x1 = 0;
    for (j=0; j<=lh_minus_one; j++) {
      x0 = x0 + x_in[i+j] * h0[lh_minus_one-j];
      x1 = x1 + x_in[i+j] * h1[lh_minus_one-j];
    }
    x_outl[ind] = x0;
    x_outh[ind++] = x1;
  }
}


void mdwt_allocate(int m, int n, int lh, double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  *xdummy       = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low  = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_high = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *h0           = (double *) rwt_calloc(lh,            sizeof(double));
  *h1           = (double *) rwt_calloc(lh,            sizeof(double));
}


void mdwt_free(double **xdummy, double **y_dummy_low, double **y_dummy_high, double **h0, double **h1) {
  rwt_free(*xdummy);
  rwt_free(*y_dummy_low);
  rwt_free(*y_dummy_high);
  rwt_free(*h0);
  rwt_free(*h1);
}


/* h0 <- reversed h
   h1 <- forward h, even values are sign reversed */
void mdwt_coefficients(int lh, double *h, double **h0, double **h1) {
  int i;
  for (i=0; i<lh; i++) {
    (*h0)[i] = h[(lh-i)-1];
    (*h1)[i] = h[i];
  }
  for (i=0; i<lh; i+=2)
    (*h1)[i] = -((*h1)[i]);
}


void MDWT(double *x, int m, int n, double *h, int lh, int L, double *y) {
  double  *h0, *h1, *y_dummy_low, *y_dummy_high, *xdummy;
  int i, actual_L, lh_minus_one;
  int upsampled_rows, upsampled_columns, pass_rows, pass_columns, idx_rows, idx_columns;
  
  mdwt_allocate(m, n, lh, &xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
  mdwt_coefficients(lh, h, &h0, &h1);

  /* analysis lowpass and highpass */
  if (n==1) {
    n = m;
    m = 1;
  }
  
  lh_minus_one = lh - 1;
  upsampled_rows = 2*m;
  upsampled_columns = 2*n;
 
  //mexPrintf("new signal. n is %d\n", n);
 
  /* main loop */
  for (actual_L=1; actual_L<=L; actual_L++) {
    if (m==1)
      upsampled_rows = 1;
    else{
      upsampled_rows = upsampled_rows/2;
      pass_rows = upsampled_rows/2;     
    }
    upsampled_columns = upsampled_columns/2;
    pass_columns = upsampled_columns/2;

    //mexPrintf("1d: %d %d\n", upsampled_columns, pass_columns);
    /* go by rows */
    for (idx_rows=0; idx_rows<upsampled_rows; idx_rows++) {            /* loop over rows */
      /* store in dummy variable */
// Why do we copy in and out of dummy vars?
      for (i=0; i<upsampled_columns; i++)
	if (actual_L==1)  
	  xdummy[i] = mat(x, idx_rows, i, m);  
	else 
	  xdummy[i] = mat(y, idx_rows, i, m);  
      /* perform filtering lowpass and highpass*/
      fpsconv(xdummy, upsampled_columns, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
      /* restore dummy variables in matrices */
      idx_columns = pass_columns;
      for  (i=0; i<pass_columns; i++) {    
	mat(y, idx_rows, i, m) = y_dummy_low[i];  
	mat(y, idx_rows, idx_columns++, m) = y_dummy_high[i];  
      } 
    }  
    
    /* go by columns in case of a 2D signal*/
    if (m>1) {
      //mexPrintf("2d: %d %d %d %d\n", upsampled_rows, pass_rows, upsampled_columns, pass_columns);
      for (idx_columns=0; idx_columns<upsampled_columns; idx_columns++) { /* loop over columns */
	/* store in dummy variables */
	for (i=0; i<upsampled_rows; i++)
	  xdummy[i] = mat(y, i, idx_columns, m);  
	/* perform filtering lowpass and highpass*/
	fpsconv(xdummy, upsampled_rows, h0, h1, lh_minus_one, y_dummy_low, y_dummy_high); 
	/* restore dummy variables in matrix */
	idx_rows = pass_rows;
	for (i=0; i<pass_rows; i++) {
	  mat(y, i, idx_columns, m) = y_dummy_low[i];  
	  mat(y, idx_rows++, idx_columns, m) = y_dummy_high[i];  
	}
      }
    }
  }
  mdwt_free(&xdummy, &y_dummy_low, &y_dummy_high, &h0, &h1);
}

