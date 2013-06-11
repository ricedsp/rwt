/*
File Name: MIRDWT.c
Last Modification Date:	06/14/95	16:22:45
Current Version: MIRDWT.c	2.4
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

MATLAB description:
%function x = mirdwt(y_low,y_high,h,L);
% 
% function computes the inverse redundant discrete wavelet transform y for a
% 1D or  2D input signal. redundant means here that the subsampling after
% each stage of the forward transform has been omitted. y_low contains the
% lowpass and y_low the highpass components as computed, e.g., by mrdwt. In
% case of a 2D signal the ordering in y_high is [lh hl hh lh hl ... ] (first
% letter refers to row, second to column filtering).  
%
%    Input:
%       y_low   : lowpass component
%       y_high   : highpass components
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(y_low) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%   
%    Output:
%	x    : finite length 1D or 2D signal
%
% see also: mdwt, midwt, mrdwt

*/

/*! \file irdwt.c
    \brief Implementation of the inverse redundant discrete wavelet transform

*/

#include "rwt_platform.h"

void bpconv(double *x_out, int lx, double *g0, double *g1, int lh, double *x_inl, double *x_inh) {
  int i, j;
  double x0;

  for (i=lh-2; i > -1; i--){
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }
  for (i=0; i<lx; i++){
    x0 = 0;
    for (j=0; j<lh; j++)
      x0 = x0 + x_inl[j+i]*g0[lh-1-j] +
	x_inh[j+i]*g1[lh-1-j];
    x_out[i] = x0;
  }
}


void irdwt_allocate(int m, int n, int lh, double **x_high, double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, 
  double **y_dummy_low_high, double **y_dummy_high_low, double **y_dummy_high_high, double **g0, double **g1) {
  *x_high            = (double *) rwt_calloc(m*n,           sizeof(double));
  *x_dummy_low       = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *x_dummy_high      = (double *) rwt_calloc(max(m,n),      sizeof(double));
  *y_dummy_low_low   = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_low_high  = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_high_low  = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *y_dummy_high_high = (double *) rwt_calloc(max(m,n)+lh-1, sizeof(double));
  *g0                = (double *) rwt_calloc(lh,            sizeof(double));
  *g1                = (double *) rwt_calloc(lh,            sizeof(double));
}


void irdwt_free(double **x_dummy_low, double **x_dummy_high, double **y_dummy_low_low, double **y_dummy_low_high, 
  double **y_dummy_high_low, double **y_dummy_high_high, double **g0, double **g1) {
  rwt_free(*x_dummy_low);
  rwt_free(*x_dummy_high);
  rwt_free(*y_dummy_low_low);
  rwt_free(*y_dummy_low_high);
  rwt_free(*y_dummy_high_low);
  rwt_free(*y_dummy_high_high);
  rwt_free(*g0);
  rwt_free(*g1);
}

/* not the same as idwt_coefficients */
void irdwt_coefficients(int lh, double *h, double **g0, double **g1) {
  int i;
  for (i=0; i<lh; i++){
    (*g0)[i] = h[i]/2;
    (*g1)[i] = h[lh-i-1]/2;
  }
  for (i=1; i<=lh; i+=2)
    (*g1)[i] = -((*g1)[i]);
}


void irdwt(double *x, int m, int n, double *h, int lh, int L, double *y_low, double *y_high) {
  double  *g0, *g1, *y_dummy_low_low, *y_dummy_low_high, *y_dummy_high_low;
  double *y_dummy_high_high, *x_dummy_low , *x_dummy_high, *x_high;
  long i;
  int actual_L, actual_m, actual_n, c_o_a, ir, n_c, n_cb, lh_minus_one;
  int ic, n_r, n_rb, c_o_a_p2n, sample_f;

  irdwt_allocate(m, n, lh, &x_high, &x_dummy_low, &x_dummy_high, &y_dummy_low_low, 
    &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
  irdwt_coefficients(lh, h, &g0, &g1);
  
  if (n==1){
    n = m;
    m = 1;
  }
  /* analysis lowpass and highpass */
  
  lh_minus_one = lh - 1;
  /* 2^L */
  sample_f = 1;
  for (i=1; i<L; i++)
    sample_f = sample_f*2;
  actual_m = m/sample_f;
  actual_n = n/sample_f;
  /* restore y_low in x */
  for (i=0;i<m*n;i++)
    x[i] = y_low[i];
  
  /* main loop */
  for (actual_L=L; actual_L >= 1; actual_L--){
    /* actual (level dependent) column offset */
    if (m==1)
      c_o_a = n*(actual_L-1);
    else
      c_o_a = 3*n*(actual_L-1);
    c_o_a_p2n = c_o_a + 2*n;
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      n_rb = m/actual_m;                 /* # of row blocks per column */
      for (ic=0; ic<n; ic++){            /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++){    /* loop within one column */
	  /* store in dummy variables */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    y_dummy_low_low[i+lh_minus_one] = mat(x, ir, ic, m);  
	    y_dummy_low_high[i+lh_minus_one] = mat(y_high, ir, c_o_a+ic, m);  
	    y_dummy_high_low[i+lh_minus_one] = mat(y_high, ir,c_o_a+n+ic, m);  
	    y_dummy_high_high[i+lh_minus_one] = mat(y_high, ir, c_o_a_p2n+ic, m);   
	  }
	  /* perform filtering and adding: first LL/LH, then HL/HH */
	  bpconv(x_dummy_low, actual_m, g0, g1, lh, y_dummy_low_low, y_dummy_low_high); 
	  bpconv(x_dummy_high, actual_m, g0, g1, lh, y_dummy_high_low, y_dummy_high_high); 
	  /* store dummy variables in matrices */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    mat(x, ir, ic, m) = x_dummy_low[i];  
	    mat(x_high, ir, ic, m) = x_dummy_high[i];  
	  }
	}
      }
    }
    
    /* go by rows */
    n_cb = n/actual_n;                 /* # of column blocks per row */
    for (ir=0; ir<m; ir++){            /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++){    /* loop within one row */      
	/* store in dummy variable */
	ic = -sample_f + n_c;
	for  (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  y_dummy_low_low[i+lh_minus_one] = mat(x, ir, ic, m);  
	  if (m>1)
	    y_dummy_high_high[i+lh_minus_one] = mat(x_high, ir, ic, m);  
	  else
	    y_dummy_high_high[i+lh_minus_one] = mat(y_high, ir, c_o_a+ic, m);  
	} 
	/* perform filtering lowpass/highpass */
	bpconv(x_dummy_low, actual_n, g0, g1, lh, y_dummy_low_low, y_dummy_high_high); 
	/* restore dummy variables in matrices */
	ic = -sample_f + n_c;
	for (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  mat(x, ir, ic, m) = x_dummy_low[i];  
	}
      }
    }
    sample_f = sample_f/2;
    actual_m = actual_m*2;
    actual_n = actual_n*2;
  }
  irdwt_free(&x_dummy_low, &x_dummy_high, &y_dummy_low_low, &y_dummy_low_high, &y_dummy_high_low, &y_dummy_high_high, &g0, &g1);
}

