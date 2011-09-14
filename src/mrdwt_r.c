/*
File Name: MRDWT.c
Last Modification Date:	09/21/95	15:42:59
Current Version: MRDWT.c	2.4
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

MATLAB description:
%[yl,yh] = mrdwt(x,h,L);
% 
% function computes the redundant discrete wavelet transform y for a 1D or
% 2D input signal . redundant means here that the subsampling after each
% stage is omitted. yl contains the lowpass and yl the highpass
% components. In case of a 2D signal the ordering in yh is [lh hl hh lh hl
% ... ] (first letter refers to row, second to column filtering). 
%
%    Input:
%	x    : finite length 1D or 2D signal (implicitely periodized)
%       h    : scaling filter
%       L    : number of levels. in case of a 1D signal length(x) must be
%              divisible by 2^L; in case of a 2D signal the row and the
%              column dimension must be divisible by 2^L.
%   
%    Output:
%       yl   : lowpass component
%       yh   : highpass components
%
% see also: mdwt, midwt, mirdwt


*/

#include <math.h>
#include <stdio.h>

/*#define mat(a, i, j) (a[m*(j)+i]) */
#define mat(a, i, j) (*(a + (m*(j)+i))) 
#define max(a, b) ((a) > (b) ? (a) : (b))

#ifdef __STDC__
MRDWT(double *x, int m, int n, double *h, int lh, int L,
      double *yl, double *yh)
#else
MRDWT(x, m, n, h, lh, L, yl, yh)
double *x, *h, *yl, *yh;
int m, n, lh, L;
#endif
{
  double *tmp;
  double  *h0, *h1, *ydummyll, *ydummylh, *ydummyhl;
  double *ydummyhh, *xdummyl , *xdummyh;
  long i, j;
  int actual_L, actual_m, actual_n, c_o_a, ir, n_c, n_cb, n_c_o;
  int ic, n_r, n_rb, n_r_o, c_o_a_p2n, sample_f;
  xdummyl = (double *)mxCalloc(max(m,n)+lh-1,sizeof(double));
  xdummyh = (double *)mxCalloc(max(m,n)+lh-1,sizeof(double));
  ydummyll = (double *)mxCalloc(max(m,n),sizeof(double));
  ydummylh = (double *)mxCalloc(max(m,n),sizeof(double));
  ydummyhl = (double *)mxCalloc(max(m,n),sizeof(double));
  ydummyhh = (double *)mxCalloc(max(m,n),sizeof(double));
  h0 = (double *)mxCalloc(lh,sizeof(double));
  h1 = (double *)mxCalloc(lh,sizeof(double));

  if (n==1){
    n = m;
    m = 1;
  }  
  /* analysis lowpass and highpass */
  for (i=0; i<lh; i++){
    h0[i] = h[lh-i-1];
    h1[i] =h[i];
  }
  for (i=0; i<lh; i+=2)
    h1[i] = -h1[i];
  
  actual_m = 2*m;
  actual_n = 2*n;
  for (i=0; i<m*n; i++)
    yl[i] = x[i];
  
  /* main loop */
  sample_f = 1;
  for (actual_L=1; actual_L <= L; actual_L++){
    actual_m = actual_m/2;
    actual_n = actual_n/2;
    /* actual (level dependent) column offset */
    if (m==1)
      c_o_a = n*(actual_L-1);
    else
      c_o_a = 3*n*(actual_L-1);
    c_o_a_p2n = c_o_a + 2*n;
    
    /* go by rows */
    n_cb = n/actual_n;                 /* # of column blocks per row */
    for (ir=0; ir<m; ir++){            /* loop over rows */
      for (n_c=0; n_c<n_cb; n_c++){    /* loop within one row */      
	/* store in dummy variable */
	ic = -sample_f + n_c;
	for (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  xdummyl[i] = mat(yl, ir, ic);  
	}
	/* perform filtering lowpass/highpass */
	fpconv(xdummyl, actual_n, h0, h1, lh, ydummyll, ydummyhh); 
	/* restore dummy variables in matrices */
	ic = -sample_f + n_c;
	for  (i=0; i<actual_n; i++){    
	  ic = ic + sample_f;
	  mat(yl, ir, ic) = ydummyll[i];  
	  mat(yh, ir, c_o_a+ic) = ydummyhh[i];  
	} 
      }
    }
      
    /* go by columns in case of a 2D signal*/
    if (m>1){
      n_rb = m/actual_m;                 /* # of row blocks per column */
      for (ic=0; ic<n; ic++){            /* loop over column */
	for (n_r=0; n_r<n_rb; n_r++){    /* loop within one column */
	  /* store in dummy variables */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    xdummyl[i] = mat(yl, ir, ic);  
	    xdummyh[i] = mat(yh, ir,c_o_a+ic);  
	  }
	  /* perform filtering: first LL/LH, then HL/HH */
	  fpconv(xdummyl, actual_m, h0, h1, lh, ydummyll, ydummylh); 
	  fpconv(xdummyh, actual_m, h0, h1, lh, ydummyhl, ydummyhh); 
	  /* restore dummy variables in matrices */
	  ir = -sample_f + n_r;
	  for (i=0; i<actual_m; i++){    
	    ir = ir + sample_f;
	    mat(yl, ir, ic) = ydummyll[i];  
	    mat(yh, ir, c_o_a+ic) = ydummylh[i];  
	    mat(yh, ir,c_o_a+n+ic) = ydummyhl[i];  
	    mat(yh, ir, c_o_a_p2n+ic) = ydummyhh[i];  
	  }
	}
      }
    }
    sample_f = sample_f*2;
  }
}

#ifdef __STDC__
fpconv(double *x_in, int lx, double *h0, double *h1, int lh,
       double *x_outl, double *x_outh)
#else
fpconv(x_in, lx, h0, h1, lh, x_outl, x_outh)
double *x_in, *h0, *h1, *x_outl, *x_outh;
int lx, lh;
#endif
{
  int i, j;
  double x0, x1;

  for (i=lx; i < lx+lh-1; i++)
    x_in[i] = x_in[i-lx];
  for (i=0; i<lx; i++){
    x0 = 0;
    x1 = 0;
    for (j=0; j<lh; j++){
      x0 = x0 + x_in[j+i]*h0[lh-1-j];
      x1 = x1 + x_in[j+i]*h1[lh-1-j];
    }
    x_outl[i] = x0;
    x_outh[i] = x1;
  }
}
