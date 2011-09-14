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

#include <math.h>
#include <stdio.h>

#define max(A,B) (A > B ? A : B)
#define mat(a, i, j) (*(a + (m*(j)+i)))  /* macro for matrix indices */

#ifdef __STDC__
MDWT(double *x, int m, int n, double *h, int lh, int L, double *y)
#else
MDWT(x, m, n, h, lh, L, y)
double *x, *h, *y;
int m, n, lh, L;
#endif
{
  double  *h0, *h1, *ydummyl, *ydummyh, *xdummy;
  long i, j;
  int actual_L, actual_m, actual_n, r_o_a, c_o_a, ir, ic, lhm1;
  xdummy = (double *)mxCalloc(max(m,n)+lh-1,sizeof(double));
  ydummyl = (double *)mxCalloc(max(m,n),sizeof(double));
  ydummyh = (double *)mxCalloc(max(m,n),sizeof(double));
  h0 = (double *)mxCalloc(lh,sizeof(double));
  h1 = (double *)mxCalloc(lh,sizeof(double));
  
  
  /* analysis lowpass and highpass */
  if (n==1){
    n = m;
    m = 1;
  }
  for (i=0; i<lh; i++){
    h0[i] = h[lh-i-1];
    h1[i] =h[i];
  }
  for (i=0; i<lh; i+=2)
    h1[i] = -h1[i];
  
  lhm1 = lh - 1;
  actual_m = 2*m;
  actual_n = 2*n;
  
  /* main loop */
  for (actual_L=1; actual_L <= L; actual_L++){
    if (m==1)
      actual_m = 1;
    else{
      actual_m = actual_m/2;
      r_o_a = actual_m/2;     
    }
    actual_n = actual_n/2;
    c_o_a = actual_n/2;
    
    /* go by rows */
    for (ir=0; ir<actual_m; ir++){            /* loop over rows */
      /* store in dummy variable */
      for (i=0; i<actual_n; i++)
	if (actual_L==1)  
	  xdummy[i] = mat(x, ir, i);  
	else 
	  xdummy[i] = mat(y, ir, i);  
      /* perform filtering lowpass and highpass*/
      fpsconv(xdummy, actual_n, h0, h1, lhm1, ydummyl, ydummyh); 
      /* restore dummy variables in matrices */
      ic = c_o_a;
      for  (i=0; i<c_o_a; i++){    
	mat(y, ir, i) = ydummyl[i];  
	mat(y, ir, ic++) = ydummyh[i];  
      } 
    }  
    
    /* go by columns in case of a 2D signal*/
    if (m>1){
      for (ic=0; ic<actual_n; ic++){            /* loop over column */
	/* store in dummy variables */
	for (i=0; i<actual_m; i++)
	  xdummy[i] = mat(y, i, ic);  
	/* perform filtering lowpass and highpass*/
	fpsconv(xdummy, actual_m, h0, h1, lhm1, ydummyl, ydummyh); 
	/* restore dummy variables in matrix */
	ir = r_o_a;
	for (i=0; i<r_o_a; i++){    
	  mat(y, i, ic) = ydummyl[i];  
	  mat(y, ir++, ic) = ydummyh[i];  
	}
      }
    }
  }
}

#ifdef __STDC__
fpsconv(double *x_in, int lx, double *h0, double *h1, int lhm1, 
	double *x_outl, double *x_outh)
#else
fpsconv(x_in, lx, h0, h1, lhm1, x_outl, x_outh)
double *x_in, *h0, *h1, *x_outl, *x_outh;
int lx, lhm1;
#endif

{
  int i, j, ind;
  double x0, x1;

  for (i=lx; i < lx+lhm1; i++)
    x_in[i] = *(x_in+(i-lx));
  ind = 0;
  for (i=0; i<(lx); i+=2){
    x0 = 0;
    x1 = 0;
    for (j=0; j<=lhm1; j++){
      x0 = x0 + x_in[i+j]*h0[lhm1-j];
      x1 = x1 + x_in[i+j]*h1[lhm1-j];
    }
    x_outl[ind] = x0;
    x_outh[ind++] = x1;
  }
}
