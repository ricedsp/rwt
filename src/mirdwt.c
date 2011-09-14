/*
File Name: mirdwt.c
Last Modification Date:	%G%	%U%
Current Version: %M%	%I%
File Creation Date: Wed Oct 12 08:44:43 1994
Author: Markus Lang  <lang@jazz.rice.edu>

Copyright: All software, documentation, and related files in this distribution
           are Copyright (c) 1994  Rice University

Permission is granted for use and non-profit distribution providing that this
notice be clearly maintained. The right to distribute any portion for profit
or as part of any commercial product is specifically reserved for the author.

Change History: Fixed code such that the result has the same dimension as the 
                input for 1D problems. Also, added some standard error checking.
		Jan Erik Odegard <odegard@ece.rice.edu> Wed Jun 14 1995

*/

#include <math.h>
/*#include <malloc.h>*/
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

#define max(A,B) (A > B ? A : B)
#define min(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])

{
  double *x, *h,  *yl, *yh, *Lf, *Lr;
  int m, n, mh, nh, h_col, h_row, lh, L, i, po2, j;
  double mtest, ntest;

  /* check for correct # of input variables */
  if (nrhs>4){
    mexErrMsgTxt("There are at most 4 input parameters allowed!");
    return;
  }
  if (nrhs<3){
    mexErrMsgTxt("There are at least 3 input parameters required!");
    return;
  }
  yl = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]); 
  yh = mxGetPr(prhs[1]);
  nh = mxGetN(prhs[1]); 
  mh = mxGetM(prhs[1]); 
  h = mxGetPr(prhs[2]);
  h_col = mxGetN(prhs[2]); 
  h_row = mxGetM(prhs[2]); 
  if (h_col>h_row)
    lh = h_col;
  else  
    lh = h_row;
  if (nrhs == 4){
    L = (int) *mxGetPr(prhs[3]);
    if (L < 0)
      mexErrMsgTxt("The number of levels, L, must be a non-negative integer");
  }
  else /* Estimate L */ {
    i=n;j=0;
    while (even(i)){
      i=(i>>1);
      j++;
    }
    L=m;i=0;
    while (even(L)){
      L=(L>>1);
      i++;
    }
    if(min(m,n) == 1)
      L = max(i,j);
    else
      L = min(i,j);
    if (L==0){
      mexErrMsgTxt("Maximum number of levels is zero; no decomposition can be performed!");
      return;
    }
  }
  /* check for consistency of rows and columns of yl, yh */
  if (min(m,n) > 1){
    if((m != mh) | (3*n*L != nh)){
      mexErrMsgTxt("Dimensions of first two input matrices not consistent!");
      return;
    }
  }
  else{
    if((m != mh) | (n*L != nh)){
      mexErrMsgTxt("Dimensions of first two input vectors not consistent!");{
	return;
      }
    }
  }
  /* Check the ROW dimension of input */
  if(m > 1){
    mtest = (double) m/pow(2.0, (double) L);
    if (!isint(mtest))
      mexErrMsgTxt("The matrix row dimension must be of size m*2^(L)");
  }
  /* Check the COLUMN dimension of input */
  if(n > 1){
    ntest = (double) n/pow(2.0, (double) L);
    if (!isint(ntest))
      mexErrMsgTxt("The matrix column dimension must be of size n*2^(L)");
  }
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  x = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  Lr = mxGetPr(plhs[1]);
  *Lr = L;
  MIRDWT(x, m, n, h, lh, L, yl, yh);
}
