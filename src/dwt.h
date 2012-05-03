/*
 * WHAT AM I?
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

#define NORMAL_DWT 1
#define REDUNDANT_DWT 2
#define INVERSE_DWT 3
#define INVERSE_REDUNDANT_DWT 4

void dwtInit(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[],int dwtType)
{
  double *x, *h,  *y, *yl, *yh, *Lf, *Lr;
  int m, n, h_col, h_row, lh, L, i, po2, j, dim;
  double mtest, ntest;
  
  /* check for correct # of input variables */
  if (nrhs>3){
    mexErrMsgTxt("There are at most 3 input parameters allowed!");
    return;
  }
  if (nrhs<2){
    mexErrMsgTxt("There are at least 2 input parameters required!");
    return;
  }
  
  /* buffer overflow will occur if matrix isn't 1-D or 2-D */
  dim = mxGetNumberOfDimensions(prhs[0]);
  if (dim > 2){
    mexErrMsgTxt("Matrix must have fewer than 3 dimensions!");
    return;
  }
  
 
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]); 
  h = mxGetPr(prhs[1]);
  h_col = mxGetN(prhs[1]); 
  h_row = mxGetM(prhs[1]); 
  if (h_col>h_row)
    lh = h_col;
  else  
    lh = h_row;
  if (nrhs == 3){
    L = (int) *mxGetPr(prhs[2]);
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
  
  switch(dwtType) {
  case NORMAL_DWT:
    x = mxGetPr(prhs[0]);
    y = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    Lr = mxGetPr(plhs[1]);
    *Lr = L;  
    MDWT(x, m, n, h, lh, L, y);
    break;
  case REDUNDANT_DWT:
    x = mxGetPr(prhs[0]);
    yl = mxGetPr(plhs[0]);
    if (min(m,n) == 1)
        plhs[1] = mxCreateDoubleMatrix(m,L*n,mxREAL);
    else
        plhs[1] = mxCreateDoubleMatrix(m,3*L*n,mxREAL);
    yh = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    Lr = mxGetPr(plhs[2]);
    *Lr = L;
    MRDWT(x, m, n, h, lh, L, yl, yh);  
    break;
  }
}

