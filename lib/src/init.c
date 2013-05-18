/*! \file init.c
    \brief Parse input from MATLAB and do some sanity checking

*/

#include "rwt_init.h"
#include <math.h>


/*!
 * Checks for correct # of input variables based on type of transform.
 *
 * @param nrhs number of items on right hand side of matlab call
 * @param transform_type 
 *
 */
int rwt_check_parameter_count(int nrhs, transform_t transform_type) {
  if (transform_type == INVERSE_REDUNDANT_DWT) {
    if (nrhs > 4) {
      mexErrMsgTxt("There are at most 4 input parameters allowed!");
      return 1;
    }
    if (nrhs < 3) {
      mexErrMsgTxt("There are at least 3 input parameters required!");
      return 1;
    }
  }
  else {
    if (nrhs > 3) {
      mexErrMsgTxt("There are at most 3 input parameters allowed!");
      return 1;
    }
    if (nrhs < 2) {
      mexErrMsgTxt("There are at least 2 input parameters required!");
      return 1;
    }
  }
  return 0;
}


/*!
 * For the inverse redundant transform check that the dimensions of the low and high inputs match
 *
 * @param prhs
 * @param params
 *
 */
int rwt_check_yl_matches_yh(const mxArray *prhs[], int nrows, int ncols, int levels) {
  int mh = mxGetM(prhs[1]);
  int nh = mxGetN(prhs[1]);
  if (min(nrows, ncols) > 1){
    if ((nrows != mh) | (3 * ncols * levels != nh)) {
      return 0;
    }
  }
  else {
    if ((nrows != mh) | (ncols * levels != nh)) {
      return 0;
    }
  }
  return 1;
}


/*!
 * Find L, the number of levels
 *
 * @param m the number of rows in the input
 * @param n the number of columns in the input
 *
 */
int rwt_find_L(int m, int n) {
  int i, j, L;
  i = n ; j = 0;
  while (even(i)) {
    i = (i >> 1);
    j++;
  }
  L = m; i = 0;
  while (even(L)) {
    L = (L >> 1);
    i++;
  }
  if (min(m, n) == 1)
    L = max(i, j);
  else
    L = min(i, j);
  if (L == 0) {
    mexErrMsgTxt("Maximum number of levels is zero; no decomposition can be performed!");
    return -1;
  }
  else return L;
}


/*!
 * Check that length is divisble by 2^L
 *
 * @param length pass in the number of rows or number if columns
 * @param L the number of levels
 *
 */
int dimensionCheck(int length, int L) {
  double test = (double) length / pow(2.0, (double) L);
  if (!isint(test)) {
    mexErrMsgTxt("The matrix dimensions must be of size m*2^(L) by n*2^(L)");
    return 1;
  }
  return 0;
}


/*!
 * Parse input from MATLAB and do some sanity checking
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 * @param transform_type which transform are we setting up to do
 *
 */
rwt_init_params rwt_matlab_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t transform_type) {
  rwt_init_params params;
  int argNumL;

  /* check for correct # of input variables */
  if (rwt_check_parameter_count(nrhs, transform_type) != 0) return;

  /* buffer overflow will occur if matrix isn't 1-D or 2-D */
  if (mxGetNumberOfDimensions(prhs[0]) > 2) {
    mexErrMsgTxt("Matrix must have fewer than 3 dimensions!");
    return;
  }

  /* Get input matrix row and column number */
  params.nrows = mxGetM(prhs[0]);
  params.ncols = mxGetN(prhs[0]);

  /* Read L from command line or compute L */
  argNumL = (transform_type == INVERSE_REDUNDANT_DWT) ? 3 : 2;
  if ((argNumL + 1) == nrhs)
    params.levels = (int) *mxGetPr(prhs[argNumL]);
  else
    params.levels = rwt_find_L(params.nrows, params.ncols);

  if (params.levels < 0) {
    mexErrMsgTxt("The number of levels, L, must be a non-negative integer");
    return;
  }

  /* Check input dimensions */
  if ((params.nrows > 1 && dimensionCheck(params.nrows, params.levels)) || (params.ncols > 1 && dimensionCheck(params.ncols, params.levels)))
    return;

  if (transform_type == INVERSE_REDUNDANT_DWT) {
    params.scalings = mxGetPr(prhs[2]);
    params.lh = max(mxGetM(prhs[2]), mxGetN(prhs[2]));
    if (!rwt_check_yl_matches_yh(prhs, params.nrows, params.ncols, params.levels)) {
      mexErrMsgTxt("Dimensions of first two input matrices not consistent!");
      return;
    }
  }
  else {
    params.scalings = mxGetPr(prhs[1]);
    params.lh = max(mxGetM(prhs[1]), mxGetN(prhs[1]));
  }

  /*! Create the first item in the output array as a double matrix with the same dimensions as the input */
  plhs[0] = mxCreateDoubleMatrix(params.nrows, params.ncols, mxREAL);

  return params;
}


