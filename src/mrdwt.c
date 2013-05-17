/*
File Name: mrdwt.c
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

/*! \file mrdwt.c
    \brief MATLAB gateway for the redundant discrete wavelet transform

    This file is used to produce a MATLAB MEX binary for the redundant discrete wavelet transform
*/

#include "mex.h"
#include "matrix.h"
#include "rwt_init.h"
#include "rwt_transforms.h"

/*!
 * Matlab MEX definition for the redundant discrete wavelet transform.
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *x, *yl, *yh;
  rwt_init_params params = rwt_matlab_init(nlhs, plhs, nrhs, prhs, REDUNDANT_DWT);
  if (min(params.nrows, params.ncols) == 1)
    plhs[1] = mxCreateDoubleMatrix(params.nrows, params.levels*params.ncols, mxREAL);
  else
    plhs[1] = mxCreateDoubleMatrix(params.nrows, 3*params.levels*params.ncols, mxREAL);
  x = mxGetPr(prhs[0]);
  yl = mxGetPr(plhs[0]);
  yh = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[2]) = params.levels;
  RDWT(x, params.nrows, params.ncols, params.scalings, params.lh, params.levels, yl, yh);
}

