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

#include "mex.h"
#include "matrix.h"
#include "dwt_init.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *x, *yl, *yh;
  rwt_init_params params = dwtInit(nlhs, plhs, nrhs, prhs, REDUNDANT_DWT);
  x = mxGetPr(prhs[0]);
  yl = mxGetPr(plhs[0]);
  if (min(params.nrows, params.ncols) == 1)
    plhs[1] = mxCreateDoubleMatrix(params.nrows, params.levels*params.ncols, mxREAL);
  else
    plhs[1] = mxCreateDoubleMatrix(params.nrows, 3*params.levels*params.ncols, mxREAL);
  yh = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[2]) = params.levels;
  MRDWT(x, params.nrows, params.ncols, params.scalings, params.lh, params.levels, yl, yh);
}

