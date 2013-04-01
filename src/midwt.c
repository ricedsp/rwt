/*
File Name: midwt.c
Last Modification Date:	06/14/95	12:55:58
Current Version: midwt.c	1.4
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
  double *x, *y, *Lr;
  rwt_init_params params = dwtInit(nlhs, plhs, nrhs, prhs, INVERSE_DWT);
  y = mxGetPr(prhs[0]);
  x = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  Lr = mxGetPr(plhs[1]);
  *Lr = params.L;
  MIDWT(x, params.m, params.n, params.h, params.lh, params.L, y);
}

