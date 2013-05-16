/*
File Name: mdwt.c
Last Modification Date:	06/14/95	12:56:43
Current Version: mdwt.c	1.5
File Creation Date: Wed Oct 12 08:44:43 1994
Author: Markus Lang  <lang@jazz.rice.edu>

Copyright: All software, documentation, and related files in this distribution
           are Copyright (c) 1994 Rice University

Permission is granted for use and non-profit distribution providing that this
notice be clearly maintained. The right to distribute any portion for profit
or as part of any commercial product is specifically reserved for the author.

Change History: Fixed code such that the result has the same dimension as the 
                input for 1D problems. Also, added some standard error checking.
		Jan Erik Odegard <odegard@ece.rice.edu> Wed Jun 14 1995

MATLAB gateway for MDWT.c, discrete wavelet transform
*/

#include "mex.h"
#include "matrix.h"
#include "dwt_init.h"

/*!
 * Matlab MEX definition for the discrete wavelet transform.
 *
 * @param nlhs number of items on left hand side of matlab call
 * @param plhs pointer to left hand side data structure
 * @param nrhs number of items on right hand side of matlab call
 * @param prhs pointer to right hand side data structure
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  rwt_init_params params = dwtInit(nlhs, plhs, nrhs, prhs, NORMAL_DWT);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[1]) = params.levels; /*! The second returned item is the number of levels */
  MDWT(mxGetPr(prhs[0]), params.nrows, params.ncols, params.scalings, params.lh, params.levels, mxGetPr(plhs[0]));
}

