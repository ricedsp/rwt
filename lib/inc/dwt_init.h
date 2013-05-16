/*
 * dwt_init.h
 *
 *  Created on: Jul 5, 2012
 *      Author: robert
 */

#ifndef DWT_INIT_H_
#define DWT_INIT_H_

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#include "limits.h"

typedef struct {
  int nrows;        /*!< The number of rows in the input matrix  */
  int ncols;        /*!< The number of columns in the input matrix  */
  int lh;
  int levels;
  double *scalings; /*!< Wavelet scaling coefficients */
} rwt_init_params;

#define max(A,B) (A > B ? A : B)
#define min(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)

typedef enum {NORMAL_DWT, REDUNDANT_DWT, INVERSE_DWT, INVERSE_REDUNDANT_DWT} transform_t;

#ifdef __cplusplus
extern "C" {
#endif

int MDWT(double *x, int m, int n, double *h, int lh, int L, double *y);
int MIDWT(double *x, int m, int n, double *h, int lh, int L, double *y);
int MRDWT(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);
int MIRDWT(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);

int dwtInputCheck(int nrhs, transform_t dwtType);
int dwtEstimateL(int n, int m);
rwt_init_params dwtInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t dwtType);

#ifdef __cplusplus
}
#endif

#endif /* DWT_INIT_H_ */
