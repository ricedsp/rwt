/*! \file rwt_init.h
    \brief Header for matlab init functions in init.c
*/
#ifndef RWT_INIT_H_
#define RWT_INIT_H_

#include "rwt_platform.h"

#ifdef MATLAB_MEX_FILE
  #include "mex.h"
  #ifndef HAVE_OCTAVE
    #include "matrix.h"
  #endif
#endif

typedef struct {
  int nrows;        /*!< The number of rows in the input matrix. Output matrix will match.  */
  int ncols;        /*!< The number of columns in the input matrix. Output matrix will match. */
  int levels;       /*!< L, the number of levels for the transform. */
  int lh;           /*!< Length of h / the number of scaling coefficients */
  double *scalings; /*!< Wavelet scaling coefficients */
} rwt_init_params;

typedef enum {NORMAL_DWT, REDUNDANT_DWT, INVERSE_DWT, INVERSE_REDUNDANT_DWT} transform_t;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MATLAB_MEX_FILE
  rwt_init_params rwt_matlab_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t dwtType);
#endif

#ifdef __cplusplus
}
#endif

#endif /* RWT_INIT_H_ */
