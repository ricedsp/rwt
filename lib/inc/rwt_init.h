/*! \file rwt_init.h
    \brief Header for matlab init functions in init.c
*/
#ifndef RWT_INIT_H_
#define RWT_INIT_H_

#include "matrix.h"
#include "rwt_common.h"

typedef struct {
  int nrows;        /*!< The number of rows in the input matrix  */
  int ncols;        /*!< The number of columns in the input matrix  */
  int lh;
  int levels;
  double *scalings; /*!< Wavelet scaling coefficients */
} rwt_init_params;

typedef enum {NORMAL_DWT, REDUNDANT_DWT, INVERSE_DWT, INVERSE_REDUNDANT_DWT} transform_t;

#ifdef __cplusplus
extern "C" {
#endif

rwt_init_params rwt_matlab_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t dwtType);

#ifdef __cplusplus
}
#endif

#endif /* RWT_INIT_H_ */
