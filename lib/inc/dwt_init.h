#ifndef DWT_INIT_H_
#define DWT_INIT_H_

#include "matrix.h"
#include "dwt_common.h"

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

rwt_init_params dwtInit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], transform_t dwtType);

#ifdef __cplusplus
}
#endif

#endif /* DWT_INIT_H_ */
