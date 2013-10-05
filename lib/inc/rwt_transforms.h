/*! \file rwt_transforms.h
    \brief Function prototypes for the transform implementations
*/
#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

void   dwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y);
void  idwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *y);
void  rdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *yl, double *yh);
void irdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *yl, double *yh);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFORMS_H_ */
