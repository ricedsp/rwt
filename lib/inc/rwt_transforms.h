/*! \file rwt_transforms.h
    \brief Function prototypes for the transform implementations
*/
#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

void dwt(double *x, int m, int n, double *h, int lh, int L, double *y);
void idwt(double *x, int m, int n, double *h, int lh, int L, double *y);
void rdwt(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);
void irdwt(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFORMS_H_ */
