#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

int MDWT(double *x, int m, int n, double *h, int lh, int L, double *y);
int MIDWT(double *x, int m, int n, double *h, int lh, int L, double *y);
int MRDWT(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);
int MIRDWT(double *x, int m, int n, double *h, int lh, int L, double *yl, double *yh);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFORMS_H_ */
