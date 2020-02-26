/*! \file rwt_transforms.h
    \brief Function prototypes for the transform implementations
*/
#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include <math.h>
#include "rwt_init.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! dwt and rdwt take an input x and store the result in y or yl and yh
 *  idwt and irdwt take an input y or yl and yh and store the result in x
 *  In all cases it is expected that the output array has already been
 *  allocated prior to calling the transform function.
 */
void dwt_double(const double *x, double *y,const rwt_init_params * parms);
void dwt_float(const float *x, float * y,const rwt_init_params * parms);
void idwt_double(double *x, const double *y,const rwt_init_params * parms);
void idwt_float(float *x, const float *y,const rwt_init_params * parms);
void  rdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *yl, double *yh);
void irdwt(double *x, size_t nrows, size_t ncols, double *h, int ncoeff, int levels, double *yl, double *yh);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFORMS_H_ */
