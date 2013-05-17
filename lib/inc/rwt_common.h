/*! \file rwt_common.h
    \brief Define some macros used throughout the code
*/
#ifndef RWT_COMMON_H_
#define RWT_COMMON_H_

#define max(A,B) (A > B ? A : B)
#define min(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)

#endif /* RWT_COMMON_H_ */
