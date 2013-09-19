/*! \file rwt_platform.h
    \brief Header for memory allocation wrapper functions in platform.c
*/
#ifndef RWT_PLATFORM_H
#define RWT_PLATFORM_H

#include "rwt_common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/*! For MATLAB we address 2d inputs and outputs in column-major order */
/*! For Python we address 2d inputs and outputs in row-major order */
#ifdef MATLAB_MEX_FILE
  #include "matrix.h"
  #include "mex.h"
  #define mat(a, i, j, m, n) (*(a + (m*(j)+i)))
  #define mat_offset(a, i, j, m, n) (m*(j)+i)
  #define rwt_printf(fmt, ...) mexPrintf(fmt, ##__VA_ARGS__)
  #define offset_row(offset, m, n) (offset % m)
  #define offset_col(offset, m, n) ((offset - (offset % m)) / m)
#else
  #define mat(a, i, j, m, n) (*(a + (n*(i)+j)))
  #define mat_offset(a, i, j, m, n) (n*(i)+j)
  #define rwt_printf(fmt, ...) printf(fmt, ##__VA_ARGS__)
  #define offset_row(offset, m, n) ((offset - (offset % n)) / n)
  #define offset_col(offset, m, n) (offset % n)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void *rwt_malloc(size_t size);
void *rwt_calloc(size_t num, size_t size);
void rwt_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
