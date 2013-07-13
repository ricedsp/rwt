/*! \file rwt_platform.h
    \brief Header for memory allocation wrapper functions in platform.c
*/
#ifndef RWT_PLATFORM_H
#define RWT_PLATFORM_H

#include "rwt_common.h"
#include <stdlib.h>

/*! For MATLAB we need to address inputs and outputs in column-major order */
#ifdef MATLAB_MEX_FILE
  #define mat(a, i, j, m, n) (*(a + (m*(j)+i)))
#else
  #define mat(a, i, j, m, n) (*(a + (n*(i)+j)))
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
