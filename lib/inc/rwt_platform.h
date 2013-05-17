/*! \file rwt_platform.h
    \brief Header for memory allocation wrapper functions in platform.c
*/
#ifndef RWT_PLATFORM_H
#define RWT_PLATFORM_H

#include "rwt_common.h"
#include <stdlib.h>

#define mat(a, i, j, m) (*(a + (m*(j)+i)))  /* macro for matrix indices */

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
