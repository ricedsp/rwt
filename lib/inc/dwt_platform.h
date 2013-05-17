#ifndef platform_h
#define platform_h

#include "dwt_common.h"
#include <stdlib.h>

#define mat(a, i, j, m) (*(a + (m*(j)+i)))  /* macro for matrix indices */

#ifndef EXTERN_C
  #ifdef __cplusplus
    #define EXTERN_C extern "C"
  #else
    #define EXTERN_C extern
  #endif
#endif

EXTERN_C void *rwt_malloc(size_t size);
EXTERN_C void *rwt_calloc(size_t num, size_t size);
EXTERN_C void rwt_free(void *ptr);

#endif
