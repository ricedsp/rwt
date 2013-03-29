#include "platform.h"
#ifdef MATLAB_MEX_FILE
  #include "matrix.h"
  void *rwt_malloc(size_t size) {
    return mxMalloc(size);
  }
  void *rwt_calloc(size_t num, size_t size) {
    return mxCalloc(num, size);
  }
  void rwt_free(void *ptr) {
    mxFree(ptr);
  }
#else
  void *rwt_malloc(size_t size) {
    return malloc(size);
  }
  void *rwt_calloc(size_t num, size_t size) {
    return calloc(num, size);
  }
  void rwt_free(void *ptr) {
    free(ptr);
  }
#endif
