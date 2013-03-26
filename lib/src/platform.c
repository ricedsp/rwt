#include "platform.h"
#ifdef MATLAB_MEX_FILE
 #include "matrix.h"
 void *rmalloc(size_t size) {
   return mxMalloc(size);
 }
 void *rcalloc(size_t num, size_t size) {
   return mxCalloc(num, size);
 }
 void rfree(void *ptr) {
   return mxFree(ptr);
 }
#else
 void *rmalloc(size_t size) {
   return malloc(size);
 }
 void *rcalloc(size_t num, size_t size) {
   return calloc(num, size);
 }
 void rfree(void *ptr) {
   free(ptr);
 }
#endif
