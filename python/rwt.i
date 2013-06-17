%module rwt

%{
  #define SWIG_FILE_WITH_INIT
  #include "../lib/inc/rwt_transforms.h"
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* x, int m, int n)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* h, int lh)};
%apply (double* INPLACE_ARRAY2) {(double* y)};
%apply (double* INPLACE_ARRAY2) {(double* yl)};
%apply (double* INPLACE_ARRAY2) {(double* yh)};

%include "../lib/inc/rwt_transforms.h"

