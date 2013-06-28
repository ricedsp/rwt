%module rwt

%rename(_c_dwt) dwt;
%rename(_c_idwt) idwt;
%rename(_c_rdwt) rdwt;
%rename(_c_irdwt) irdwt;

%{
  #define SWIG_FILE_WITH_INIT
  #include "../lib/inc/rwt_transforms.h"
%}

%include "numpy.i"

%init %{
  import_array();
%}

void _c_dwt_1(  double* INPLACE_ARRAY1, int DIM1,           double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY1, int DIM1);
void _c_dwt_2(  double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY2, int DIM1, int DIM2);
void _c_idwt_1( double* INPLACE_ARRAY1, int DIM1,           double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY1, int DIM1);
void _c_idwt_2( double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY2, int DIM1, int DIM2);
void _c_rdwt_1( double* INPLACE_ARRAY1, int DIM1,           double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY1, int DIM1, double* INPLACE_ARRAY1, int DIM1);
void _c_rdwt_2( double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY2, int DIM1, int DIM2);
void _c_irdwt_1(double* INPLACE_ARRAY1, int DIM1,           double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY1, int DIM1, double* INPLACE_ARRAY1, int DIM1);
void _c_irdwt_2(double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1, int L, double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY2, int DIM1, int DIM2);

void py_swig_test1(double* INPLACE_ARRAY1, int DIM1, double* INPLACE_ARRAY1, int DIM1, double* INPLACE_ARRAY1, int DIM1);
void py_swig_test2(double* INPLACE_ARRAY2, int DIM1, int DIM2, double* INPLACE_ARRAY1, int DIM1, double* INPLACE_ARRAY2, int DIM1, int DIM2);

%inline %{

void _c_dwt_1(double *x, int m, double *h, int lh, int L, double *y, int toss1) {
  dwt(x, m, 1, h, lh, L, y);
}

void _c_idwt_1(double *x, int m, double *h, int lh, int L, double *y, int toss1) {
  idwt(x, m, 1, h, lh, L, y);
}

void _c_rdwt_1(double *x, int m, double *h, int lh, int L, double *yl, int toss1, double *yh, int toss2) {
  rdwt(x, m, 1, h, lh, L, yl, yh);
}

void _c_irdwt_1(double *x, int m, double *h, int lh, int L, double *yl, int toss1, double *yh, int toss2) {
  irdwt(x, m, 1, h, lh, L, yl, yh);
}

void _c_dwt_2(double *x, int m, int n, double *h, int lh, int L, double *y, int toss1, int toss2) {
  dwt(x, m, n, h, lh, L, y);
}

void _c_idwt_2(double *x, int m, int n, double *h, int lh, int L, double *y, int toss1, int toss2) {
  idwt(x, m, n, h, lh, L, y);
}

void _c_rdwt_2(double *x, int m, int n, double *h, int lh, int L, double *yl, int toss1, int toss2, double *yh, int toss3, int toss4) {
  rdwt(x, m, n, h, lh, L, yl, yh);
}

void _c_irdwt_2(double *x, int m, int n, double *h, int lh, int L, double *yl, int toss1, int toss2, double *yh, int toss3, int toss4) {
  irdwt(x, m, n, h, lh, L, yl, yh);
}

void py_swig_test2(double *x, int m, int n, double *h, int lh, double *z, int toss1, int toss2) {
  swig_test(x, m, n, h, lh, z);
}

void py_swig_test1(double *x, int m, double *h, int lh, double *z, int toss1) {
  swig_test(x, m, 1, h, lh, z);
}

%}

%pythoncode%{

def pst(x, h):
  z = x
  dim = len(x.shape)
  if (dim == 1):
    _rwt.py_swig_test1(x, h, z)
    return z
  if (dim == 2):
    _rwt.py_swig_test2(x, h, z)
    return z

%}
