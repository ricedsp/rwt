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

%}


%pythoncode%{

from numpy import *

def dwt(x, h, L):
  y = x
  dim = len(x.shape)
  if (dim == 1):
    _rwt._c_dwt_1(x, h, L, y)
  if (dim == 2):
    _rwt._c_dwt_2(x, h, L, y)
  return y, L

def idwt(y, h, L):
  x = y
  dim = len(x.shape)
  if (dim == 1):
    _rwt._c_idwt_1(x, h, L, y)
  if (dim == 2):
    _rwt._c_idwt_2(x, h, L, y)
  return x, L

def rdwt(x, h, L):
  yl = yh = x
  dim = len(x.shape)
  if (dim == 1):
    _rwt._c_rdwt_1(x, h, L, yl, yh)
  if (dim == 2):
    _rwt._c_rdwt_2(x, h, L, yl, yh)
  return yl, yh, L

def irdwt(y, h, L):
  x = y
  dim = len(x.shape)
  if (dim == 1):
    _rwt._c_irdwt_1(x, h, L, yl, yh)
  if (dim == 2):
    _rwt._c_irdwt_2(x, h, L, yl, yh)
  return x, L

def daubcqf(n, dtype = 'min'):
  if (n % 2 != 0):
    raise Exception("No Daubechies filter exists for ODD length")
  k = n / 2
  a = p = q = 1
  h_0 = array([1, 1])
  for j in range(1, k):
    a = -a * 0.25 * (j + k - 1) / j
    h_0 = hstack((0, h_0)) + hstack((h_0, 0))
    p = hstack((0, -p)) + hstack((p, 0))
    p = hstack((0, -p)) + hstack((p, 0))
    q = hstack((0, q, 0)) + a*p
  q = sort(roots(q))
  qt = q[0:k-1]
  if (dtype == 'mid'):
    raise Exception("Not implemented") # MATLAB code for this is opaque
    #if (k % 2 == 1):
    #  qt = hstack((q[0:n-2:4], q[1:n-2:4]))
    #else:
    #  qt = hstack((q[0], q[3:k-1:4], q[4:k-1:4], q[n-4:k-1:-4], q[n-5:k-1:-4]))
  h_0 = convolve(h_0, real(poly(qt)))
  h_0 = sqrt(2)*h_0 / sum(h_0)
  if (dtype == 'max'):
    h_0 = flipud(h_0)
  if (abs(sum(power(h_0, 2))) -1 > 1e-4):
    raise Exception("Numerically unstable for this value of n")
  h_1 = copy(flipud(h_0))
  h_1[0:n-1:2] = -h_1[0:n-1:2]
  return h_0, h_1

%}
