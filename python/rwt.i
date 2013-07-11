%module rwt

%rename(_c_dwt)     dwt;
%rename(_c_idwt)   idwt;
%rename(_c_rdwt)   rdwt;
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

%pythoncode %{

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

def hard_th(y, thld):
  return (abx(x) > thld) * x

def soft_th(y, thld):
  x = abs(y)
  return sign(y) * (x >= thld) * (x - thld)

def makesig(signame, n):
  t = array(range(1, n + 1)) / float(n)
  if (signame == 'HeaviSine'):
    y = 4 * sin(4 * pi * t)
    return y - sign(t - .3) - sign(.72 - t)
  if (signame == 'Bumps'):
    pos = array([.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81])
    hgt = array([4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2])
    wth = array([.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005])
    y = zeros(n)
    for j in range(0, pos.size):
      y = y + hgt[j] / pow((1 + abs((t - pos[j]) / wth[j])), 4)
    return y
  if (signame == 'Blocks'):
    pos = array([.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81])
    hgt = array([4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2)])
    y = zeros(n)
    for j in range(0, pos.size):
      y = y + (1 + sign(t - pos[j])) * (hgt[j]/2)
    return y
  if (signame == 'Doppler'):
    return sqrt(t * (1-t)) * sin((2 * pi * 1.05) / (t+.05))
  if (signame == 'Ramp'):
    return t - (t >= .37)
  if (signame == 'Cusp'):
    return sqrt(abs(t - .37))
  if (signame == 'Sing'):
    k = floor(n * .37)
    return 1 / abs(t - (k + .5)/n)
  if (signame == 'HiSine'):
    return sin(pi * (n * .6902) * t)
  if (signame == 'LoSine'):
    return sin(pi * (n * .3333) * t)
  if (signame == 'LinChirp'):
    return sin(pi * t * ((n * .125) * t))
  if (signame == 'TwoChirp'):
    return sin(pi * t * (n * t)) + sin((pi / 3) * t * (n * t))
  if (signame == 'QuadChirp'):
    return sin((pi/3) * t * (n * pow(t,2)))
  if (signame == 'MishMash'):
    y = sin((pi/3) * t * (n * pow(t,2)))
    y = y + sin(pi * (n * .6902) * t)
    return y + sin(pi * t * (n * .125 * t))
  if (signame == 'WernerSorrows'):
    y = sin(pi * t * (n/2 * pow(t, 2)))
    y = y + sin(pi * (n * .6902) * t)
    y = y + sin(pi * t * (n * t))
    pos = array([.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81])
    hgt = array([4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2])
    wth = array([.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005])
    for j in range(0, pos.size):
      y = y + hgt[j] / pow((1 + abs((t - pos[j]) / wth[j])), 4)
    return y
  if (signame == 'Leopold'):
    return (t == floor(.37 * n)/n) * 1

%}
