function [y,L] = mdwt(x,h,L)
%    [y,L] = mdwt(x,h[,L[,transdims]]);
%
%    Function computes the discrete wavelet transform y for a 1D or 2D input
%    signal x using the scaling filter h.
%
%    Input:
%       x : input matrix or nd array ( transform will be done on leading 1 or 2 dims)
%       h : scaling filter
%       L : number of levels. In the case of a 1D transform, size(x,1) must be
%           divisible by 2^L; for a 2D transform, size(x,2) must also be
%           divisible by 2^L. The default is the maximal possible L.
%       transdims: 1 or 2 dimensional transform, default is 2 (if size allows )
%
%    Output:
%       y : the wavelet transform of the signal
%           (see example to understand the coefficients)
%       L : number of decomposition levels
%
%    1D Example:
%       x = makesig('LinChirp',8);
%       h = daubcqf(4,'min');
%       L = 2;
%       [y,L] = mdwt(x,h,L)
%
%    1D Example's  output and explanation:
%
%       y = [1.1097 0.8767 0.8204 -0.5201 -0.0339 0.1001 0.2201 -0.1401]
%       L = 2
%
%    The coefficients in output y are arranged as follows
%
%       y(1) and y(2) : Scaling coefficients (lowest frequency)
%       y(3) and y(4) : Band pass wavelet coefficients
%       y(5) to y(8)  : Finest scale wavelet coefficients (highest frequency)
%
%    2D Example:
%
%       load test_image
%       h = daubcqf(4,'min');
%       L = 1;
%       [y,L] = mdwt(test_image,h,L);
%
%    2D Example's  output and explanation:
%
%       The coefficients in y are arranged as follows.
%
%              .------------------.
%              |         |        |
%              |    4    |   2    |
%              |         |        |
%              |   L,L   |   H,L  |
%              |         |        |
%              --------------------
%              |         |        |
%              |    3    |   1    |
%              |         |        |
%              |   L,H   |  H,H   |
%              |         |        |
%              `------------------'
%
%       where
%            1 : High pass vertically and high pass horizontally
%            2 : Low pass vertically and high pass horizontally
%            3 : High pass vertically and low  pass horizontally
%            4 : Low pass vertically and Low pass horizontally
%                (scaling coefficients)
%
%    4D Tensor with 2D transforms across leading  Example:
%
%       x=randn(64,64,3,2);
%       h=[1,1]*sqrt(.5);
%       y=mdwt(x,h);
%       y22 = mdwt(x(:,:,2,2),h);
%       assert (norm(y(:,:,2,2) - y22 ) < 1e-9)
%
%    See also: midwt, mrdwt, mirdwt
%
%Author: Markus Lang  <lang@jazz.rice.edu>
if exist('OCTAVE_VERSION', 'builtin')
  x = x * 1.0;
  if (exist('L'))
    [y,L] = omdwt(x,h,L);
  else
    [y,L] = omdwt(x,h);
  end
else
  error('You must compile wavelet toolbox before use')
end
