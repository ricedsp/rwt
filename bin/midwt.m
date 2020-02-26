function [y,L] = midwt(x,h,L)
%    [x,L] = midwt(y,h[,L[,transdims]])
%
%    Compute the inverse discrete wavelet transform x for a 1D or
%    2D input signal (or tensor of signals) y using the scaling filter h.
%
%    Input:
%       y : input wavelet domain coefficients (transform is done over leading dims)
%           (see function mdwt to find the structure of y)
%       h : scaling filter
%       L : number of levels. In the case of a 1D transform, size(x,1) must be
%           divisible by 2^L; for a 2D transform, size(x,2) must also be
%           divisible by 2^L. The default is the maximal possible L.
%       transdims: 1 or 2 dimensional transform, default is 2 (if size allows )
%
%    Output:
%       x : periodic reconstructed signal
%       L : number of decomposition levels
%
%    1D Example:
%       xin = makesig('LinChirp',8);
%       h = daubcqf(4,'min');
%       L = 1;
%       [y,L] = mdwt(xin,h,L);
%       [x,L] = midwt(y,h,L)
%
%    1D Example's  output:
%
%       x = 0.0491 0.1951 0.4276 0.7071 0.9415 0.9808 0.6716 0.0000
%       L = 1
%
%    See also: mdwt, mrdwt, mirdwt
%
%Author: Markus Lang  <lang@jazz.rice.edu>
if exist('OCTAVE_VERSION', 'builtin')
  if (exist('L'))
    [y,L] = omidwt(x,h,L);
  else
    [y,L] = omidwt(x,h);
  end
else
  error('You must compile wavelet toolbox before use')
end
