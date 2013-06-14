function [h_0,h_1] = daubcqf(N,TYPE)
%    [h_0,h_1] = daubcqf(N,TYPE); 
%
%    Function computes the Daubechies' scaling and wavelet filters
%    (normalized to sqrt(2)).
%
%    Input: 
%       N    : Length of filter (must be even)
%       TYPE : Optional parameter that distinguishes the minimum phase,
%              maximum phase and mid-phase solutions ('min', 'max', or
%              'mid'). If no argument is specified, the minimum phase
%              solution is used.
%
%    Output: 
%       h_0 : Minimal phase Daubechies' scaling filter 
%       h_1 : Minimal phase Daubechies' wavelet filter 
%
%    Example:
%       N = 4;
%       TYPE = 'min';
%       [h_0,h_1] = daubcqf(N,TYPE)
%       h_0 = 0.4830 0.8365 0.2241 -0.1294
%       h_1 = 0.1294 0.2241 -0.8365 0.4830
%
%    Reference: "Orthonormal Bases of Compactly Supported Wavelets",
%                CPAM, Oct.89 
%

%File Name: daubcqf.m
%Last Modification Date: 01/02/96	15:12:57
%Current Version: daubcqf.m	2.4
%File Creation Date: 10/10/88
%Author: Ramesh Gopinath  <ramesh@dsp.rice.edu>
%
%Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
%Created by Ramesh Gopinath, Department of ECE, Rice University. 
%
%This software is distributed and licensed to you on a non-exclusive 
%basis, free-of-charge. Redistribution and use in source and binary forms, 
%with or without modification, are permitted provided that the following 
%conditions are met:
%
%1. Redistribution of source code must retain the above copyright notice, 
%   this list of conditions and the following disclaimer.
%2. Redistribution in binary form must reproduce the above copyright notice, 
%   this list of conditions and the following disclaimer in the 
%   documentation and/or other materials provided with the distribution.
%3. All advertising materials mentioning features or use of this software 
%   must display the following acknowledgment: This product includes 
%   software developed by Rice University, Houston, Texas and its contributors.
%4. Neither the name of the University nor the names of its contributors 
%   may be used to endorse or promote products derived from this software 
%   without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS, 
%AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
%BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
%FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY 
%OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
%OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
%WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
%OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE 
%USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%For information on commercial licenses, contact Rice University's Office of 
%Technology Transfer at techtran@rice.edu or (713) 348-6173

if(nargin < 2),
  TYPE = 'min';
end;
if(rem(N,2) ~= 0),
  error('No Daubechies filter exists for ODD length');
end;
K = N/2;
a = 1;
p = 1;
q = 1;
h_0 = [1 1];
for j  = 1:K-1,
  a = -a * 0.25 * (j + K - 1)/j;
  h_0 = [0 h_0] + [h_0 0];
  p = [0 -p] + [p 0];
  p = [0 -p] + [p 0];
  q = [0 q 0] + a*p;
end;
q = sort(roots(q));
qt = q(1:K-1);
if TYPE=='mid',
  if rem(K,2)==1,  
    qt = q([1:4:N-2 2:4:N-2]);
  else
    qt = q([1 4:4:K-1 5:4:K-1 N-3:-4:K N-4:-4:K]);
  end;
end;
h_0 = conv(h_0,real(poly(qt)));
h_0 = sqrt(2)*h_0/sum(h_0); 	%Normalize to sqrt(2);
if(TYPE=='max'),
  h_0 = fliplr(h_0);
end;
if(abs(sum(h_0 .^ 2))-1 > 1e-4) 
  error('Numerically unstable for this value of "N".');
end;
h_1 = rot90(h_0,2);
h_1(1:2:N)=-h_1(1:2:N);
