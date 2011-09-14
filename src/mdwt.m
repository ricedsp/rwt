function [y,L] = mdwt(x,h,L);
%    [y,L] = mdwt(x,h,L);
%
%    Function computes the discrete wavelet transform y for a 1D or 2D input
%    signal x using the scaling filter h.
%
%    Input:
%	x : finite length 1D or 2D signal (implicitly periodized)
%       h : scaling filter
%       L : number of levels. In the case of a 1D signal, length(x) must be
%           divisible by 2^L; in the case of a 2D signal, the row and the
%           column dimension must be divisible by 2^L. If no argument is
%           specified, a full DWT is returned for maximal possible L.
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
%
%
%
%    See also: midwt, mrdwt, mirdwt
%

%File Name: mdwt.m
%Last Modification Date: 08/07/95	15:13:25
%Current Version: mdwt.m	2.4
%File Creation Date: Wed Oct 19 10:51:58 1994
%Author: Markus Lang  <lang@jazz.rice.edu>
%
%Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
%Created by Markus Lang, Department of ECE, Rice University. 
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
%
%Change History:
% 
%Modification #1
%Mon Aug  7 11:42:11 CDT 1995
%Rebecca Hindman <hindman@ece.rice.edu>
%Added L to function line so that it can be displayed as an output
%
%Change History:
% 
%Modification #1
%Thu Mar  2 13:07:11 CDT 2000
%Ramesh Neelamani<neelsh@ece.rice.edu>
%Revamped the help file
%



