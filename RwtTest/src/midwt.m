function [y,L] = midwt(x,h,L);
%    [x,L] = midwt(y,h,L);
% 
%    Function computes the inverse discrete wavelet transform x for a 1D or
%    2D input signal y using the scaling filter h.
%
%    Input:
%	y : finite length 1D or 2D input signal (implicitly periodized)
%           (see function mdwt to find the structure of y)
%       h : scaling filter
%       L : number of levels. In the case of a 1D signal, length(x) must be
%           divisible by 2^L; in the case of a 2D signal, the row and the
%           column dimension must be divisible by 2^L.  If no argument is
%           specified, a full inverse DWT is returned for maximal possible
%           L.
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

%
%
%File Name: midwt.m
%Last Modification Date: 08/07/95	15:13:52
%Current Version: midwt.m	2.4
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
%Mon Aug  7 11:52:33 CDT 1995
%Rebecca Hindman <hindman@ece.rice.edu>
%Added L to function line so that it can be displayed as an output
% 
%Thu Mar  2 13:07:11 CDT 2000
%Ramesh Neelamani<neelsh@ece.rice.edu>
%Revamped the help file
%




