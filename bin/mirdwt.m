function [x,L] = mirdwt(yl,yh,h,L);
%    function [x,L] = mirdwt(yl,yh,h,L);
% 
%    Function computes the inverse redundant discrete wavelet
%    transform x  for a 1D or 2D input signal. (Redundant means here
%    that the sub-sampling after each stage of the forward transform
%    has been omitted.) yl contains the lowpass and yl the highpass
%    components as computed, e.g., by mrdwt. In the case of a 2D
%    signal, the ordering in
%    yh is [lh hl hh lh hl ... ] (first letter refers to row, second
%    to column filtering).  
%
%    Input:
%       yl : lowpass component
%       yh : highpass components
%       h  : scaling filter
%       L  : number of levels. In the case of a 1D signal, 
%            length(yl) must  be divisible by 2^L;
%            in the case of a 2D signal, the row and
%            the column dimension must be divisible by 2^L.
%   
%    Output:
%	     x : finite length 1D or 2D signal
%	     L : number of levels
%
%  HERE'S AN EASY WAY TO RUN THE EXAMPLES:
%  Cut-and-paste the example you want to run to a new file 
%  called ex.m, for example. Delete out the % at the beginning 
%  of each line in ex.m (Can use search-and-replace in your editor
%  to replace it with a space). Type 'ex' in matlab and hit return.
%
%
%    Example 1:
%    xin = makesig('Leopold',8);
%    h = daubcqf(4,'min');
%    L = 1;
%    [yl,yh,L] = mrdwt(xin,h,L);
%    [x,L] = mirdwt(yl,yh,h,L)
%    x = 0.0000 1.0000 0.0000 -0.0000 0 0 0 -0.0000
%    L = 1
%  
%    Example 2:  
%    load lena;
%    h = daubcqf(4,'min');
%    L = 2;
%    [ll_lev2,yh,L] = mrdwt(lena,h,L); % lena is a 256x256 matrix
%    N = 256;
%    lh_lev1 = yh(:,1:N); 
%    hl_lev1 = yh(:,N+1:2*N); 
%    hh_lev1 = yh(:,2*N+1:3*N);
%    lh_lev2 = yh(:,3*N+1:4*N); 
%    hl_lev2 = yh(:,4*N+1:5*N); 
%    hh_lev2 = yh(:,5*N+1:6*N);
%    figure; colormap(gray); imagesc(lena); title('Original Image');
%    figure; colormap(gray); imagesc(ll_lev2); title('LL Level 2');
%    figure; colormap(gray); imagesc(hh_lev2); title('HH Level 2');
%    figure; colormap(gray); imagesc(hl_lev2); title('HL Level 2');
%    figure; colormap(gray); imagesc(lh_lev2); title('LH Level 2');
%    figure; colormap(gray); imagesc(hh_lev1); title('HH Level 1');
%    figure; colormap(gray); imagesc(hl_lev2); title('HL Level 1');
%    figure; colormap(gray); imagesc(lh_lev2); title('LH Level 1');
%    [lena_Hat,L] = mirdwt(ll_lev2,yh,h,L);
%    figure; colormap(gray); imagesc(lena_Hat); 
%                            title('Reconstructed Image');
%
%    See also: mdwt, midwt, mrdwt
%
%    Warning! min(size(yl))/2^L should be greater than length(h)
%

%File Name: mirdwt.m
%Last Modification Date: 08/07/95	15:14:21
%Current Version: mirdwt.m	2.4
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
%Mon Aug  7 15:09:51 CDT 1995
%Rebecca Hindman <hindman@ece.rice.edu>
%Added L to function line so that it can be displayed as an output
% 
%Modification #2  
%Thursday Mar 2 2000
% Added Example 2
% Felix Fernandes <felixf@rice.edu>
















