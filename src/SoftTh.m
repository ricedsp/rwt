function  x = SoftTh(y,thld)
%    x = SoftTh(y,thld); 
%
%    SOFTTH soft thresholds the input signal y with the threshold value
%    thld.
%
%    Input:  
%       y    : 1D or 2D signal to be thresholded
%       thld : Threshold value
%
%    Output: 
%       x : Soft thresholded output (x = sign(y)(|y|-thld)_+)
%
%  HERE'S AN EASY WAY TO RUN THE EXAMPLES:
%  Cut-and-paste the example you want to run to a new file 
%  called ex.m, for example. Delete out the % at the beginning 
%  of each line in ex.m (Can use search-and-replace in your editor
%  to replace it with a space). Type 'ex' in matlab and hit return.
%
%
%    Example:
%       y = makesig('Doppler',8);
%       thld = 0.2;
%       x = SoftTh(y,thld)
%       x = 0 0 0 -0.0703 0 0.2001 0.0483 0 
%
%    See also: HardTh
%
%    Reference: 
%       "De-noising via Soft-Thresholding" Tech. Rept. Statistics,
%       Stanford, 1992. D.L. Donoho.
%

%File Name: SoftTh.m
%Last Modification Date: 8/15/95	17:49:48
%Current Version: SoftTh.m	2.4
%File Creation Date: Mon Mar  7 10:38:45 1994
%Author: Haitao Guo  <harry@jazz.rice.edu>
%
%Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
%Created by Haitao Guo, Department of ECE, Rice University. 
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

x = abs(y);
x = sign(y).*(x >= thld).*(x - thld); 
