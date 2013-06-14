function option = setopt(opt_par,default);
%    option = setopt(opt_par,default); 
%
%    SETOPT can modify a default option vector with user specified options.
%
%    Input: 
%       opt_par : Users desired option vector
%       default : Program default option vector
%
%    Output: 
%       option : New option vector
%
%    Example:
%       opt_par = [1 2 3 4];
%       default = [1 1 1 1];
%       option = setopt(opt_par,default)
%       option = 1     2     3     4
%

%File Name: setopt.m
%Last Modification Date: 10/26/95	13:06:18
%Current Version: setopt.m	2.4
%File Creation Date: Thu Nov 11 10:40:22 1993
%Author: Jan Erik Odegard  <odegard@ece.rice.edu>
%
%Copyright (c) 2000 RICE UNIVERSITY. All rights reserved.
%Created by Jan Erik Odegard, Department of ECE, Rice University. 
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

if (nargin < 2) 
  error('You need two inputs');
end;
len = length(opt_par);
option = zeros(size(default));
option(1:len) = opt_par(1:len);
option = option + (option == 0).*default; 
