function [x,N] = makesig(SigName,N)
% [x,N] = makesig(SigName,N) Creates artificial test signal identical to the
%     standard test signals proposed and used by D. Donoho and I. Johnstone
%     in WaveLab (- a matlab toolbox developed by Donoho et al. the statistics
%     department at Stanford University).
%
%    Input:  SigName - Name of the desired signal (Default 'all')
%                        'AllSig' (Returns a matrix with all the signals)
%                        'HeaviSine'
%                        'Bumps'
%                        'Blocks'
%                        'Doppler'
%                        'Ramp'
%                        'Cusp'
%                        'Sing'
%                        'HiSine'
%                        'LoSine'
%                        'LinChirp'
%                        'TwoChirp'
%                        'QuadChirp'
%                        'MishMash'
%                        'Werner Sorrows' (Heisenberg)
%                        'Leopold' (Kronecker)
%            N       - Length in samples of the desired signal (Default 512)
%
%    Output: x   - vector/matrix of test signals
%            N   - length of signal returned
%
%    See also: 
%
%    References:
%            WaveLab can be accessed at
%            www_url: http://playfair.stanford.edu/~wavelab/
%            Also see various articles by D.L. Donoho et al. at
%            web_url: http://playfair.stanford.edu/

%File Name: makesig.m
%Last Modification Date: 08/30/95	15:52:03
%Current Version: makesig.m	2.4
%File Creation Date: Thu Jun  8 10:31:11 1995
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
%Change History: This m-file is a copy of the  code provided with WaveLab
%                customized to be consistent with RWT.
%                Jan Erik Odegard <odegard@ece.rice.edu> Thu Jun  8 1995
%

if(nargin < 1)
  SigName = 'AllSig';
  N = 512;
elseif(nargin == 1)
  N = 512;
end;
t = (1:N) ./N;
x = [];
y = [];
if(strcmp(SigName,'HeaviSine') | strcmp(SigName,'AllSig')),
  y = 4.*sin(4*pi.*t);
  y = y - sign(t - .3) - sign(.72 - t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Bumps') | strcmp(SigName,'AllSig')),
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
  wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
  y = zeros(size(t));
  for j =1:length(pos)
    y = y + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
  end 
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Blocks') | strcmp(SigName,'AllSig')),
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
  y = zeros(size(t));
  for j=1:length(pos)
    y = y + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
  end
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Doppler') | strcmp(SigName,'AllSig')),
  y = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Ramp') | strcmp(SigName,'AllSig')),
  y = t - (t >= .37);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Cusp') | strcmp(SigName,'AllSig')),
  y = sqrt(abs(t - .37));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Sing') | strcmp(SigName,'AllSig')),
  k = floor(N * .37);
  y = 1 ./abs(t - (k+.5)/N);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'HiSine') | strcmp(SigName,'AllSig')),
  y = sin( pi * (N * .6902) .* t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'LoSine') | strcmp(SigName,'AllSig')),
  y = sin( pi * (N * .3333) .* t);
end;
x = [x;y];
y = [];
if(strcmp(SigName,'LinChirp') | strcmp(SigName,'AllSig')),
  y = sin(pi .* t .* ((N .* .125) .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'TwoChirp') | strcmp(SigName,'AllSig')),
  y = sin(pi .* t .* (N .* t)) + sin((pi/3) .* t .* (N .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'QuadChirp') | strcmp(SigName,'AllSig')),
  y = sin( (pi/3) .* t .* (N .* t.^2));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'MishMash') | strcmp(SigName,'AllSig')),  
  % QuadChirp + LinChirp + HiSine
  y = sin( (pi/3) .* t .* (N .* t.^2)) ;
  y = y +  sin( pi * (N * .6902) .* t);
  y = y +  sin(pi .* t .* (N .* .125 .* t));
end;
x = [x;y];
y = [];
if(strcmp(SigName,'WernerSorrows') | strcmp(SigName,'AllSig')),
  y = sin( pi .* t .* (N/2 .* t.^2)) ;
  y = y +  sin( pi * (N * .6902) .* t);
  y = y +  sin(pi .* t .* (N .* t));
  pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
  hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
  wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
  for j =1:length(pos)
    y = y + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
  end 
end;
x = [x;y];
y = [];
if(strcmp(SigName,'Leopold') | strcmp(SigName,'AllSig')),
  y = (t == floor(.37 * N)/N); 		% Kronecker
end;
x = [x;y];
y = [];

%  disp(sprintf('MakeSignal: I don*t recognize << %s>>',SigName))
%  disp('Allowable SigNames are:')
%  disp('AllSig'),
%  disp('HeaviSine'),
%  disp('Bumps'),
%  disp('Blocks'),
%  disp('Doppler'),
%  disp('Ramp'),
%  disp('Cusp'),
%  disp('Crease'),
%  disp('Sing'),
%  disp('HiSine'),
%  disp('LoSine'),
%  disp('LinChirp'),
%  disp('TwoChirp'),
%  disp('QuadChirp'),
%  disp('MishMash'),
%  disp('WernerSorrows'),
%  disp('Leopold'),
%end
