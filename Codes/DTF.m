function [gamma2] = DTF(ts,low_freq,high_freq,p,fs)
% DTF - perform Directed Transfer Function analysis among multi-channel time series. 
%
% Usage: gamma2 = DTF(ts,low_freq,high_freq,p,fs)
%
% Input: ts - the time series where each column is the temporal data from a
%           single channel
%           low_freq - the lowest frequency bound to compute the DTF
%           high_freq - the highest frequency bound to compute the DTF
%           p - the order of the MVAR model
%           fs - the sampling frequency 
%
% Output: gamma2 - the computed DTF values
%
% Description: This program calculates the DTF values for a given frequency
%              range for the input time series. The output is in the form
%              gamma2(a,b,c) where a = the sink channel, b = the source
%              channel, c = the frequency index.
%
% ARfit Package:
% The ARfit package is used in DTF computation. 
% See below for detailed description of the ARfit package:
% A. Neumaier and T. Schneider, 2001: Estimation of parameters and eigenmodes of 
% multivariate autoregressive models. ACM Trans. Math. Softw., 27, 27?57.
% T. Schneider and A. Neumaier, 2001: Algorithm 808: ARfit-A Matlab package for the 
% estimation of parameters and eigenmodes of multivariate autoregressive models. 
% ACM Trans. Math. Softw., 27, 58?65. 
% http://www.gps.caltech.edu/~tapio/arfit/' 
%
% Program Authors: Lei Ding and Christopher Wilke, University of Minnesota, USA
%
% User feedback welcome: e-mail: econnect@umn.edu
%

% License
% ==============================================================
% This program is part of the eConnectome.
% 
% Copyright (C) 2010 Regents of the University of Minnesota. All rights reserved.
% Correspondence: binhe@umn.edu
% Web: econnectome.umn.edu
%
% This program is free software for academic research: you can redistribute it and/or modify
% it for non-commercial uses, under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see http://www.gnu.org/copyleft/gpl.html.
%
% This program is for research purposes only. This program
% CAN NOT be used for commercial purposes. This program 
% SHOULD NOT be used for medical purposes. The authors 
% WILL NOT be responsible for using the program in medical
% conditions.
% ==========================================

% Revision Logs
% ==========================================
%
% Yakang Dai, 01-Mar-2010 15:20:30
% Release Version 1.0 beta 
%
% ==========================================

% Default sampling rate is 400 Hz
if nargin < 5
    fs = 400;
end

% The number of frequencies to compute the DTF over
tot_range = [low_freq:high_freq];
nfre = length(tot_range);

% The number of channels in the time series
nchan = size(ts,2);

% The sampling period
dt = 1/fs;

% Create the MVAR matrix for the time series
[w,A] = arfit(ts,p,p);

% Rearrange the format of the MVAR matrix
B = [];
B(:,:,1) = -eye(nchan);
for i=1:nchan
    for j=1:nchan
%         dc(i,j) = sum(A(i,j:nchan:nchan*p).*A(i,j:nchan:nchan*p));
        B(i,j,2:p+1) = A(i,j:nchan:nchan*p);
    end
end

% Calculate the non-normalized DTF value
theta2 = [];
for k = 1:nfre
    Af = zeros(nchan,nchan);
    fre = tot_range(k);
    for i = 1:nchan
        for j = 1:nchan
            for h = 1:p+1
                Af(i,j) = Af(i,j)-B(i,j,h)*exp(-pi*fre*dt*(h-1)*2i);
            end
        end
    end
    dett2 = det(Af);
    dett2 = dett2.*conj(dett2);
    for i = 1:nchan
        for j = 1:nchan
            Apf = Af;
            Apf(:,i) = [];
            Apf(j,:) = [];
            det2 = det(Apf);
            det2 = det2.*conj(det2);
            theta2(i,j,k) = det2/dett2;
        end
    end
end

% Calculate the normalized DTF values
gamma2 = [];
for k=1:nfre
    for i=1:nchan
        for j=1:nchan
            gamma2(i,j,k) = theta2(i,j,k) / sum(theta2(i,:,k),2);
        end
    end
end
