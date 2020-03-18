function [gamma2] = DTFvalue(A,low_freq,high_freq,fs)
% DTFvalue - calculate DTF values given a coefficient matrix from a MVAAR model
%
% Usage: gamma2 = DTFvalue(A,low_freq,high_freq,fs)
%
% Authors: Lei Ding and Christopher Wilke
% 
% Input: A - MVAR model parameters
%        low_freq - lowest frequency to perform DTF analysis
%        high_freq - highest frequency to perform DTF analysis
%        fs - the sampling frequency 
%
% Output: gamma2 - the DTF values 
%
% Description: This function performs DTF analysis given the parameters 
%              from a MVAAR model. The output is in the form gamma2(i,j,k) 
%              where i is the sink channel, j is the source channel and k
%              is the frequency index.
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
% Yakang Dai, 20-June-2010 00:32:30
% Release Version 1.0  
%
% ==========================================

% Number of channels
nchan = size(A,1);

% Sampling frequency
% Fs = 500;
Fs=fs;
% Sampling period
dt = 1/Fs;

% Total range over which to perform DTF
tot_range = [low_freq:high_freq];
% Number of frequencies
nfre = length(tot_range);

% Order of the mvaar model
p = size(A,3) - 1;

% Non-normalized DTF
theta2 = [];
for k = 1:nfre
    Af = zeros(nchan,nchan);
    fre = tot_range(k);
    for i = 1:nchan
        for j = 1:nchan
            for h = 1:p+1
                Af(i,j) = Af(i,j) - A(i,j,h)*exp(-pi*fre*dt*(h-1)*2i);
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

% Normalized DTF
gamma2 = [];
for k = 1:nfre
    for i = 1:nchan
        for j = 1:nchan
            gamma2(i,j,k) = theta2(i,j,k) / sum(theta2(i,:,k),2);
        end
    end
end
