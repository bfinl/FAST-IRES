function [new_gamma2] = DTFsigvalues(ts, low_freq, high_freq, p, fs, shufftimes, siglevel, handle)
% DTFsigvalues - compute statistical significance values for relative DTF values
%
% Input: ts - the time series where each column is the temporal data from a single trial.
%           low_freq - the lowest frequency to perform DTF analysis.
%           high_freq - the highest frequency to perform DTF analysis.
%           p - the order of the model.
%           fs - the sampling rate of the data.
%           siglevel - the significance level, default is 0.05.
%           shufftimes - the shuffling times, default is 1000.
%           handle - the handle of the uicontrol for displaying computation
%                         progress, the progress will be displayed in the command window 
%                         of MATLAB if handle = []. 
%
% Output: new_gamma2 - the significant points from the surrogate DTF analysis
%
% Description: This function generates surrogate datasets by phase shifting
%              the original time series and then performs DTF analysis on
%              these new time series for statistical testing. The output is
%              in the form gamma2_sig(a,b,c) where a = the sink channel, 
%              b = the source channel, and c = the frequency index.
%
% Program Author: Christopher Wilke, University of Minnesota, USA
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

% Number of shuffled datasets to create
if isempty(shufftimes) 
    nreps = 1000;
else
    nreps = shufftimes;
end

if isempty(siglevel)
    tvalue = 0.05;
else
    tvalue = siglevel;
end

% Number of frequencies
tot_range = [low_freq:1:high_freq];
nfreq = length(tot_range);

% Number of channels in the time series
nchan = size(ts,2);

sig_size = floor(tvalue * nreps)+1;
new_gamma = zeros(sig_size-1,nchan,nchan,nfreq);

% Number of surrogate datasets to generate
for i=1:nreps

    % Display the progress of the function on a uicontrol
    rate = round(100 * i / nreps);
    progress = ['Completing ' num2str(rate) '%'];
    if ~isempty(handle)
        set(handle, 'string',progress);
        drawnow;
    else
        fprintf('%s \n',['Completing ' num2str(rate) '%']);
    end
    
    % Generate a surrogate time series
    for j=1:nchan
        Y = fft(ts(:,j));
        Pyy = sqrt(Y.*conj(Y));
        Phyy = Y./Pyy;
        index = 1:size(ts,1);
        index = surrogate(index);
        Y = Pyy.*Phyy(index);
        newts(:,j) = real(ifft(Y));
    end
    
    % Compute the DTF value for each surrogate time series
    gamma2 = DTF(newts,low_freq,high_freq,p,fs);
     
    % Save the surrogate DTF values
    new_gamma(sig_size,:,:,:) = gamma2;
    new_gamma = sort(new_gamma,'descend');
    new_gamma(sig_size,:,:,:) = [];

end

% take the surrogated DTF values at a certain signficance
new_gamma2 = zeros(nchan,nchan,nfreq);
for i = 1:nchan
    for j = 1:nchan
        for k = 1:nfreq
            new_gamma2(i,j,k) = new_gamma(sig_size-1,i,j,k);
        end
    end
end

if ~isempty(handle)
    set(handle, 'string', 'Done');
    drawnow;
else
    fprintf('%s \n', 'Done');
end