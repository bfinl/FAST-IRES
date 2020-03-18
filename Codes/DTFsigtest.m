function [gamma2] = DTFsigtest(gamma2,gamma2_sig)
% DTFsigtest - perform significance test for the calculated DTF values.
%
% Usage: sig_gamma2 = DTFsigtest(gamma2,gamma2_sig)
%
% Input: gamma2 - the calculated DTF values (nchan x nchan x nfreq)
%           gamma2_sig - the significance cut-off values (nchan x nchan x nfreq)
%
% Output: gamma2 - the significant DTF values obtained through surrogate
%              data testing
%
% Description: This function takes the calculated DTF values and the
%              distribution of the DTF surrogate dataset and calculates the
%              channels which have significant information outflow.
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

% Generate the significant cut-off
nchan = size(gamma2,1);
nfreq = size(gamma2,3);

for i=1:nchan
    for j=1:nchan
        for k=1:nfreq
           if i~=j
                if gamma2(i,j,k) < gamma2_sig(i,j,k)
                    gamma2(i,j,k) = 0;
                end
           else
               gamma2(i,j,k) = 0;
           end
        end
    end
end