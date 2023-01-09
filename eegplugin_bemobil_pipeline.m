% eegplugin_bemobil_pipeline() - BeMoBIL pipeline plugin for eeglab 
% 
% The BeMoBIL Pipeline is an open-source MATLAB toolbox for fully synchronized, automatic, transparent, and replicable
% import, processing and visualization of mobile brain/body imaging and other EEG data. It includes wrappers for EEGLAB
% functions, the use of various EEGLAB plugins, and comes with additional new functionalities. All parameters are
% configurable in central scripts and everything is stored in the EEG.etc struct. Additionally, analytics plots are
% generated for each step. A comprehensive guide to installing, using, and understanding the pipeline can be found in
% our wiki on github!
%
% https://github.com/BeMoBIL/bemobil-pipeline
%
% Usage:
%   >> eegplugin_bemobil_pipeline(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Klug, M., Jeung, S., Wunderlich, A., Gehrke, L., Protzak, J., Djebbara, Z., Argubi-Wollesen, A., Wollesen, B., & Gramann, K. (2022). 
%
% See also: eeglab()
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_bemobil_pipeline(fig, try_strings, catch_strings)

    %global EEG
    vers = '2.0.0';
    if nargin < 3
        error('eegplugin_bemobil_pipeline requires 3 arguments');
    end
    
    % add bemobil pipeline folder to path
    % -----------------------
    p = fileparts(which('eegplugin_bemobil_pipeline'));
    addpath(genpath(p));
    
end