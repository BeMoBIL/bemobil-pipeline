% eegplugin_bemobil_pipeline() - BeMoBIL pipeline plugin for eeglab
%   BeMoBIL pipeline is a collection of functionalities
%   for processing EEG/MoBI data collected at BeMoBIL.
%
% Usage:
%   >> eegplugin_bemobil_pipeline(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Lukas Gehrke, TU Berlin
%
% See also: eeglab()

% Copyright (C) 2016, Lukas Gehrke
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

function eegplugin_bemobil_pipeline(fig, try_strings, catch_strings)

p = fileparts(which('run_bemobil_pipeline'));
if isempty(p)
    p = fileparts(which('eeglab'));
    p = [p filesep 'plugins' filesep 'bemobil_pipeline'];
end
addpath(p);

highlevelmenu = findobj(fig,'Label','Tools');

menu = uimenu(highlevelmenu, 'label', 'BeMoBIL Pipeline', 'separator', 'on', 'userdata', 'startup:on');

load_data_menu = uimenu(menu, 'label', '1. Load data');
uimenu( load_data_menu, 'label', 'From MoBILab', 'CallBack', 'disp(''Starting MoBILab...''); runmobilab;');
uimenu( load_data_menu, 'label', 'Convert .xdf to .set', 'CallBack', 'disp(''Starting data conversion plugin...''); pop_convert_xdf_to_set();');
uimenu( load_data_menu, 'label', 'From .set file(s)', 'CallBack', 'disp(''Loading data for eeglab...''); ALLEEG = pop_loadset(); eeglab redraw; disp(''Done.'')');
uimenu( load_data_menu, 'label', 'Merge multiple .set files and save the created file (optional step)', 'CallBack', 'disp(''Starting merging datasets...''); [ ALLEEG EEG CURRENTSET ] = bemobil_merge( ALLEEG, EEG, CURRENTSET ); eeglab redraw; disp(''Done.'')');

preprocess_menu = uimenu(menu, 'label', '2. Preprocessing', 'CallBack', 'disp(''Starting preprocessing...''); [ ALLEEG EEG CURRENTSET ] = pop_bemobil_preprocess(ALLEEG, EEG, CURRENTSET); eeglab redraw;');

channel_data_cleaning_menu = uimenu(menu, 'label', '3. Data cleaning (channel level)');
uimenu( channel_data_cleaning_menu, 'label', '1. Reject irrelevant experiment segments', 'CallBack', 'disp(''Starting segments GUI...''); EEG = pop_bemobil_segment(EEG); eeglab redraw;');
uimenu( channel_data_cleaning_menu, 'label', 'CleanLine (optional)', 'CallBack', 'disp(''Starting CleanLine...''); EEG = pop_cleanline(EEG); eeglab redraw;');
uimenu( channel_data_cleaning_menu, 'label', 'Clean_Rawdata (ASR) (optional)', 'CallBack', 'disp(''Starting Clean_Rawdata...''); EEG = pop_clean_rawdata(EEG); eeglab redraw;');
uimenu( channel_data_cleaning_menu, 'label', '2. Manual channel rejection', 'CallBack', 'disp(''Select channels...''); EEG = pop_select(EEG); eeglab redraw;');
uimenu( channel_data_cleaning_menu, 'label', '3. Manual time domain cleaning', 'CallBack', 'disp(''Starting data cleaning on channel level...''); pop_eegplot(EEG); eeglab redraw;');

signal_decomposition_menu = uimenu( menu, 'label', '4. Signal decomposition', 'CallBack', 'disp(''Opening AMICA GUI...''); [ALLEEG EEG CURRENTSET] = pop_bemobil_signal_decomposition(ALLEEG, EEG, CURRENTSET); eeglab redraw;');

comp_data_cleaning_menu = uimenu( menu, 'label', '5. Data cleaning (component level)', 'CallBack', 'disp(''Starting data cleaning on component level...''); pop_eegplot(EEG, 0);');

dipfit_menu = uimenu( menu, 'label', '6. Fit single/dual dipoles');
uimenu( dipfit_menu, 'label', '1. Interpolation & Av Ref', 'CallBack', 'disp(''Starting interpolation and average referencing...''); EEG = bemobil_interp_avref(EEG); eeglab redraw;');
uimenu( dipfit_menu, 'label', '2. Dipfit', 'CallBack', 'disp(''Starting dipfit GUI...''); EEG = pop_bemobil_dipfit(EEG); eeglab redraw;');

finalize_single_subject_menu = uimenu( menu, 'label', '7. Finalize single subject processing', 'CallBack', 'disp(''Started GUI to copy weights...''); EEG = pop_editset(EEG); eeglab redraw;');

batch_script_menu = uimenu( menu, 'label', '8. Batch processing', 'CallBack', 'disp(''Opening batch template script...''); edit sample_bemobil_batch.m;');
%uimenu( finalize_single_subject_menu, 'label', '2. Epoching');

end