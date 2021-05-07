% pop_bemobil_preprocess() - This pipeline is based on Makotos preprocessing pipeline:
%   https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline with further
%   discussion and changes. Specifically it is adapted to data recorded at
%   BeMoBIl.
%
% Usage:
%   >>  [ EEG ] = pop_bemobil_preprocess( EEG )
%
% Inputs:
%   EEG                   - eeglab EEG struct
%    
% Outputs:
%   EEG     - Preprocessed EEGLAB EEG structure
%
% See also: 
%   BEMOBIL_PREPROCESS, EEGLAB

function [ ALLEEG EEG CURRENTSET ] = pop_bemobil_preprocess( ALLEEG, EEG, CURRENTSET )

if nargin < 1
	help pop_bemobil_preprocess;
	return;
end;	

% add preprocess GUI
uigeom_cols = 2;
uigeom_rows = 6;

uigom = {
    {uigeom_cols uigeom_rows [0 0] [1 1]} 
    {uigeom_cols uigeom_rows [1 0] [1 1]}

    {uigeom_cols uigeom_rows [0 1] [1 1]} 
    {uigeom_cols uigeom_rows [1 1] [1 1]}

    {uigeom_cols uigeom_rows [0 2] [1 1]} 
    {uigeom_cols uigeom_rows [1 2] [1 1]}

    {uigeom_cols uigeom_rows [0 3] [1 1]} 
    {uigeom_cols uigeom_rows [1 3] [1 1]}

    {uigeom_cols uigeom_rows [0 4] [1 1]} 
    {uigeom_cols uigeom_rows [1 4] [1 1]}

    {uigeom_cols uigeom_rows [0 5] [1 1]} 
    {uigeom_cols uigeom_rows [1 5] [1 1]}

    };

uilist = {
    { 'style', 'text', 'string', 'Channel locations filename'}
    { 'style', 'edit', 'string', 'channel_locations.elc' 'tag' 'chanlocs'}

    { 'style', 'text', 'string', 'Channel(s) to remove'}
    { 'style', 'edit', 'string', 'N29 N30 N31' 'tag' 'remove_chans'}

    { 'style', 'text', 'string', 'EOG channels'}
    { 'style', 'edit', 'string', 'G16 G32' 'tag' 'heog_chans'}

    { 'style', 'text', 'string', 'Resample frequency'}
    { 'style', 'edit', 'string', '250' 'tag' 'resample_freq'}

    { 'style', 'text', 'string', 'Filter(FIR) low cutoff frequency'}
    { 'style', 'edit', 'string', '1' 'tag' 'locutoff'}

    { 'style', 'text', 'string', 'Filter(FIR) high cutoff frequency'}
    { 'style', 'edit', 'string', '124' 'tag' 'highcutoff'}

    };

[ tmp1 tmp2 strhalt structout ] = inputgui('geom', uigom, 'uilist', uilist, ...
    'title', 'Preprocessing', 'helpcom', 'pophelp(''pop_bemobil_preprocess'');');

if isempty(strhalt)
    disp('Canceling preprocessing...')
    return;
end;

% close menu figure after OK
close(gcbf);

% prepare values from input gui
structout.remove_chans = strsplit(structout.remove_chans, ' ')';
structout.heog_chans = strsplit(structout.heog_chans, ' ')';
structout.resample_freq = str2double(structout.resample_freq);
structout.locutoff = str2double(structout.locutoff);
structout.highcutoff = str2double(structout.highcutoff);

% run processing function with values from gui input
[ ALLEEG EEG CURRENTSET ] = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, structout.chanlocs, structout.remove_chans,...
    structout.heog_chans, structout.locutoff, structout.highcutoff, ...
    structout.resample_freq);
end

