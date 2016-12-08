% pop_bemobil_dipfit() - GUI for providing parameters for dipole fitting
%
% Usage:
%   >>  EEG = pop_bemobil_dipfit( EEG );
%
% Inputs:
%   EEG     - EEGLAB EEG structure
%    
% Outputs:
%   EEG     - EEGLAB EEG structure with fitted dipole in EEG.dipfit
%
% See also: 
%   BEMOBIL_DIPFIT, EEGLAB

function [ EEG ] = pop_bemobil_dipfit( EEG )

if nargin < 1
	help pop_bemobil_dipfit;
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
    { 'style', 'text', 'string', 'Headmodel'}
    { 'style', 'edit', 'string', 'BEM' 'tag' 'headmodel'}

    { 'style', 'text', 'string', 'Channels to include'}
    { 'style', 'edit', 'string', '[1:15 17:31 33:128 156]' 'tag' 'channels_to_include'}

    { 'style', 'text', 'string', 'Components to fit'}
    { 'style', 'edit', 'string', '[]' 'tag' 'components_to_fit'}

    { 'style', 'text', 'string', '%RV threshold'}
    { 'style', 'edit', 'string', '15' 'tag' 'RV_threshold'}

    { 'style', 'text', 'string', 'Remove dipoles outside head'}
    { 'style', 'checkbox', 'value', 1, 'tag' 'remove_outside_head'}

    { 'style', 'text', 'string', 'Fit bilateral dipoles'}
    { 'style', 'checkbox', 'value', 1, 'tag' 'fit_bilateral_dipoles'}

    };

[ tmp1 tmp2 strhalt structout ] = inputgui('geom', uigom, 'uilist', uilist, ...
    'title', 'dipfit', 'helpcom', 'pophelp(''pop_bemobil_dipfit'');');

if isempty(strhalt)
    disp('Canceling dipfit...')
    return;
end;

% close menu figure after OK
close(gcbf);

% prepare values from input gui
structout.RV_threshold = str2double(structout.RV_threshold);
structout.channels_to_include = str2num(structout.channels_to_include);

if structout.remove_outside_head
    structout.remove_outside_head = 'on';
end
if structout.fit_bilateral_dipoles
    structout.fit_bilateral_dipoles = 2;
end

EEG = bemobil_dipfit( EEG, structout.headmodel,...
    structout.channels_to_include, structout.components_to_fit,...
    structout.RV_threshold, structout.remove_outside_head,...
    structout.fit_bilateral_dipoles);
end

