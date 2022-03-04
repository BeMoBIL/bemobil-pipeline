% bemobil_dipfit() - Prepares data for dipole fitting and runs the dipole fitting procedure
%
% Usage:
%   >> [ALLEEG, EEG, CURRENTSET] = bemobil_dipfit( EEG , ALLEEG, CURRENTSET, warping_channel_names, RV_threshold,...
%     remove_outside_head, number_of_dipoles, out_filename, out_filepath)

% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   warping_channel_names   - cell array of channel number and corresponding standard 5/10 system
%       electrode names (e.g. {33,'C5'; 64,'C6'}). Electrodes should cover the head sphere, but not
%       too many ( e.g. each one frontal, occipital, temporal left/right and central).
%   RV_threshold            - number percentage of residual variance accepted, default is '15'
%   remove_outside_head     - 'on' or 'off' to remove dipoles located outside the head
%   number_of_dipoles       - '1' or '2', 2 meaning bilateral dipole fitting
%   out_filename            - output filename (OPTIONAL ARGUMENT)
%   out_filepath            - output filepath (OPTIONAL ARGUMENT - File will only be saved on disk
%       if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   Currentset              - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%   EEGLAB, coregister, pop_dipfit_settings, pop_multifit
% 
% Author: Marius Klug, 2021

function [ALLEEG, EEG, CURRENTSET] = bemobil_dipfit( EEG , ALLEEG, CURRENTSET, warping_channel_names, RV_threshold,...
    remove_outside_head, number_of_dipoles, out_filename, out_filepath)

% only save a file on disk if both a name and a path are provided
save_file_on_disk = (exist('out_filename', 'var') && exist('out_filepath', 'var'));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(out_filepath);
    if ismember(out_filename, {dir_files.name})
        warning([out_filename ' file already exists in: ' out_filepath '. File will be overwritten...']);
    end
end

% load some standard data for dipfit
dipfitdefs;

% if strcmp(headmodel, 'evy') % add later when available
% elseif strcmp(headmodel, 'klaus') % add later when available
% end

% use standard BEM headmodel
% change relevant EEG electrode labels so that it matches standard 5/10 template

% build the command for editing channel labels with variable channel editing counts
if ~isempty(warping_channel_names)
    
    if isnumeric(warping_channel_names) && all(size(warping_channel_names) == [1 9])
        % transform is provided directly
        transform = warping_channel_names;
    else
        command = 'EEG = pop_chanedit(EEG';
        for channel = 1:length(warping_channel_names)
            command = [command ', ''changefield'',{' num2str(warping_channel_names{channel,1}) ' ''labels'' ''' warping_channel_names{channel,2} '''}'];
        end
        command = [command ');'];

        % run it
        eval(command);


        % warp the locations to the standard head model
        [newlocs transform] = coregister(EEG.chanlocs, template_models(2).chanfile, 'warp', 'auto', 'manual', 'off');
    end
else
    transform = [0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485];
end

% do the dipole fitting (default headmodel)
dipfit_base_path = fileparts(which('dipfitdefs'));
EEG = pop_dipfit_settings( EEG, 'hdmfile',fullfile(dipfit_base_path, 'standard_BEM', 'standard_vol.mat'),...
    'coordformat','MNI',...
    'mrifile', fullfile(dipfit_base_path, 'standard_BEM', 'standard_mri.mat'),...
    'chanfile', fullfile(dipfit_base_path, 'standard_BEM', 'elec', 'standard_1005.elc'),...
    'coord_transform',transform ,...
    'chansel',1:EEG.nbchan);


EEG = pop_multifit(EEG, 1:size(EEG.icaweights,1) ,'threshold',RV_threshold,...
    'dipoles', number_of_dipoles, 'rmout', remove_outside_head);

% save dipfit info in EEG.etc

EEG.etc.dipfit.warping_channel_names = warping_channel_names;
EEG.etc.dipfit.RV_threshold = RV_threshold;
EEG.etc.dipfit.remove_outside_head = remove_outside_head;
EEG.etc.dipfit.number_of_dipoles = number_of_dipoles;

% new data set in EEGLAB
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',out_filename,'filepath', out_filepath);
    disp('...done');
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);