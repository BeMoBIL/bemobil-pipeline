% convert_xdf_to_set() - GUI to convert an XDF file or a folder containing
% multiple XDF files to EEGLAB .set files
%
% Usage:
%   >>  EEG = bemobil_preprocess( filename, filepath );
%
% Inputs:
%   filename     - filename of a *.xdf file
%    
% Outputs:
%   .set file    - EEGLAB compatible .set file
%
% See also: 
%   LOAD_XDF, POP_LOAD_XDF, POP_SAVESET, EEGLAB

function pop_convert_xdf_to_set( input_args )

uilist = {
    { 'style' 'text' 'string' 'Select a single .xdf file or a folder containing several .xdf files!' }, ...
    { 'style' 'pushbutton' 'string' 'Select single file' 'Callback' @load_single_file} ...
    { 'style' 'pushbutton' 'string' 'Select folder' 'Callback' @load_folder}};
uigom = { 1 1 1 };

inputgui(uigom, uilist, 'pophelp(''pop_bemobil_import_data'')', 'Load an XDF file', 'title', 'Load .xdf file(s)');

function load_single_file(source, event)

    % close menu figure
    close(gcbf);
    % ask user
    [filename, filepath] = uigetfile('*.xdf;*.xdfz', 'Choose an XDF file -- pop_loadxdf()');
    drawnow;
    if filename == 0
        disp('Exiting XDF to set conversion.');
        return;
    end;
    convert_xdf_to_set(filename, filepath);
end

function load_folder(source, event)

    % close menu figure
    close(gcbf);
    % prompt user
    directoryname = uigetdir;
    drawnow;
    if directoryname == 0
        disp('Exiting XDF to set conversion.');
        return;
    end;
    convert_xdf_to_set([], directoryname);
end
end