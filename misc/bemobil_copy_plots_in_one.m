% bemobil_copy_plots_in_one() - Copy all png plots from the individual subject folders into one top-level folder as well
% as one semi-top-level folder for the processing steps for easier inspection.
%
% Usage:
%   >> bemobil_copy_plots_in_one(bemobil_config, newfoldername)

% Inputs:
%   bemobil_config          - configuration struct with all necessary information. See EEG_processing_example file
%                               that comes with the pipeline!
%   newfoldername           - OPTIONAL: name of the new folder. Default = 'allPlots'
%
% Outputs:
%   NONE, just copied files on disk
%
% Author: Marius Klug, 2021

function bemobil_copy_plots_in_one(bemobil_config, newfoldername)

if ~exist('bemobil_config','var') || isempty(bemobil_config)
    error('Bemobil config input required!')
end

if ~exist('newfoldername','var') || isempty(newfoldername)
    newfoldername = 'allPlots';
end

%% copy preprocessing
[~,message] = mkdir(fullfile(bemobil_config.study_folder,newfoldername));

if strcmp(message,'Directory already exists.')
    disp('Removing old folder...')
    rmdir(fullfile(bemobil_config.study_folder,newfoldername),'s')
    mkdir(fullfile(bemobil_config.study_folder,newfoldername))
end

disp('Copying preprocessing plots...')

basefilepath = fullfile(bemobil_config.study_folder,bemobil_config.EEG_preprocessing_data_folder);
allfolders = dir(basefilepath);
copyfilepath = fullfile(basefilepath,newfoldername);
[~,message] = mkdir(copyfilepath);
if strcmp(message,'Directory already exists.')
    disp('Removing old folder...')
    rmdir(copyfilepath,'s')
    mkdir(copyfilepath)
end

for i_folder = 3:length(allfolders)
    disp(allfolders(i_folder).name)
    allsubfiles = dir(fullfile(basefilepath,allfolders(i_folder).name));
    for i_file = 3:length(allsubfiles)
        if contains(allsubfiles(i_file).name,'png')
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(copyfilepath,allsubfiles(i_file).name))
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(bemobil_config.study_folder,newfoldername,allsubfiles(i_file).name))
        end
    end
end

disp('Copying preprocessing plots done!')

%% copy amica

disp('Copying AMICA plots...')

basefilepath = fullfile(bemobil_config.study_folder,bemobil_config.spatial_filters_folder,bemobil_config.spatial_filters_folder_AMICA);
allfolders = dir(basefilepath);
copyfilepath = fullfile(basefilepath,newfoldername);
[~,message] = mkdir(copyfilepath);
if strcmp(message,'Directory already exists.')
    disp('Removing old folder...')
    rmdir(copyfilepath,'s')
    mkdir(copyfilepath)
end

for i_folder = 3:length(allfolders)
    disp(allfolders(i_folder).name)
    allsubfiles = dir(fullfile(basefilepath,allfolders(i_folder).name));
    for i_file = 3:length(allsubfiles)
        if contains(allsubfiles(i_file).name,'png')
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(copyfilepath,allsubfiles(i_file).name))
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(bemobil_config.study_folder,newfoldername,allsubfiles(i_file).name))
        end
    end
end

disp('Copying AMICA plots done!')

%% copy final single subject plots

disp('Copying final single subject plots...')

basefilepath = fullfile(bemobil_config.study_folder,bemobil_config.single_subject_analysis_folder);
allfolders = dir(basefilepath);
copyfilepath = fullfile(basefilepath,newfoldername);
[~,message] = mkdir(copyfilepath);
if strcmp(message,'Directory already exists.')
    disp('Removing old folder...')
    rmdir(copyfilepath,'s')
    mkdir(copyfilepath)
end

for i_folder = 3:length(allfolders)
    disp(allfolders(i_folder).name)
    allsubfiles = dir(fullfile(basefilepath,allfolders(i_folder).name));
    for i_file = 3:length(allsubfiles)
        if contains(allsubfiles(i_file).name,'png')
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(copyfilepath,allsubfiles(i_file).name))
            copyfile(fullfile(allsubfiles(i_file).folder,allsubfiles(i_file).name),fullfile(bemobil_config.study_folder,newfoldername,allsubfiles(i_file).name))
        end
    end
end

disp('Copying final single subject plots done!')

disp('All copying done!')