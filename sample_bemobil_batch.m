%% README

% This is a template batch script. If you want to adapt it, which is likely,
% copy the content of this script to a new script and save it in one of your workplace folders.
% DO NOT OVERWRITE THIS SCRIPT.

close all
clear all
%% Filename informations

study_folder = 'P:\Lukas_Gehrke\studies\Spot_Rotation\';
subjects_folder = 'data\level_0\';
number_of_subjects = 11;
eeg_data_filenames ={'control_body.set' 'control_joy.set' 'test_body.set' 'test_joy.set'};

%% Parameters

load_from_xdf = false;

% processing
channel_locations_filename = 'channel_locations.elc';
channels_to_remove = ['N29' 'N30' 'N31'];
eog_channels  = ['G16' 'G32'];
locutoff = 1;
highcutoff = 124;
resample_freq = 250;
%EEG = bemobil_preprocess(EEG, channel_locations, channels_to_remove,...
%eog_channels, locutoff, highcutoff, resample_freq);

% segmenting
keep_or_remove = [''];
start_segment = ['start1'];
end_segment = ['end1'];
%EEG = bemobil_segment(EEG, keep_or_remove, start_segment, end_segment)

% CleanLine
line_frequencies = [50 100];
%EEG = pop_cleanline(EEG, 'LineFrequencies', line_frequencies);

% ASR
arg_flatline = -1;
arg_highpass = -1;
arg_channel = -1;
arg_noisy = -1;
arg_burst = -1;
arg_window = -1;
%EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window)

% manual channel rejection
%EEG = pop_select(EEG);

% manual time domain cleaning: channels
%EEG = pop_eegplot(EEG);

% AMICA
iteration = 1;
amica = 1;
num_models = 1;
max_threads = 4;
other_algorithm = '';
%EEG = bemobil_signal_decomposition(EEG, iteration, amica, num_models,
%max_threads, other_algorithm);

% manual time domain cleaning: components
%EEG = pop_eegplot(EEG, 0);

% interpolation and average referencing
%EEG = bemobil_interp_avref(EEG);

% warp electrodemontage and run dipfit
headmodel = 'BEM';
channels_to_include = [1:15 17:31 33:128 157];
components_to_fit = [];
RV_threshold = 15;
remove_outside_head = 'on';
fit_bilateral_dipoles = 2;
%EEG = bemobil_dipfit(EEG, headmodel,...
%channels_to_include, components_to_fit, RV_threshold,...
%remove_outside_head, fit_bilateral_dipoles);

% finalize single subject (copy back ICA weights to original dataset)

%% Main Loop

current_dir = pwd;
datapath = strcat(study_folder, subjects_folder);
cd(datapath);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for subject = 1:number_of_subjects
    cd(num2str(subject))
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    
    if load_from_xdf
        convert_xdf_to_set([],pwd);
    end
    
    EEG = pop_loadset('filename', eeg_data_filenames, 'filepath', pwd);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);
    
    bemobil_merge(ALLEEG,EEG,CURRENTSET,[1:length(ALLEEG)]);
    
    % preprocessing
    EEG = bemobil_preprocess(ALLEEG, EEG, CURRENTSET, channel_locations_filename, channels_to_remove,...
        eog_channels, locutoff, highcutoff, resample_freq);
    
    % remove unwanted/irrelevant segments
    %     EEG = bemobil_segment(EEG, keep_or_remove, start_segment, end_segment);
    %
    %     % optional automatic cleaning
    %     %EEG = pop_cleanline(EEG, 'LineFrequencies', line_frequencies);
    %     %EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window)
    %
    %     % manual channel rejection
    %     EEG = pop_select(EEG);
    %
    %     % manual time domain cleaning: components
    %     EEG = pop_eegplot(EEG, 0);
    %
    %     % interpolation and average referencing
    %     EEG = bemobil_interp_avref(EEG);
    %
    %     % dipfit
    %     EEG = bemobil_dipfit(EEG, headmodel, channels_to_include,...
    %         components_to_fit, RV_threshold, remove_outside_head,...
    %         fit_bilateral_dipoles);
    
    % finish dataset
    
end

cd(current_dir);

%% Epoching

% this is subject to EEG event structure of the individual study


