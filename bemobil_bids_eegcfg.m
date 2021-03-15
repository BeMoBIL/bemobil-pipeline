% eeg specific information 
%--------------------------------------------------------------------------
% eeg configuration for data2bids
eegcfg                              = cfg;
eegcfg.datatype                     = 'eeg';
eegcfg.method                       = 'convert';

% full path to eloc file 
eegcfg.elec                         = fullfile(participantDir, bemobil_config.channel_locations_filename);

% coordinate system
eegcfg.coordsystem.EEGCoordinateSystem      = 'todo';
eegcfg.coordsystem.EEGCoordinateUnits       = 'mm';
