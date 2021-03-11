% information needed to construct file paths and names 
cfg.sub                                     = num2str(participantNr);
cfg.dataset                                 = fullfile(participantDir, sortedFileNames{di});
cfg.ses                                     = bemobil_config.filenames{si}; 
cfg.run                                     = di; 

% remove session label in uni-session case
if numel(bemobil_config.filenames) == 1
    cfg = rmfield(cfg, 'ses');
end

% remove session label in uni-run case
if numel(sortedFileNames) == 1
    cfg = rmfield(cfg, 'run');
end

% participant information
cfg.participants.age                        = 28;        
cfg.participants.sex                        = 'F'; 
cfg.participants.handedness                 = 'R'; % warning : this is only a placeholder for now - more to be implemented