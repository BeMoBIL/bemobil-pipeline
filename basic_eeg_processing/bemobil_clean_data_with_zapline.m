function [EEG, plothandles] = bemobil_clean_data_with_zapline(EEG, zaplineConfig)
% super small wrapper for zapline that makes sure all info is stored and uses a struct as input
%
% see also: clean_data_with_zapline

zaplineFields = fieldnames(zaplineConfig);

zaplineFields(contains(zaplineFields,'linefreqs')) = [];

varargin = {};
for i_fieldname = 1:length(zaplineFields)
    
    varargin{1+(i_fieldname-1)*2} = zaplineFields{i_fieldname};
    varargin{2+(i_fieldname-1)*2} = zaplineConfig.(zaplineFields{i_fieldname});
    
end

% run zapline subfunction with applied parameters
[EEG.data, EEG.etc.zapline.NremoveFinal, EEG.etc.zapline.Scores, EEG.etc.zapline.config, plothandles] =...
    clean_data_with_zapline(EEG.data, EEG.srate, zaplineConfig.linefreqs, varargin{:});
