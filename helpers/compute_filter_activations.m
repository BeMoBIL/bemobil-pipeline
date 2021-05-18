function [ EEG ] = compute_filter_activations(EEG)
%Computes activations of EEG component activity

    %Project channel data to source space
    data_dimensions = size(EEG.data); 
    
    if ndims(EEG.data) > 2
        %epoched data. Flatten and project
        EEG.data = reshape(EEG.data, [data_dimensions(1), data_dimensions(2) * data_dimensions(3)]);
        EEG.icaact = (EEG.data' *  EEG.spocweights)';
        EEG.icaact = reshape(EEG.icaact, [size(EEG.icaact,1), data_dimensions(2), data_dimensions(3)]);
        EEG.data = reshape(EEG.data, data_dimensions);
    else
        EEG.icaact = (EEG.data' * EEG.spocweights)';
% 		EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end

end

