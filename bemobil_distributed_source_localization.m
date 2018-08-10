function weights = bemobil_distributed_source_localization(EEG, pattern, plot_threshold)

figure;topoplot(pattern,EEG.chanlocs);title('Pattern');

ica_weights_combined = EEG.icaweights * EEG.icasphere;
% ica_weights_combined = EEG.icawinv';

weights = abs(ica_weights_combined * pattern);

weights = weights/max(weights);
weights(weights<1e-10) = 1e-10;

relevant_pattern_indices = weights>plot_threshold;
relevant_patterns = find(relevant_pattern_indices);

bemobil_plot_patterns(EEG.icawinv(:,relevant_patterns), EEG.chanlocs, 'weights', weights(relevant_patterns), 'titles' , strcat(repmat({'IC '},[length(relevant_patterns) 1]),cellstr(num2str(relevant_patterns))));

dipolelength = weights(1); 
if dipolelength < 0.1, dipolelength = 0; end;
fprintf('Plotting dipole 1 ')
pop_dipplot( EEG,relevant_patterns(1),'mri','P:\\Marius\\toolboxes\\eeglab14_1_0b\\plugins\\dipfit2.3\\standard_BEM\\standard_mri.mat','projlines','on','dipolelength',0,'dipolesize',weights(relevant_patterns(1))*50,'spheres','on');

for weight = relevant_patterns(2:end)
    fprintf('%d ',weight)
    
    dipolelength = weights(weight); if dipolelength < 0.1, dipolelength = 0; end;
    
    pop_dipplot( EEG,weight ,'mri','P:\\Marius\\toolboxes\\eeglab14_1_0b\\plugins\\dipfit2.3\\standard_BEM\\standard_mri.mat','projlines','on','dipolelength',0,'dipolesize',weights(weight)*50,'holdon','on','spheres','on');

end
fprintf('\n')