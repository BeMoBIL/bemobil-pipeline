% inspect xdf stream sampling rates

% streamsTrue = load_xdf('P:\Marius\Test Mocap\10\test_body_irregularRate_2.xdf', 'HandleJitterRemoval', true, 'Verbose', true);
% streamsFalse = load_xdf('P:\Marius\Test Mocap\10\test_body_irregularRate_2.xdf', 'HandleJitterRemoval', false, 'Verbose', true);

streamsTrue = load_xdf('P:\Lukas_Gehrke\studies\Spot_Rotation\data\level_0\5\test_body.xdf', 'HandleJitterRemoval', true, 'Verbose', true);
% streamsFalse = load_xdf('P:\Lukas_Gehrke\studies\Spot_Rotation\data\level_0\4\test_body.xdf', 'HandleJitterRemoval', false, 'Verbose', true);

streams = streamsTrue;
for stream = 1:length(streams)
    figure;
    plot(streams{stream}.time_stamps(1:end-1),diff(streams{stream}.time_stamps));
    
    effectiveSamplingRate = streams{stream}.info.effective_srate;
    
    varianceEffectiveSamplingRate = var(diff(streams{stream}.time_stamps));
    
    title(['Stream: ' streams{stream}.info.name ', Eff SR: ' num2str(effectiveSamplingRate) ' , Variance: ' num2str(varianceEffectiveSamplingRate)]);
end