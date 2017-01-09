% inspect mobilab stream sampling rates

streams = size(mobilab.allStreams);

for stream = 1:streams(2)
    figure;
    plot(mobilab.allStreams.item{stream}.timeStamp(1:end-1),diff(mobilab.allStreams.item{stream}.timeStamp));
    
    effectiveSamplingRate = length(mobilab.allStreams.item{stream}.timeStamp) / (mobilab.allStreams.item{stream}.timeStamp(end) - mobilab.allStreams.item{stream}.timeStamp(1));
    
    varianceEffectiveSamplingRate = var(diff(mobilab.allStreams.item{stream}.timeStamp));
    
    title(['Stream: ' mobilab.allStreams.item{stream}.name ', Eff SR: ' num2str(effectiveSamplingRate) ' , Variance: ' num2str(varianceEffectiveSamplingRate)]);
end