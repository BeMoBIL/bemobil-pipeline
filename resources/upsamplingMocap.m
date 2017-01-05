% upsampling von mocap daten
mocapStream = mobilab.allStreams.item{10}; % head velocity values
effectiveSamplingRate = 1/mean(diff(mocapStream.timeStamp));

newSamplingRate = 250;

newTimestamps = 1/newSamplingRate:1/newSamplingRate:mocapStream.timeStamp(end);

% length(newTimestamps)/length(mocapStream.timeStamp)
% mocapStream.timeStamp(end)*effectiveSamplingRate - length(mocapStream.timeStamp)

x = mocapStream.timeStamp;
y = mocapStream.data(:,4); % yaw values
xx = newTimestamps;
newData = spline(x,y,xx);