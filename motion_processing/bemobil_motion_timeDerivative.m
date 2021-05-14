% Computes the time derivatives of motion capture data. It smooths the signals after each order of derivation to
% minimize cumulative precision errors. As smoother it uses a Hamming windowed zero-phase FIR lowpass filter. Each new
% derivative is stored in a new object.
%
% Input:
%       EEG dataset containing motion channels with Euler angles data in radians (6 channels)
%
% Output:
%       EEG dataset containing motion channels and their first and second derivatives (18 channels altogether)

function EEG_motion_out = bemobil_motion_timeDerivative(EEG_motion_in)

% checking for Quaternionvalues

disp('Computing time derivatives...');

for channel = 1:EEG_motion_in.nbchan
    
    % checking for already present euls
    if any(~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat')))
        error('Dataset contains quaternion data. Please transfor to eul before deriving')
    end
end

dt = 1/EEG_motion_in.srate;
tmpData = EEG_motion_in.data;

for channel=1:size(tmpData,1)
    
    % deriving
    tmpData(channel,1:end-1) = diff(tmpData(channel,:),1)/dt;
    tmpData(channel,end) = tmpData(channel,end-1);
    
    % check if channel is eul angles and if so,
    % correct for turns over pi or -pi respectively
    if contains(lower(EEG_motion_in.chanlocs(channel).labels),'eul') && ~contains(lower(EEG_motion_in.chanlocs(channel).labels),'_derivative')
        
        dataChannel = tmpData(channel,:);
        
        assert(max(abs(dataChannel))<2*pi/dt,'Data must be in radian!')
        
        % if there is a jump of more than 180 degrees (pi) between the frames it is assumed to be a turn in the other
        % direction and thus 360 degrees are subtracted. this means turning rates of more than half a circle per frame
        % are not possible.
        dataChannel(dataChannel > pi/dt) = dataChannel(dataChannel > pi/dt) - 2*pi/dt;
        dataChannel(dataChannel < -pi/dt) = dataChannel(dataChannel < -pi/dt) + 2*pi/dt;
        
        tmpData(channel,:) = dataChannel;
        
    end
end

% creating the new object and filling the derived data
EEG_motion_out = EEG_motion_in;
EEG_motion_out.data = tmpData;

for channel = 1:EEG_motion_out.nbchan
    
    EEG_motion_out.chanlocs(channel).labels = [EEG_motion_out.chanlocs(channel).labels '_derivative'];
    
end