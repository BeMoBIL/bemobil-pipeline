function EEG_mocap_out = bemobil_mocap_timeDerivative(EEG_mocap_in)
% Computes the time derivatives of motion capture data. It smooths
% the signals after each order of derivation to minimize cumulative
% precision errors. As smoother it uses a Hamming windowed zero-phase FIR lowpass filter.
% Each new derivative is stored in a new object.
%
% Input arguments:
%       order:     maximum order of derivation, default: 3 (1 = velocity,
%                  2 = acceleration, 3 = jerk)
%       cutOff:    lowpass filter cutoff, default: 18 Hz.
%
% Output argument:
%       cobj:      handle to the object containing the latest order of
%                  derivation


% checking for Quaternionvalues

for channel = 1:EEG_mocap_in.nbchan
    
    % checking for already present eulers
    if any(~cellfun(@isempty,strfind(lower({EEG_mocap_in.chanlocs.labels}),'quat')))
        error('Dataset contains quaternion data. Please transfor to euler before deriving')
    end
end

dt = 1/EEG_mocap_in.srate;
tmpData = EEG_mocap_in.data;

for channel=1:size(tmpData,1)
    
    % deriving
    tmpData(channel,1:end-1) = diff(tmpData(channel,:),1)/dt;
    tmpData(channel,end) = tmpData(channel,end-1);
    
    % check if channel is Euler angles and if so,
    % correct for turns over pi or -pi respectively
    if contains(lower(EEG_mocap_in.chanlocs(channel).labels),'eul')
        
        dataChannel = tmpData(channel,:);
        
        % if there is a jump of almost 360 degrees between the frames it is assumed to be actually a smooth transition
        % and thus 360 degrees are subtracted. this means turning rates of more than half a circle per frame are not
        % possible
        dataChannel(dataChannel > 180/dt) = dataChannel(dataChannel > 180/dt) - 2*180/dt;
        dataChannel(dataChannel < -180/dt) = dataChannel(dataChannel < -180/dt) + 2*180/dt;
        
        tmpData(channel,:) = dataChannel;
        
    end
end

% creating the new object and filling the derived data
EEG_mocap_out = EEG_mocap_in;
EEG_mocap_out.data = tmpData;

for channel = 1:EEG_mocap_out.nbchan
    
    EEG_mocap_out.chanlocs(channel).labels = [EEG_mocap_out.chanlocs(channel).labels '_derivative'];
    
end