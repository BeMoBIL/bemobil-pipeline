% Cleans Euler angle values of a rigid body dataset. Changes in angles above a threshold are assumed to be outliers and
% interpolated using linear interpolation.
%
% Input:
%       EEG_motion                   - EEGLAB dataset containing motion channels with 6 DOF
%       threshold                   - (OPTIONAL) median absolute distance (MAD) multiplier for the detection of bad 
%                                       samples. Samples that are more than threshold * 1.4826* mad(velocities,1) away
%                                       from the median or exactly 0 are assumed to be bad and interpolated                                       
%
% Output:
%       EEG_motion                   - EEGLAB dataset with cleaned Euler angles


function EEG_motion = bemobil_motion_clean_euler(EEG_motion,threshold)

if ~exist('threshold','var') || isempty(threshold)
    threshold = 20;
end

% make sure all euler values are not exceeding pi due to some weirdness in the upsampling

eul_indices = ~cellfun(@isempty,strfind(lower({EEG_motion.chanlocs.labels}),'eul'));
euldata = EEG_motion.data(eul_indices,:);
euldata(euldata>pi)=pi;
euldata(euldata<-pi)=-pi;

% interpolate irregular angular jumps 

for i_dim = 1:3
%     i_dim

    velocities = diff(euldata(i_dim,:),1,2);

    this_threshold = nanmedian(velocities) + threshold*1.4826*mad(velocities,1);
    
    disp(['Cleaning dimension ' num2str(i_dim) ', allowing a maximum of ' num2str(this_threshold*EEG_motion.srate) ' rad/sec as angular velocity...'])
    
    if this_threshold < 0.05 || this_threshold > 0.5
        warning(['threshold for cleaning was ' num2str(this_threshold*EEG_motion.srate) ', which is not suitable for cleaning. Using default of '...
            num2str(0.25*EEG_motion.srate) '!'])
        this_threshold = 0.25;
    end
        
    idx = [false ((velocities < -this_threshold) & ~(velocities < -2*pi - -this_threshold)) |...
        ((velocities > this_threshold) & ~(velocities > 2*pi - this_threshold)) |...
        velocities == 0 ];
    
    % interpolate iteratively until no samples are above the threshold
    iter = 0;
    while sum(idx) > 0 && iter < 20
        iter = iter+1;
%         sum(idx)
        euldata(i_dim,idx) = NaN;

%         disp('Unwrapping data...')
        euldata(i_dim,:) = unwrap(euldata(i_dim,:));
%         disp('...done. Interpolating data.')
        euldata(i_dim,:) = interp1(EEG_motion.times(~idx),euldata(i_dim,~idx),EEG_motion.times,'linear');
        euldata(i_dim,:) = wrapToPi(euldata(i_dim,:));
%         disp('...done.')
    
        % make sure all euler values are not exceeding pi due to some weirdness in the interpolation

        eul_indices = ~cellfun(@isempty,strfind(lower({EEG_motion.chanlocs.labels}),'eul'));
        euldata = EEG_motion.data(eul_indices,:);
        euldata(euldata>pi)=pi;
        euldata(euldata<-pi)=-pi;

        velocities = diff(euldata(i_dim,:),1,2);

        idx = [false ((velocities < -this_threshold) & ~(velocities < -2*pi - -this_threshold)) |...
        ((velocities > this_threshold) & ~(velocities > 2*pi - this_threshold)) |...
        velocities == 0 ];
    end
end

EEG_motion.data(eul_indices,:) = euldata;