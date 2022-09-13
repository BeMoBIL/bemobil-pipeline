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


function EEG_motion = bemobil_clean_motion(EEG_motion,threshold)

if ~exist('threshold','var') || isempty(threshold)
    threshold = 20;
end

% make sure all euler values are not exceeding pi due to some weirdness in the upsampling

eul_indices = find(~cellfun(@isempty,strfind(lower({EEG_motion.chanlocs.labels}),'eul')));
euldata = EEG_motion.data(eul_indices,:);

if sum(euldata(:)>pi) > 0
    warning(['Euler angles data contained ' num2str(sum(euldata(:)>pi)) ' samples > pi, with an average value of '...
        num2str(nanmean(euldata(euldata>pi))) ', correcting them all to pi!'])
    euldata(euldata>pi)=pi;
end
if sum(euldata(:)<-pi) > 0
    warning(['Euler angles data contained ' num2str(sum(euldata(:)<-pi)) ' samples < -pi, with an average value of '...
        num2str(nanmean(euldata(euldata>pi))) ', correcting them all to -pi!'])
    euldata(euldata<-pi)=-pi;
end

% interpolate irregular angular jumps 

for i_dim = 1:3
%     i_dim

    velocities = diff(euldata(i_dim,:),1,2);

    this_threshold = nanmedian(velocities) + threshold*1.4826*mad(velocities,1);
    
    disp(['Cleaning and interpolating ' EEG_motion.chanlocs(eul_indices(i_dim)).labels ', allowing a maximum of ' num2str(this_threshold*EEG_motion.srate) ' rad/sec as angular velocity...'])
    
    minthresh = 0.04;
    maxthresh = 0.4;
    if this_threshold < minthresh
        warning(['threshold for cleaning was ' num2str(this_threshold*EEG_motion.srate) 'rad/sec as angular velocity, which is not suitable for cleaning. Using minimum of '...
            num2str(minthresh*EEG_motion.srate) '!'])
        this_threshold = minthresh;
    elseif this_threshold > maxthresh
        warning(['threshold for cleaning was ' num2str(this_threshold*EEG_motion.srate) 'rad/sec as angular velocity, which is not suitable for cleaning. Using maximum of '...
            num2str(maxthresh*EEG_motion.srate) '!'])
        this_threshold = maxthresh;
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
        euldata(i_dim,:) = fillmissing(euldata(i_dim,:),'pchip','SamplePoints',EEG_motion.times,'endvalues','nearest');
        euldata(i_dim,:) = wrapToPi(euldata(i_dim,:));
%         disp('...done.')
    
        % make sure all euler values are not exceeding pi due to some weirdness in the interpolation
        if sum(euldata(:)>pi) > 0
            warning(['Euler angles data contained ' num2str(sum(euldata(:)>pi)) ' samples > pi after interpolation, with an average value of '...
                num2str(nanmean(euldata(euldata>pi))) ', correcting them all to pi!'])
            euldata(euldata>pi)=pi;
        end
        if sum(euldata(:)<-pi) > 0
            warning(['Euler angles data contained ' num2str(sum(euldata(:)<-pi)) ' samples < -pi after interpolation, with an average value of '...
                num2str(nanmean(euldata(euldata>pi))) ', correcting them all to -pi!'])
            euldata(euldata<-pi)=-pi;
        end

        velocities = diff(euldata(i_dim,:),1,2);

        idx = [false ((velocities < -this_threshold) & ~(velocities < -2*pi - -this_threshold)) |...
        ((velocities > this_threshold) & ~(velocities > 2*pi - this_threshold)) |...
        velocities == 0 ];
    end
end

% interpolate non euler data
indices = find(cellfun(@isempty,strfind(lower({EEG_motion.chanlocs.labels}),'eul')));

for i_chan = indices
    
    disp(['Interpolating ' EEG_motion.chanlocs(i_chan).labels '...'])
    EEG_motion.data(i_chan,:) = fillmissing(EEG_motion.data(i_chan,:),'pchip','SamplePoints',EEG_motion.times,'endvalues','nearest');
    
end

EEG_motion.data(eul_indices,:) = euldata;