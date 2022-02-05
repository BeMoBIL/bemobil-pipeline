
function motionOut = bemobil_bids_motionconvert(motionIn, objects, trackSysConfig)
% This function performs minimal preprocessing and channel sorting for motion data 
%--------------------------------------------------------------------------

disp(['Coverting motion data from tracking system ' trackSysConfig.name ' to BIDS format'])

% Finding quaternion data
%--------------------------------------------------------------------------
% method 1 : keyword + components
% trackSysConfig.quaternions.keyword              = '_quat'; 
% trackSysConfig.quaternions.components           = {'w', 'x', 'y', 'z'};   % components are assumed to follow an "_", e.g., "quat_w"
%                                                                           % if this rule is violated, use channel_names option
%
% method 2 : channel_names 
% trackSysConfig.quaternions.channel_names        = {'headRigid_rotW', 'rightHand_rotW';, ... 
%                                                    'headRigid_rotX', 'rightHand_rotX';, ...
%                                                    'headRigid_rotY', 'rightHand_rotY';, ...
%                                                    'headRigid_rotZ', 'rightHand_rotZ'}; 
%
% if nothing is found using these methods, no stream is processed as
% quaternion
% search by channel names is prioritized

% quaternion [w,x,y,z] components, in this order
% w is the non-axial component
quaternionComponents    = {'w','x','y','z'}; % default values 
quaternionKeyword       = 'quat';
quaternionNames         = {};
eulerComponents         = {'x','y','z'};

if isfield(trackSysConfig, 'quaternions')
    if isfield(trackSysConfig.quaternions, 'channel_names')
        quaternionNames    = trackSysConfig.quaternions.channel_names;
    end
    
    if isfield(trackSysConfig.quaternions, 'keyword')
        quaternionKeyword = trackSysConfig.quaternions.keyword; 
    end
    
    if isfield(trackSysConfig.quaternions, 'components')
        quaternionComponents = trackSysConfig.quaternions.components; 
    end
    
    if isfield(trackSysConfig.quaternions, 'output_order')
        eulerComponents         = trackSysConfig.euler_components;
    end
end

% Finding position data
%--------------------------------------------------------------------------
% method 1 : keyword + components
% trackSysConfig.positions.keyword              = '_pos_'; 
% trackSysConfig.positions.components           = {'x', 'y', 'z'}; 
%
% method 2 : channel_names 
% trackSysConfig.quaternions.channel_names        = {'headRigid_posX', 'rightHandX';, ... 
%                                                    'headRigid_posY', 'rightHandY';, ...
%                                                    'headRigid_posZ', 'rightHandZ'}; 
%
% if nothing is found using these methods, no stream is processed as
% quaternion
% search by channel names is prioritized

% cartesian [x,y,z] components, in this order
% 1-D or 2-D cases are not implemented yet
cartCoordinates    = {'x','y','z'}; % default values
cartKeyword       = 'pos';
cartNames         = {}; 

if isfield(trackSysConfig, 'positions')
    if isfield(trackSysConfig.positions, 'channel_names')
        cartNames      = trackSysConfig.positions.channel_names;
    end
    
    if isfield(trackSysConfig.positions, 'keyword')
        cartKeyword    = trackSysConfig.positions.keyword; 
    end
    
    if isfield(trackSysConfig.positions, 'components')
        cartCoordinates = trackSysConfig.positions.components; 
    end
    
end

% missing value (how tracking loss is represented in the stream)
if isfield(trackSysConfig, 'missing_values')
    switch trackSysConfig.missing_values
        case '0'
            missingval = 0;
        case 'NaN'
            missingval = NaN;
        otherwise
            warning(['Unrecognized value for field "missing_values" in tracking system ' trackSysConfig.name ': it should be "0" or "NaN" formatted as string.'])
            warning('Taking default value NaN for missing samples.')
            missingval = NaN;
    end
else
    missingval = NaN;
end

newCell = {}; 
% check if object input is a nested cell
for Oi = 1:numel(objects)
   if iscell(objects(Oi))
       newCell = [newCell objects{Oi}];
   else
       newCell{end + 1} = objects(Oi); 
   end
end
objects = newCell; 

% iterate over different objects 
%--------------------------------------------------------------------------
motionStreamAll    = cell(numel(motionIn), 1);

for iM = 1:numel(motionIn)

    motionStream                  = motionIn{iM};
    labelsPre                     = [motionStream.label];
    motionStream.label            = [];
    motionStream.hdr.label        = [];
    motionStream.hdr.chantype     = [];
    motionStream.hdr.chanunit     = [];
    
    dataPre                 = motionStream.trial{1};
    dataPost                = [];
    oi                      = 0;
    
    for ni = 1:numel(objects)
        
        % find all channels from the given object 
        objectChans         = find(contains(labelsPre, objects{ni})); 
        
        % check first if the object exists at all and if not, skip
        if isempty(objectChans)
            continue;
        else
            oi = oi + 1;
        end
        
        quaternionIndices = NaN(1,4);
        
        % find quaternion channels
        %------------------------------------------------------------------
        quatFound = false; 
        if ~isempty(quaternionNames)
            try
                for qi = 1:4
                    quaternionIndices(qi) = find(contains(labelsPre, quaternionNames{qi,:},'IgnoreCase', true));
                end
                quatFound = true; 
            catch
                try
                    for qi = 1:4
                        quaternionIndices(qi) = find(strcmpi(labelsPre, quaternionNames{qi,:}));
                    end
                    quatFound = true;
                catch
                    warning('Quaternion names for the tracking system were specified but could not be used.')
                end
            end
        end
        
        if ~quatFound
            try
                for qi = 1:4
                    quaternionIndices(qi) = find(contains(labelsPre, objects{ni}) & contains(labelsPre, ['_' quaternionComponents{qi}], 'IgnoreCase', true) & contains(labelsPre, quaternionKeyword, 'IgnoreCase', true));
                end
                quatFound = true; 
            catch
                warning('No quaternion data found')
            end
        end
        
        cartIndices = NaN(1,3);
        
        % find position channels
        %------------------------------------------------------------------
        cartFound = false; 
        if ~isempty(cartNames)
            try
                for ci = 1:3
                    cartIndices(ci) = find(contains(labelsPre, cartNames{ci,:},'IgnoreCase', true));
                end
                cartFound = true; 
            catch
                try
                    for ci = 1:3
                        cartIndices(ci) = find(strcmpi(labelsPre, cartNames{ci,:}));
                    end
                    cartFound = true;
                catch
                    warning('Cartesian coordinates for the tracking system were specified but could not be used.')
                end
            end
        end
        
        if ~cartFound
            try
                for ci = 1:3
                    cartIndices(ci) = find(contains(labelsPre, objects{ni}) & contains(labelsPre, ['_' cartCoordinates{ci}], 'IgnoreCase', true) & contains(labelsPre, cartKeyword, 'IgnoreCase', true));
                end
                cartFound = true;
            catch
                warning('No cartesian position data found')
            end
        end
        
        if quatFound
            % convert from quaternions to euler angles
            orientationInQuaternion    = dataPre(quaternionIndices,:)';
            orientationInEuler         = util_quat2eul(orientationInQuaternion); % the BeMoBIL util script
            orientationInEuler         = orientationInEuler';
            occindices                 = find(orientationInQuaternion(1,:) == missingval);
            orientationInEuler(:,occindices) = nan;
            orientationInEuler  = fillmissing(orientationInEuler', 'nearest')';
            
            % unwrap euler angles
            orientationInEuler  = unwrap(orientationInEuler, [], 2);
        else
            orientationInEuler = []; 
            quaternionIndices = []; 
        end
        
        if cartFound
            % find and fill missing values
            position                    = dataPre(cartIndices,:);
            occindices                  = find(position(1,:) == missingval);
            position(:,occindices)      = nan;
            position            = fillmissing(position', 'pchip')';
        else
           position     = []; 
           cartIndices  = [];  
        end
        
        % construct a latency channel
        latency  = motionStream.time{1};
        
        % all other channels 
        otherChans          = setdiff(objectChans, union(cartIndices, quaternionIndices));
        otherData           = dataPre(otherChans,:);
        
        % concatenate the converted data
        objectData         = [orientationInEuler; position; otherData; latency];
        dataPost           = [dataPost; objectData];
        
        % enter channel information
        for ei = 1:size(orientationInEuler,1)
            motionStream.label{end + 1}                 = [objects{ni} '_eul_' eulerComponents{ei}];
            motionStream.hdr.label{end + 1}             = [objects{ni} '_eul_' eulerComponents{ei}];
            motionStream.hdr.chantype{end + 1}          = 'ORNT';
            motionStream.hdr.chanunit{end + 1}          = 'rad';
        end
        
        for ci = 1:size(position,1)
            motionStream.label{end + 1}                 = [objects{ni} '_cart_' cartCoordinates{ci}];
            motionStream.hdr.label{end + 1}             = [objects{ni} '_cart_' cartCoordinates{ci}];
            motionStream.hdr.chantype{end + 1}          = 'POS';
            motionStream.hdr.chanunit{end + 1}          = 'm';
        end
        
        for ii = 1:size(otherData,1)
            motionStream.label{end + 1}                 = motionStream.hdr.orig.desc.channels.channel{otherChans(ii)}.label;
            motionStream.hdr.label{end + 1}             = motionStream.hdr.orig.desc.channels.channel{otherChans(ii)}.label;
            motionStream.hdr.chantype{end + 1}          = motionStream.hdr.orig.desc.channels.channel{otherChans(ii)}.type;
            if isfield( motionStream.hdr.orig.desc.channels.channel{otherChans(ii)}, 'unit')
                motionStream.hdr.chanunit{end + 1}          = motionStream.hdr.orig.desc.channels.channel{otherChans(ii)}.unit;
            else
                motionStream.hdr.chanunit{end + 1}          = 'n/a';
            end
        end
        
        motionStream.label{end + 1}                 = [objects{ni} '_latency'];
        motionStream.hdr.label{end + 1}             = [objects{ni} '_latency'];
        motionStream.hdr.chantype{end + 1}          = 'LATENCY';
        motionStream.hdr.chanunit{end + 1}          = 'seconds';
        
    end
    
    % only include streams that have data from at least one object 
    if oi > 0 
        motionStream.trial{1}     = dataPost;
        motionStream.hdr.nChans   = numel(motionStream.hdr.chantype);
        
        % add time stamps channel to the stream before concatenating
        motionStreamAll{iM}       = motionStream;
    end
    
end

%--------------------------------------------------------------------------
% find the one with the highest sampling rate 
% & check if data can be reasonably concatenated
motionsrates    = []; 
nSamples        = []; 
for iM = 1:numel(motionStreamAll)
    motionsrates(iM) = motionStreamAll{iM}.hdr.Fs; 
    nSamples(iM)     = motionStreamAll{iM}.hdr.nSamples; 
end

if all(nSamples(:) == nSamples(1))
    doResample      = ~strcmp(trackSysConfig.keep_timestamps, 'on'); 
else
    doResample      = true;
    if strcmp(trackSysConfig.keep_timestamps, 'on') 
        warning('Config field keep_timestamps was specified as "on", but number of samples among streams in tracksys do not match - resampling...')
    end
end

[~,maxind] = max(motionsrates);

% copy the header from the stream with max srate
keephdr             = motionStreamAll{maxind}.hdr;

% overwrite some fields in header with resampled information 
keephdr.nSamples            = size(motionStreamAll{maxind}.trial{1},2);
keephdr.FirstTimeStamp      = motionStreamAll{maxind}.time{1}(1);
lastTimeStamp               = motionStreamAll{maxind}.time{1}(end);
keephdr.Fs                  = keephdr.nSamples/(lastTimeStamp - keephdr.FirstTimeStamp);
keephdr.TimeStampPerSample  = (lastTimeStamp - keephdr.FirstTimeStamp)/keephdr.nSamples;

% construct evenly spaced time points
regularTime         = {linspace(keephdr.FirstTimeStamp, lastTimeStamp, (lastTimeStamp- keephdr.FirstTimeStamp)*keephdr.Fs)};

keephdr.nChans      = 0;
keephdr.label       = {};
keephdr.chantype    = {};
keephdr.chanunit    = {};

if numel(motionStreamAll)> 1 
    
    % resample all data structures by interpolation, this will also align the time axes
    for i=1:numel(motionStreamAll)
        
        % append channel information to the header
        keephdr.nChans      = keephdr.nChans + motionStreamAll{i}.hdr.nChans;
        keephdr.label       = [keephdr.label;       motionStreamAll{i}.hdr.label];
        keephdr.chantype    = [keephdr.chantype;    motionStreamAll{i}.hdr.chantype];
        keephdr.chanunit    = [keephdr.chanunit;    motionStreamAll{i}.hdr.chanunit];
        
        if doResample
            % resample
            %--------------------------------------------------------------
            ft_notice('resampling %s', motionStreamAll{i}.hdr.orig.name);
            cfg                 = [];
            cfg.time            = regularTime;
            cfg.detrend         = 'no';
            motionStreamAll{i} = ft_resampledata(cfg, motionStreamAll{i});
        end
        
    end

    % append all data structures
    motionOut = ft_appenddata([], motionStreamAll{:});
    
    % modify some fields in the header
    motionOut.hdr = keephdr;
else
    % simply return the first and only one
    motionOut = motionStreamAll{1};
end

% wrap all orientation channels back to [pi, -pi]
%--------------------------------------------------------------------------
% find orientation channels
for Ci = 1:numel(motionOut.label)
    if contains(motionOut.label{Ci}, '_eul_')
        motionOut.trial{1}(Ci,:) = wrapToPi(motionOut.trial{1}(Ci,:)); 
    end
end

end