% streams test
streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\spotrotation\vp-6\vp-6_control_body.xdf');

% select stream 
xdfstream = streams{6};

%--------------------------------------------------------------
%                Convert Motion Data to BIDS
%--------------------------------------------------------------

motionsrates = []; 
ftmotion = [];

% construct header
% hdr.Fs                  = xdfstream.info.effective_srate;
hdr.nSamplesPre         = 0;
hdr.nSamples            = length(xdfstream.time_stamps);
hdr.nTrials             = 1;
hdr.FirstTimeStamp      = xdfstream.time_stamps(1);
hdr.TimeStampPerSample  = (xdfstream.time_stamps(end)-xdfstream.time_stamps(1)) / (length(xdfstream.time_stamps) - 1);

if isfield(xdfstream.info.desc, 'channels')
    hdr.nChans    = numel(xdfstream.info.desc.channels.channel);
else
    hdr.nChans    = str2double(xdfstream.info.channel_count);
end

hdr.label       = cell(hdr.nChans, 1);
hdr.chantype    = cell(hdr.nChans, 1);
hdr.chanunit    = cell(hdr.nChans, 1);

prefix = xdfstream.info.name;
for j=1:hdr.nChans
    if isfield(xdfstream.info.desc, 'channels')
        hdr.label{j} = [prefix '_' xdfstream.info.desc.channels.channel{j}.label];
        hdr.chantype{j} = xdfstream.info.desc.channels.channel{j}.type;
        try
            hdr.chanunit{j} = xdfstream.info.desc.channels.channel{j}.unit;
        catch
            disp([hdr.label{j} ' missing unit'])
        end
    else
        % the stream does not contain continuously sampled data
        hdr.label{j} = num2str(j);
        hdr.chantype{j} = 'unknown';
        hdr.chanunit{j} = 'unknown';
    end


% keep the original header details
hdr.orig = xdfstream.info;

ftdata.trial    = {xdfstream.time_series};
ftdata.time     = {xdfstream.time_stamps};
ftdata.hdr = hdr;
ftdata.label = hdr.label;

end


% spotrotation motionconvert
quaternionComponents    = {'X','Y','Z','W'};
eulerComponents         = {'z','y','x'}; 
cartCoordinates         = {'X','Y','Z'};

motionIn = ftdata;

motion                  = motionIn;
motion.label            = [];
labelsPre               = [motion.label];


motion.hdr.label        = []; 
motion.hdr.chantype     = []; 
motion.hdr.chanunit     = []; 

dataPre                 = motion.trial{1};
dataPost                = []; 
oi                      = 0; 

dataPre                 = motion.trial{1};
dataPost                = []; 
oi                      = 0; 

for ni = 1:numel(objects)
    
    % check first if the object exists at all and if not, skip
    if isempty(find(contains(labelsPre, objects{ni}),1))
        continue; 
    else 
        oi = oi + 1; 
    end
    
    quaternionIndices = NaN(1,4); 
    
    for qi = 1:4
        quaternionIndices(qi) = find(contains(labelsPre, [objects{ni} '_quat_' quaternionComponents{qi}]));        
    end
    
    cartIndices = NaN(1,3); 
    
    for ci = 1:3
        cartIndices(ci) = find(contains(labelsPre, [objects{ni} '_' cartCoordinates{ci}]));
    end
    
    % convert from quaternions to euler angles
    orientationInQuaternion    = dataPre(quaternionIndices,:)';
    orientationInEuler         = util_quat2eul(orientationInQuaternion);    % the BeMoBIL util script
    orientationInEuler         = orientationInEuler';
    position                   = dataPre(cartIndices,:);
    
    % concatenate the converted data 
    dataPost                   = [dataPost; orientationInEuler; position];
    
    % enter channel information 
    for ei = 1:3
        motion.label{6*(oi-1) + ei}                 = [objects{ni} '_eul_' eulerComponents{ei}];
        motion.hdr.label{6*(oi-1) + ei}             = [objects{ni} '_eul_' eulerComponents{ei}];
        motion.hdr.chantype{6*(oi-1) + ei}          = 'orientation';
        % to do : convert the orientation data to radian 
        motion.hdr.chanunit{6*(oi-1) + ei}         = 'radian';
    end
    
    for ci = 1:3
        motion.label{6*(oi-1) + 3 + ci}                 = [objects{ni} '_cart_' cartCoordinates{ci}];
        motion.hdr.label{6*(oi-1) + 3 + ci}             = [objects{ni} '_cart_' cartCoordinates{ci}];
        motion.hdr.chantype{6*(oi-1) + 3 + ci}          = 'position';
        motion.hdr.chanunit{6*(oi-1) + 3 + ci}          = 'meters';
    end
    
end

motion.trial{1}     = dataPost;
motion.hdr.nChans   = numel(motion.hdr.chantype);
motionOut           = motion; 

            
