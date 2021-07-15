
%% dot indexing problem in motionconvert [solved]

% streams import spotrotation
streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\spotrotation\vp-6\vp-6_control_body.xdf');
% streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\spotrotation\vp-6\vp-6_control_joy.xdf');
streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\spotrotation\vp-7\vp-7_control_joy.xdf');
% select stream 
 % stream body
 xdfstream = streams{4};
 
 % stream joy
 xdfstream = streams{3};

 streamNr = 3; 

% load .xdf data to check what is in there
streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\data\spotrotation\vp-7\vp-7_control_joy.xdf');
streamnames     = cellfun(@(x) x.info.name, streams, 'UniformOutput', 0)'
channelnames    = cellfun(@(x) x.label, streams{streamNr}.info.desc.channels.channel, 'UniformOutput', 0)'

% visualize position streams 
figure; plot(streams{streamNr}.time_series(1:3,:)', 'LineWidth', 2)
set(gca,'FontSize',15) 
legend('X', 'Y', 'Z')
title(['Stream ' streamnames{streamNr} ' position streams'], 'Interpreter','none')

% visualize orientation streams 
figure; plot(streams{streamNr}.time_series(4:7,:)', 'LineWidth', 2)
set(gca,'FontSize',15) 
legend('A', 'B', 'C', 'D')
title(['Stream ' streamnames{streamNr} ' quaternion streams'], 'Interpreter','none')

% % streams import mobiworkshop
% streams = load_xdf('C:\Users\sgrot\Documents\Uni\01_Master\6. Semester\00_MA\Themenfindung\BIDs\MoBI Workshop\data\0_source-data\vp_24\vp_24_walk.xdf');
% 
% % select stream 
% xdfstream = streams{1};


%--------------------------------------------------------------
%                Convert Motion Data to BIDS
%--------------------------------------------------------------

motionsrates = []; 
ftmotion = [];

% construct header
% hdr.Fs                  = xdfstream.info.effective_srate; % depends on calculating .effective_srate which is not important atm
hdr.nominalSampleRate   = str2num(xdfstream.info.nominal_srate);
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
bemobil_config.rigidbody_streams        = {'headRigid', 'vizardFirstControllerRigid'};
motionStreamNames = {'headRigid', 'vizardFirstControllerRigid'}
objects = motionStreamNames;

motion                  = motionIn;
labelsPre               = [motion.label];
motion.label            = [];

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


ftmotion{1} = ftdata
            
            % construct fieldtrip data
            for iM = 1:numel(xdfmotion)
%                 ftmotion{iM} = stream2ft(xdfmotion{iM}); 
           
                % if needed, execute a custom function for any alteration to the data to address dataset specific issues
                % (quat2eul conversion, unwrapping of angles, resampling, wrapping back to [pi, -pi], and concatenating for instance)
                motion = feval(motionCustom, ftmotion{iM}, {'headRigid'}, participantNr, si, di);
                
            end

%% construct tracked points

bemobil_config.rigidbody_streams        = {'headRigid', 'vizardFirstControllerRigid'};
bemobil_config.rigidbody_names          =  {'Head', 'JoyStickRotation'}; 
bemobil_config.rigidbody_anat           = {'head', 'n/a'}; 
motioncfg.channels.placement            = [];
 

for ci = numel(motionOut.hdr.label)
                splitlabel                          = regexp(hdr.label{ci}, '_', 'split');
                motioncfg.channels.name{ci}         = motionOut.hdr.label{ci};
    
                % assign object names and anatomical positions
                for iRB = 1:numel(bemobil_config.rigidbody_streams)
                    if contains(motionOut.hdr.label{ci}, bemobil_config.rigidbody_streams{iRB})
                        motioncfg.channels.tracked_point{ci}       = bemobil_config.rigidbody_names{iRB};
                        if iscell(bemobil_config.rigidbody_anat)
                            motioncfg.channels.placement{ci}  = bemobil_config.rigidbody_anat{iRB};
                        else
                            motioncfg.channels.placement{ci} =  bemobil_config.rigidbody_anat;
                        end
                    end
                end
                
                motioncfg.channels.component{ci}    = splitlabel{end};                  % REQUIRED. Component of the representational system that the channel contains.
%                 motioncfg.channels.datafile{ci}      = ['...acq-' motioncfg.acq  '_motion.tsv'];
  
end 

for ci = numel(motionOut.hdr.label)
                motioncfg.channels.name{ci}         = motionOut.hdr.label{ci};
end 

%% count tracked points
trackedpoints = numel(bemobil_config.rigidbody_anat) - sum(strcmpi(bemobil_config.rigidbody_anat, 'n/a'))



%%

sort({'SamplingFrequencyNominal' 'SoftwareFilters'  'SamplingFrequencyEffective'})