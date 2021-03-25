% motion metadata
% copy general fields 
motioncfg       = cfg; 
motioncfg.datatype                                = 'motion';

%--------------------------------------------------------------------------
if ~exist('motionInfo', 'var')

    % data type and acquisition label
    motionInfo.acq                                     = 'Motion';
    
    % motion specific fields in json
    motionInfo.motion.Manufacturer                     = 'Undefined';
    motionInfo.motion.ManufacturersModelName           = 'Undefined';
    motionInfo.motion.RecordingType                    = 'continuous';
    
    % coordinate system
    motionInfo.coordsystem.MotionCoordinateSystem      = 'Undefined';
    motionInfo.coordsystem.MotionRotationRule          = 'Undefined';
    motionInfo.coordsystem.MotionRotationOrder         = 'Undefined';

end


% data type and acquisition label
motioncfg.acq                                     = motionInfo.acq;

% motion specific fields in json
motioncfg.motion                                  = motionInfo.motion;

% coordinate system
motioncfg.coordsystem.MotionCoordinateSystem      = motionInfo.coordsystem;

%--------------------------------------------------------------------------
% rename and fill out motion-specific fields to be used in channels_tsv
motioncfg.channels.name                 = cell(motion.hdr.nChans,1); 
motioncfg.channels.tracked_point        = cell(motion.hdr.nChans,1); 
motioncfg.channels.component            = cell(motion.hdr.nChans,1); 
motioncfg.channels.placement            = cell(motion.hdr.nChans,1); 
motioncfg.channels.datafile             = cell(motion.hdr.nChans,1); 

for ci  = 1:motion.hdr.nChans
    
    if  contains(motion.hdr.chantype{ci},'position')
        motion.hdr.chantype{ci} = 'POS'; 
        motion.hdr.chanunit{ci} = bemobil_config.bids_motion_positionunits{si}; 
    end
    
    if  contains(motion.hdr.chantype{ci},'orientation')
        motion.hdr.chantype{ci} = 'ORNT';
        motion.hdr.chanunit{ci} = bemobil_config.bids_motion_orientationunits{si}; 
    end

    splitlabel                          = regexp(motion.hdr.label{ci}, '_', 'split');
    motioncfg.channels.name{ci}         = motion.hdr.label{ci};
    
    % assign object names and anatomical positions
    for iRB = 1:numel(bemobil_config.rigidbody_streams)
        if contains(motion.hdr.label{ci}, bemobil_config.rigidbody_streams{iRB})
            motioncfg.channels.tracked_point{ci}       = bemobil_config.rigidbody_names{iRB};
            if iscell(bemobil_config.rigidbody_anat)
                motioncfg.channels.placement{ci}  = bemobil_config.rigidbody_anat{iRB};
            else
                motioncfg.channels.placement{ci} =  bemobil_config.rigidbody_anat;
            end
        end
    end
    
    motioncfg.channels.component{ci}    = splitlabel{end};                  % REQUIRED. Component of the representational system that the channel contains.   
    motioncfg.channels.datafile{ci}      = ['...acq-' motioncfg.acq  '_motion.tsv'];                     
    
end