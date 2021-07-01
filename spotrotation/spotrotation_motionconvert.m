
function motionOut = spotrotation_motionconvert(motionIn, objects, pi, si, di)

quaternionComponents    = {'X','Y','Z','W'};
eulerComponents         = {'z','y','x'}; 
cartCoordinates         = {'X','Y','Z'};
   
motion                  = motionIn; 
labelsPre               = [motion.label];
motion.label            = [];
motion.hdr.label        = []; 
motion.hdr.chantype     = []; 
motion.hdr.chanunit     = []; 

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

end