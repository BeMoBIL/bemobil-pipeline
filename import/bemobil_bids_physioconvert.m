
function physioOut = bemobil_bids_physioconvert(physioIn, objects, pi, si, di)
% Process generic physiological data 
% (resampling to the highest sampling rate among all streams of the type)

% iterate over different objects 
%--------------------------------------------------------------------------
physioStreamAll    = cell(numel(physioIn), 1);

for iP = 1:numel(physioIn)

    physioStream                  = physioIn{iP};
    labelsPre                     = [physioStream.label];
    
    dataPre                 = physioStream.trial{1};
    dataPost                = [];
    oi                      = 0;
    
    for ni = 1:numel(objects)
        
        % check first if the object exists at all and if not, skip
        if isempty(find(contains(labelsPre, objects{ni}),1))
            continue;
        else
            oi = oi + 1;
        end

        % no particular processing done here
        dataPost                   = dataPre;
        
    end
    
    % do not include streams that contain 0 objects
    if oi > 0 
        physioStream.trial{1}     = dataPost;
        physioStream.hdr.nChans   = numel(physioStream.hdr.chantype);
        physioStreamAll{iP}       = physioStream;
    end
    
    
end

%--------------------------------------------------------------------------
% find the one with the highest sampling rate
physiosrates = []; 

for iP = 1:numel(physioStreamAll)
    physiosrates(iP) = physioStreamAll{iP}.hdr.Fs; 
end

[~,maxind] = max(physiosrates);

% copy the header from the stream with max srate
keephdr             = physioStreamAll{maxind}.hdr;

% overwrite some fields in header with resampled information 
keephdr.nSamples            = size(physioStreamAll{maxind}.trial{1},2);
keephdr.FirstTimeStamp      = physioStreamAll{maxind}.time{1}(1);
lastTimeStamp               = physioStreamAll{maxind}.time{1}(end);
keephdr.Fs                  = keephdr.nSamples/(lastTimeStamp - keephdr.FirstTimeStamp);
keephdr.TimeStampPerSample  = (lastTimeStamp - keephdr.FirstTimeStamp)/keephdr.nSamples;

% construct evenly spaced time points
regularTime         = {linspace(keephdr.FirstTimeStamp, lastTimeStamp, (lastTimeStamp- keephdr.FirstTimeStamp)*keephdr.Fs)};

keephdr.nChans      = 0;
keephdr.label       = {};
keephdr.chantype    = {};
keephdr.chanunit    = {};

if numel(physioStreamAll)>1
    % resample all data structures, except the one with the max sampling rate
    % this will also align the time axes
    for i=1:numel(physioStreamAll)
        
        % append channel information to the header
        keephdr.nChans      = keephdr.nChans + physioStreamAll{i}.hdr.nChans;
        keephdr.label       = [keephdr.label;       physioStreamAll{i}.hdr.label];
        keephdr.chantype    = [keephdr.chantype;    physioStreamAll{i}.hdr.chantype];
        keephdr.chanunit    = [keephdr.chanunit;    physioStreamAll{i}.hdr.chanunit];
        
        % resample
        %------------------------------------------------------------------ 
        ft_notice('resampling %s', physioStreamAll{i}.hdr.orig.name);
        cfg                 = [];
        cfg.time            = regularTime;
        cfg.detrend         = 'no'; 
        physioStreamAll{i}  = ft_struct2double(physioStreamAll{i});
        physioStreamAll{i}  = ft_resampledata(cfg, physioStreamAll{i});
    end
    
    % append all data structures
    physioOut = ft_appenddata([], physioStreamAll{:});
    
    % modify some fields in the header
    physioOut.hdr = keephdr;
else
    % simply return the first and only one
    physioOut = physioStreamAll{1};
end


end