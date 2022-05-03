
function physioOut = bemobil_bids_physioconvert(physioIn, objects, interpPhys)
% Process generic physiological data 
% (resampling to the highest sampling rate among all streams of the type)

% iterate over different objects 
%--------------------------------------------------------------------------
physioStreamAll    = cell(numel(physioIn), 1);

for iP = 1:numel(physioIn)

    physioStream                  = physioIn{iP};
    labelsPre                     = [physioStream.hdr.orig.name];
    
    dataPre                 = double(physioStream.trial{1});
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
        
        % construct a latency channel
        latency  = physioStream.time{1};
        dataPost = [dataPost; latency];
        physioStream.label{end + 1}                 = [objects{ni} '_latency'];
        physioStream.hdr.label{end + 1}             = [objects{ni} '_latency'];
        physioStream.hdr.chantype{end + 1}          = 'LATENCY';
        physioStream.hdr.chanunit{end + 1}          = 'seconds';
        
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


if numel(physioStreamAll) == 1
    % if there is one physio stream, check the config field to decide
    % whether to resample or keep timestamps 
    doResample      = interpPhys; 
else
    doResample      = true;
    if interpPhys 
        warning('Config field phys.skip_interp was specified as true, but there are multiple streams to be concatenated - resampling via interpolation')
    end
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
        
        if doResample
            % resample
            %------------------------------------------------------------------
            ft_notice('resampling %s', physioStreamAll{i}.hdr.orig.name);
            cfg                 = [];
            cfg.time            = regularTime;
            cfg.detrend         = 'no';
            physioStreamAll{i}  = ft_struct2double(physioStreamAll{i});
            physioStreamAll{i}  = ft_resampledata(cfg, physioStreamAll{i});
            
            % remove excess latency channel now after resampling, it will crash when appending
            trial = physioStreamAll{i}.trial{1};
            trial = trial(1:end-1,:);
            physioStreamAll{i}.trial = {trial};
            physioStreamAll{i}.label = physioStreamAll{i}.label(1:end-1);
            
        end
    end
    
    % append all data structures
    physioOut = ft_appenddata([], physioStreamAll{:});
    
    % modify some fields in the header
    physioOut.hdr = keephdr;
    
    % add latency channel for later interpolation
    physioOut.label = [physioOut.label; 'physio_latency'];
    trial = physioOut.trial{1};
    trial = [trial;regularTime{1}];
    physioOut.trial = {trial};
    
else
    % simply return the first and only one
    physioOut = physioStreamAll{1};
end


end