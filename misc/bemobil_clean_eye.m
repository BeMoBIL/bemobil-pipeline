% bemobil_clean_eye cleans eye gaze data by detecting blinks based on pupil radius data, and interpolating them using
% pchip interpolation. Eye blink events are added to the events of the data sets and blink extraction information is
% stored in the etc struct field of the data set.
%
% The blink detection is modified from code written by Ravi Chacko with support from the Intramural Research Program of
% the National Institute of Mental Health (ZIAMH002886) and the Signal Processing and Instrumentation Section in the NIH
% Center for Information Technology. Please cite: Mitz, A. R., Chacko, R. V., Putnam, P. T., Rudebeck, P. H., & Murray,
% E. A. (2017). Using pupil size and heart rate to infer affective states during behavioral neurophysiology and
% neuropsychology experiments. Journal of neuroscience methods, 279, 1–12.
% https://doi.org/10.1016/j.jneumeth.2017.01.004
%
% Input:
%       EEG_eye             - EEGLAB dataset containing eye gaze data with pupil radius
%       idx_pupil           - indices of pupil channels (mean of all channels is used to detect blinks)
%       idx_clean           - OPTIONAL indices of channels to clean using pchip interpolation (default = all)
%       ms_search           - OPTIONAL searchbuffer to detect blinks in milliseconds (default = 20)
%       ms_apply            - OPTIONAL buffer to apply detected blinks in milliseconds (default = 30)
%       n                   - OPTIONAL standard deviation threshold used for blink detection (default = 3)
%       min_duration        - OPTIONAL after cleaning the eye data, blink events that lasted less than the minimal
%                               duration in seconds are removed (default = 0.1)
%       createplots         - OPTIONAL boolean whether plots should be created (default = 1)
%
% Output:
%       EEG_eye             - EEGLAB dataset containing cleaned eye data with blinkextract information stored in .etc
%       plothandles         - handle to the three created plots to allow script-based saving and closing 
%
% SEE ALSO:
%   blinkextract


function [EEG_eye, plothandles] = bemobil_clean_eye(EEG_eye,idx_pupil,idx_clean,ms_search,ms_apply,n,min_duration,createplots)

if nargin == 0
    help bemobil_clean_eye
    return
end

if ~isfield(EEG_eye,'data')
    error('Wrong data type for EEG_eye input! Needs to be an EEGLAB struct.')
end

if ~exist('idx_pupil','var') || ~isvector(idx_pupil)
    error('Missing required input "idx_pupil" of pupil radius indices.')
end

cleanallchan = 0;
if ~exist('idx_clean','var') || ~isvector(idx_clean)
    disp('Cleaning all channels.')
    idx_clean = 1:EEG_eye.nbchan;
    cleanallchan = 1;
end

if ~exist('ms_search','var') || ~isscalar(ms_search)
    ms_search = 20;
    disp(['Using default ' num2str(ms_search) ' ms search buffer to detect blinks.'])
end

if ~exist('ms_apply','var') || ~isscalar(ms_apply)
    ms_apply = 30;
    disp(['Using default ' num2str(ms_apply) ' ms apply buffer to remove blinks.'])
end

if ~exist('n','var') || ~isscalar(n)
    n = 3;
    disp(['Using ' num2str(n) ' standard deviations to detect blinks.'])
end

if ~exist('min_duration','var') || ~isscalar(min_duration)
    min_duration = 0.1;
end
    
if ~exist('createplots','var') || ~isscalar(createplots)
    createplots = 1;
end

plothandles = [];

%% plot raw

if createplots
    
    plothandles(1) = figure('color','w','position',[50 100 1500 800]);
    plot(transpose(normalize(EEG_eye.data,2,'range',[0 1.5]) + flipud([1:2:2*EEG_eye.nbchan]')), 'color', [78 165 216]/255)
    xlim([0 EEG_eye.pnts])
    ylim([0 2*EEG_eye.nbchan+1])
    yticks([1:2:2*EEG_eye.nbchan]+0.5)
    yticklabels(strrep(fliplr({EEG_eye.chanlocs.labels}),'_',' '))
    title('Raw data')
    drawnow
    
end

%% detect blinks

labels_pupil = {EEG_eye.chanlocs(idx_pupil).labels}';

disp(['Detecting blinks based on pupil radius of channel(s): [' num2str(idx_pupil) '], with label(s): ' strjoin(labels_pupil,', ')])
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')
disp('Mitz, A. R., Chacko, R. V., Putnam, P. T., Rudebeck, P. H., & Murray, E. A. (2017).')
disp('Using pupil size and heart rate to infer affective states during behavioral neurophysiology and neuropsychology experiments.')
disp('Journal of neuroscience methods, 279, 1–12. https://doi.org/10.1016/j.jneumeth.2017.01.004')
disp(' ')
disp('--------------------------- PLEASE CITE ----------------------------')
disp(' ')

[~,blinks,blink_onsets_offsets,threshold] = blinkextract(mean(EEG_eye.data(idx_pupil,:),1),round(ms_search*EEG_eye.srate/1000),...
    'applybuffer',round(ms_apply*EEG_eye.srate/1000),'n',n,'createplot',0);

% store in etc struct
EEG_eye.etc.blinkextract.blinks = blinks;
EEG_eye.etc.blinkextract.blink_onsets_offsets = blink_onsets_offsets;
EEG_eye.etc.blinkextract.threshold = threshold;
EEG_eye.etc.blinkextract.idx_pupil = idx_pupil;
EEG_eye.etc.blinkextract.idx_clean = idx_clean;
EEG_eye.etc.blinkextract.ms_search = ms_search;
EEG_eye.etc.blinkextract.ms_apply = ms_apply;
EEG_eye.etc.blinkextract.n = n;

%% interpolate blinks

PR = mean(EEG_eye.data(idx_pupil,:),1);

if ~cleanallchan
    labels_clean = {EEG_eye.chanlocs(idx_clean).labels}';
    disp(['Cleaning channel(s): [' num2str(idx_clean) '], with label(s): ' strjoin(labels_clean,', ') ' using pchip interpolation!'])
else
    disp('Cleaning all channels using pchip interpolation!')
end

EEG_eye.data(idx_clean,logical(blinks)) = deal(nan);
for idx = idx_clean
    
    EEG_eye.data(idx,:) = fillmissing(EEG_eye.data(idx,:),'pchip','SamplePoints',EEG_eye.times,'endvalues','nearest');
    
end

newPR = mean(EEG_eye.data(idx_pupil,:),1);

%% plot clean

if createplots
    
    plothandles(2) = figure('color','w','position',[50 100 1500 800]);
    
    ax1 = subplot(311);
    plot(newPR)
    title('Cleaned Pupil Radius')
    
    ax2 = subplot(312);
    plot(PR)
    hold on
    p = plot(min(newPR)+(blinks*(max(newPR)-min(newPR))),'r');
    title('Raw Pupil Radius')
    legend(p,'blinks')
    
    ax3 = subplot(313);
    plot(newPR)
    hold on
    plot(min(newPR)+(blinks*(max(newPR)-min(newPR))),'r')
    title('Cleaned Pupil Radius')
    
    linkaxes([ax1, ax2, ax3])
    drawnow
    
    
    plothandles(3) = figure('color','w','position',[50 100 1500 800]);
    plot(transpose(normalize(EEG_eye.data,2,'range',[0 1.5]) + flipud([1:2:2*EEG_eye.nbchan]')), 'color', [78 165 216]/255)
    xlim([0 EEG_eye.pnts])
    ylim([0 2*EEG_eye.nbchan+1])
    yticks([1:2:2*EEG_eye.nbchan]+0.5)
    yticklabels(strrep(fliplr({EEG_eye.chanlocs.labels}),'_',' '))
    title('Cleaned data')
    drawnow
end

%% add events in EEG_eye set

disp('Adding "blink:start" and "blink:stop" events to data set.')

for i_blink = 1:size(blink_onsets_offsets,1)
    
    i = numel(EEG_eye.event) + 1;
    EEG_eye.event(i).type = 'blink:start';
    EEG_eye.event(i).latency = blink_onsets_offsets(i_blink,1);
    EEG_eye.event(i).duration = 1/EEG_eye.srate;
    
    i = numel(EEG_eye.event) + 1;
    EEG_eye.event(i).type = 'blink:stop';
    EEG_eye.event(i).latency = blink_onsets_offsets(i_blink,2);
    EEG_eye.event(i).duration = 1/EEG_eye.srate;
end

EEG_eye = eeg_checkset(EEG_eye, 'eventconsistency');

% remove events if below a duration threshold
if min_duration > 0
    
    disp(['Removing blink events that lasted less than ' num2str(min_duration*1000) ' milliseconds.'])
    
    idx_remove=[];
    blink = 0;
    startlatency = 0;
    lastidx = 0;
    
    for i=1:length(EEG_eye.event)
        if strcmp(EEG_eye.event(i).type,'blink:start')
            blink = 1;
            startlatency = EEG_eye.event(i).latency;
            lastidx = i;
        end
        if blink && strcmp(EEG_eye.event(i).type,'blink:stop')
            if EEG_eye.event(i).latency - startlatency < min_duration*EEG_eye.srate
                idx_remove(end+1) = lastidx;
                idx_remove(end+1) = i;
            end
            blink = 0;
        end
    end
    EEG_eye.event(idx_remove) = [];
    EEG_eye = eeg_checkset(EEG_eye, 'eventconsistency');
end

disp('Blink detection and eye data cleaning done!')
