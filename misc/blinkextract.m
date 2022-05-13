function [newPR,blinks,blink_onsets_offsets,threshold,plot_handle] = blinkextract(PR,searchbuffer,varargin)
% BLINKEXTRACT returns the major and minor axis PS with blinks removed and returned in b.
%
%    V 1.0   17 July 2009    Ravi Chacko
%    V 1.1   24 Aug  2009    Andy Mitz, documentation only
%    V 1.6   25 Aug  2009    More detailed reporting during processing. Some error trapping
%    V 2.0    7 Oct  2009    Explicitly calls findpeaks2007b.m
%    V 2.7    9 Nov  2009    bl was a row vector. Changed now into a column vector to match PR
%	 V 2.8	  5 Nov  2019    Marius Klug, changes for NeSitA Project, bemobil.bpn.tu-berlin.de
%    V 3.0   20 Apr  2021    Marius Klug, added warning when no blinks were found
%
% INPUTS
%   PR				original Pupil Radius
%	searchbuffer	buffer around coarse theshold based blink detection to find more fine grained blinks (in samples)
%   n				OPTIONAL: number of standard deviations away from the mean (to establish blink threshold)
%					3 is a sufficient n (default = 3).
%	applybuffer		OPTIONAL: buffer around actually detected peaks in
%					samples (default = searchbuffer * 2/3)
%	createplot		OPTIONAL: boolean whether or not a plot of the data
%					should be created (default = 1)
%
% OUTPUTS
%   newPR					array		Pupil Radius (size) without spikes (blinks removed)
%   blinks					array		b=1 for applybuffer samples around detected blinks, 0 elsewhere
%	blink_onsets_offsets	2-D matrix	indices of blink start and end
%										points
%   threshold				scalar		threshold that was used for peak detection
%	plot_handle				handle		handle to the results plot if chosen

p = inputParser;

p.addRequired('PR',@(x) validateattributes(x,{'numeric'},{'vector'},'blinkextract_2','PR'));
p.addRequired('searchbuffer',@(x) validateattributes(x,{'numeric'},{'positive','scalar','integer'},'blinkextract_2','searchbuffer'));
p.addOptional('n',3,@(x) validateattributes(x,{'numeric'},{'positive','scalar'},'blinkextract_2','n'));
p.addOptional('applybuffer',0,@(x) validateattributes(x,{'numeric'},{'positive','scalar','integer'},'blinkextract_2','applybuffer'));
p.addOptional('createplot',1,@(x) validateattributes(x,{'numeric'},{'scalar','binary'},'createplot','applybuffer'));

p.parse(PR,searchbuffer,varargin{:});

n = p.Results.n;
applybuffer = p.Results.applybuffer;
createplot = p.Results.createplot;

if applybuffer == 0
    applybuffer = round(searchbuffer*(2/3));
end

%%
fprintf('Searching for blinks ');

% MK: added "nan" to allow for more robust estimation also when tracking was lost or subjects blinked very much.
PR(PR<0.2) = NaN;
m=nanmean(PR);
s=nanstd(PR);
threshold = (m-n*s); %threshold n std devs below the mean, line added by MK
PR(isnan(PR)) = 0; % Enable searching again

counter=2; % start with second sample since first sample should not at all be a blink
LPR=length(PR);
blinks = zeros(LPR,1);
progress = 10;
fprintf('0');
while counter < LPR % start with second sample since first sample should not at all be a blink
    if counter/LPR*100 >= progress
        fprintf('..%d',progress);
        progress=progress+10;
    end
    if PR(counter) < threshold
        % 		PR(counter)
        lbound = counter;
        while PR(counter) < (m-n*s) && counter < LPR
            counter=counter+1;
        end
        rbound = counter;
        [left,right] = link(lbound,rbound,searchbuffer,applybuffer,PR,LPR,threshold); % cuts out spike, see nested function link below
        blinks(left:right)=1;
        counter=right; 
    else
        counter = counter+1; % MK: changed +2 to +1 to not skip any samples
    end
end % while
fprintf('..100\n');

%% now actually clean the data
% this can only be done afterwards because only then we know consecutive
% detected blinks

disp('Interpolating blinks...')
newPR = PR;
eyeblink_data = diff(blinks);
blink_onsets_offsets = find(eyeblink_data~=0);

if length(blink_onsets_offsets) == 0
    warning('No blinks were found - are you sure the correct pupil radius data was used?')
    return
end

% test if we actually start and end with eyes open
if eyeblink_data(blink_onsets_offsets(1)) ~= 1
    blink_onsets_offsets(1) = [];
end
if eyeblink_data(blink_onsets_offsets(end)) ~= -1
    blink_onsets_offsets(end) = [];
end

blink_onsets_offsets = reshape(blink_onsets_offsets,[2,length(blink_onsets_offsets)/2])';

for i_blink = 1:length(blink_onsets_offsets)
    
    newPR(blink_onsets_offsets(i_blink,1)+1:blink_onsets_offsets(i_blink,2)) = ...
        linspace(newPR(blink_onsets_offsets(i_blink,1)),...
        newPR(blink_onsets_offsets(i_blink,2)+1),...
        blink_onsets_offsets(i_blink,2)-blink_onsets_offsets(i_blink,1));
end

%% plot
if createplot
	
	plot_handle = figure('color','w');
	
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
end

%% ========= Nested subfunction link() =========
% link the ends before and after blink
    function [left,right] = link(l,r,searchbuffer,applybuffer,PR,LPR,threshold)
		% MK: changed function to not use persistent variables but proper
		% input and output parameters
        
        if l-searchbuffer < 1   % Clip values that go below sample number 1
            lb = 1;
        else
            lb = l-searchbuffer;
        end
        
        if r+searchbuffer > LPR  % Clip values that go above the end of the data
            rb = LPR;
        else
            rb = r+searchbuffer;
        end
        
        diblink = diff(PR(lb:rb));
        [pks,locs] = findpeaks(abs(diblink)); % MK: use current version of findpeaks
        
        % MK changes to find actually blink start and end, not just peaks of
        % pupil size
        
        left = max(1,l-applybuffer); % left always starts at the applybuffer before threshold was reached, can't be lower than the first sample
        
        % right is after the last blink in the searchbuffer ends, plus
        % applybuffer
        pk_threshold = (max(pks)-min(pks))/3; % MK: we can ignore small pks, the large are actual blinks, the small are other pupil dilation
        
        locs(pks<pk_threshold) = [];
        
        if isempty(locs)==0
            right = l-searchbuffer+locs(length(locs));
            right = right+applybuffer; % MK added to have a buffer around blinks again after finding exact blink end
            if right > LPR    % to stay below upperbound
                right = LPR;
            end
            
            % MK: make sure it always at least stays within the [l,r]
            % boundaries, otherwise it can lead to issues when a blink is
            % happening at the end of the data
            right = max(r,right);
		else
			% if no peaks are found use default applybuffer
            right = max(LPR,r+applybuffer);
        end
        
        % MK: added to ensure that if the search starts or ends during a
        % blink it is correctly identified
        if PR(rb)<threshold
            right = min(LPR,rb+applybuffer);
        end
        
    end % end of nested subfunction link()

end % end of blinksearchbufferract()
