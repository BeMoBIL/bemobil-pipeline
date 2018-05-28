function plot = bemobil_plot_ERSPs(average_time_frequency_data,  all_subjects_time_frequency_data, STUDY, ALLEEG, scale)

% find residual variance info for title later
if ~isfield(STUDY.cluster, 'mean_rv')
    STUDY = bemobil_dipoles(STUDY,ALLEEG);
end



% find maximum or all ERSPs to plot

% grand average maximum and number of subplots
this_scale = max(max(abs(average_time_frequency_data.grand_average.ersp)));
n_subplots = 1;
n_ERSPs = 1;

% find condition maxima and if greater, use them
if  ~isempty(average_time_frequency_data.info.experiment_conditions)
    n_subplots = n_subplots+1; % difference plot
    
    for condition = 1:length(average_time_frequency_data.info.experiment_conditions)
        n_subplots = n_subplots+1;
        n_ERSPs = n_ERSPs+1;
        condmax = max(max(abs(average_time_frequency_data.(['condition_' num2str(condition)]).ersp)));
        if condmax > this_scale
            this_scale = condmax;
        end
    end
end

% apparently this is the standard way to go, EEGLAB at least does it so:
if isempty(scale) 
    this_scale = this_scale/2;
else
    this_scale = scale;
end

if isempty(average_time_frequency_data.info.latencyMeans)
    disp('No means in latency present, no vertical lines for timewarp will be plotted.')
    do_vlines = false;
else
    do_vlines = true;
end

if isempty(STUDY) || isempty(ALLEEG)
    disp('STUDY or ALLEEG data missing, no dipoles will be plotted.')
    do_dipoles = false;
else
    do_dipoles = true;
    n_subplots = n_subplots+1;
end


% start the plot
plot = figure;

% plot dipoles if possible
if do_dipoles
    subplot(n_subplots,1,1);
    std_dipplot(STUDY,ALLEEG,'clusters',average_time_frequency_data.info.cluster,'figure','off');
    set(plot,'Color','white')
end
plot_alpha=false;
% plot average ERSPs
for ERSP = 1:n_ERSPs+1
    try
        
        if ERSP == 1
            
            % if grand average plot
            fprintf('Plotting grand average data...\n')
            this_condition = average_time_frequency_data.grand_average;
            
        elseif ERSP == n_ERSPs+1
            
            % last plot is difference plot
            fprintf('Plotting difference data...\n')
            this_condition = average_time_frequency_data.difference;
            
        else
            
            % if a condition plot
            fprintf('Plotting condition_%d data...\n', ERSP-1)
            this_condition = average_time_frequency_data.(['condition_' num2str(ERSP-1)]);
        end
        
        data_to_plot_this_condition = this_condition.ersp;
        
        
        if ~isempty(this_condition.statistics)
            calculated_statistics = true;
            plot_alpha = true;
            masked_data_to_plot_this_condition = data_to_plot_this_condition.*this_condition.statistics.p_values_mask;
        else
            disp('No statistics data present in this condition. No contour lines for masking will be plotted.')
            calculated_statistics = false;
        end
        
        subplot(n_subplots,1,ERSP+1);
        
        imagesclogy(average_time_frequency_data.times,...
            average_time_frequency_data.freqs,...
            data_to_plot_this_condition,...
            [-this_scale this_scale]);%,...
%             [average_time_frequency_data.times(1):200:average_time_frequency_data.times(end)]);
        set(gca,'XMinorTick','on')
        axis xy;
        
%         xticks(0:200:xlim);
        
        % plot vlines of baseline in grey
        vline(average_time_frequency_data.info.baseline_start_end,'white',2);
% lines =average_time_frequency_data.info.baseline_start_end;
% yLimits=get(gca,'ylim');             % Row vector
% Ylimits=repmat(yLimits', 1, length(lines));
% hold on
%  handle = gca;
%  plot(handle,[lines';lines'], Ylimits,'LineWidth',2);
        
        if do_vlines
            % plot vlines of timewarp in black
            vline(average_time_frequency_data.info.latencyMeans,'black',2);
        end
        
        this_cond_title = strjoin(' ',strsplit(this_condition.condition_title,'_'));
        
        if isfield(this_condition,'n_epochs') && isfield(this_condition,'n_epochs_mean')
            this_cond_title = [this_cond_title ', ' num2str(this_condition.n_epochs) ' epochs, mean: ' num2str(this_condition.n_epochs_mean)];
        end
        
        if isfield(this_condition,'subjects')
            this_cond_title = [this_cond_title ', ' num2str(length(this_condition.subjects)) ' subjects'];
            
        end
        
        title(this_cond_title)
        
        hold on
        
        if calculated_statistics
            contour(average_time_frequency_data.times,...
                average_time_frequency_data.freqs,...
                logical(masked_data_to_plot_this_condition),...
                1,'linecolor','black', 'LineWidth',1);
        end
        
        % plot color scaling bar
        cbar;
        
    catch
        warning('Potato error message: U dun fucked up in plot ERSPs.')
    end
    
end

ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');

set(get(ax,'Title'),'Visible','on')

plot_title = {['Cluster: ' num2str(average_time_frequency_data.info.cluster) ', ' num2str(length(average_time_frequency_data.info.subjects_used))...
    ' subjects, ' num2str(length(average_time_frequency_data.info.ICs_used)) ' ICs'] };

if do_dipoles
    plot_title{end+1} = ['Mean residual variance: ' num2str(round(STUDY.cluster(average_time_frequency_data.info.cluster).mean_rv,3)*100) '%'];
end

if plot_alpha
    plot_title{end+1} = ['Significance level : ' num2str(average_time_frequency_data.alpha) ', ' num2str(average_time_frequency_data.n_permutes) ' permutes'];
end

plot_title{end+1} = ['Trial normalization : ' num2str(average_time_frequency_data.info.trial_normalization) ', auto epoch cleaning: ' num2str(average_time_frequency_data.info.do_auto_epoch_cleaning)];

title(plot_title);

disp('done')


% plot single subject ERSPs
% disp('Plotting single subject ERSPs...')
%
% for subject = 1:length(all_subjects_time_frequency_data)
%
%
%
%     for ERSP = 1:n_ERSPs+1
%         try
%
%             if ERSP == 1
%
%                 % if grand average plot
%                 fprintf('Plotting grand average data...\n')
%                 this_condition = average_time_frequency_data.grand_average;
%
%             elseif ERSP == n_ERSPs+1
%
%                 % last plot is difference plot
%                 fprintf('Plotting difference data...\n')
%                 this_condition = average_time_frequency_data.difference;
%
%             else
%
%                 % if a condition plot
%                 fprintf('Plotting condition_%d data...\n', ERSP-1)
%                 this_condition = average_time_frequency_data.(['condition_' num2str(ERSP-1)]);
%             end
%
%             data_to_plot_this_condition = this_condition.ersp;
%
%
%             if ~isempty(this_condition.statistics)
%                 error('statistics not implemented yet')
%                 calculated_statistics = true;
%                 masked_data_to_plot_this_condition = data_to_plot_this_condition.*squeeze(p_values_mask(cluster,plot,:,:));
%             else
%                 disp('No statistics data present in this condition. No contour lines for masking will be plotted.')
%                 calculated_statistics = false;
%             end
%
%             subplot(n_subplots,1,ERSP+1);
%
%             imagesclogy(average_time_frequency_data.times,...
%                 average_time_frequency_data.freqs,...
%                 data_to_plot_this_condition,...
%                 [-scale scale]);
%
%             axis xy;
%
%             if do_vlines
%                 % plot vlines of timewarp in black
%                 vline(average_time_frequency_data.info.latencyMeans,'black');
%             end
%
%             % plot vlines of baseline in grey
%             vline(average_time_frequency_data.info.baseline_start_end,'white');
%
%             if isfield(this_condition,'n_epochs') && isfield(this_condition,'n_epochs_means')
%                 title([strjoin(' ',strsplit(this_condition.condition_title,'_')) ', ' num2str(this_condition.n_epochs) ' epochs, mean: ' num2str(this_condition.n_epochs_means)]);
%             else
%                 title(strjoin(' ',strsplit(this_condition.condition_title,'_')));
%             end
%             hold on
%
%             if calculated_statistics
%                 contour(times,...
%                     freqs,...
%                     logical(masked_data_to_plot_this_condition),...
%                     1,'linecolor','black');
%             end
%
%             % plot color scaling bar
%             cbar;
%
%         catch
%             warning('Potato error message: U dun fucked up in plot ERSPs.')
%         end
%
%     end
%
%     ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
%
%     set(get(ax,'Title'),'Visible','on')
%
%     plot_title = {['Cluster: ' num2str(average_time_frequency_data.info.cluster) ', ' num2str(length(average_time_frequency_data.info.subjects_used))...
%         ' subjects, ' num2str(length(average_time_frequency_data.info.ICs_used)) ' ICs'] };
%
%     if do_dipoles
%         plot_title{end+1} = ['Mean RV: ' num2str(round(STUDY.cluster(average_time_frequency_data.info.cluster).mean_rv,3)*100) '%'];
%     end
%
%     if calculated_statistics
%         error('nope')
%         plot_title{end+1} = ['Significance level : ' num2str(alpha) ];
%     end
%
%     title(plot_title);
%
%     disp('done')
%
% end
%
% disp('...done')


function [lineHandles textHandles] = vline(x, lineType, lineWidth, label, textPosition, axesHandle, textColor)

% Draws vectical lines (specified by a vector of x)
%
% Except the first argument 'y', everything else are optional.
% To preserve input argument order, enter [] for unused input arguments
%
% lineType:     same as the ones used in plot(...)
% label:        simple strings applies to all lines
%               cell strings applies to each line
% textPosition: [x, y] bewteen 0 and 1 relative to the 'lower' end
%
% Author:      Hoi Wong (wonghoi.ee@gmail.com)
% Date:        02/14/2008
%
% Acknowledgement:  It was based on vline() written by Brandon Kuczenski.

if( ~isvector(x) )
    error('x must be a vector');
else
    x=x(:);     % Make sure it's column vector
end

if( ~exist('label', 'var') )        label = [];         end
if( ~exist('lineType', 'var') )     lineType = [];      end
if( ~exist('lineWidth', 'var') )    lineWidth = 1;      end
if( ~exist('axesHandle', 'var') )   axesHandle = gca;   end
if( ~exist('textPosition', 'var') ) textPosition = [];  end
if( ~exist('textColor', 'var') )    textColor = [];     end

if( isempty(axesHandle) )           axesHandle = gca;   end

if( isempty(textPosition) )
    textPositionX = 0.02;
    textPositionY = 0.02;
elseif( isscalar(textPosition) )
    textPositionX = 0.02;
    textPositionY = textPosition;
elseif( length(textPosition)~=2 )
    error('Invalid textPosition');
else
    textPositionX = textPosition(1);
    textPositionY = textPosition(2);
end
clear textPosition;     % Avoid typos for similarly named variables

holdState = ishold(axesHandle);
hold(axesHandle, 'on');

yLimits=get(axesHandle,'ylim');             % Row vector
Ylimits=repmat(yLimits', 1, length(x));

% Example: for horizontal lines
% X = [2 5        Y = [3 3
%      2 5];           4 4];

if( isempty(lineType) )
    lineHandles = plot(axesHandle, [x';x'], Ylimits,'LineWidth',lineWidth);
else
    lineHandles = plot(axesHandle, [x';x'], Ylimits, lineType, 'LineWidth',lineWidth);
end

if( ~isempty(label) )
    xLimits = get(axesHandle,'xlim');
    xLowerLimit = xLimits(2);
    xUpperLimit = xLimits(1);
    xRange      = xUpperLimit - xLowerLimit;
    xPosition   = x+textPositionX*xRange;
    
    yUpperLimit = yLimits(2);
    yLowerLimit = yLimits(1);
    yRange      = yUpperLimit - yLowerLimit;
    yPosition   = yLowerLimit+textPositionY*yRange;
    Yposition   = repmat(yPosition, length(x), 1);
    
    % textHandle is a vector correspond to the line
    textHandles = text(xPosition, Yposition, label, 'Parent', axesHandle);
    
    for k=1:length(x)
        % Set the text colors to be identical to line colors
        if isempty(textColor)
            set( textHandles(k), 'color', get(lineHandles(k), 'color') );
        else
            set( textHandles(k), 'color', textColor );
        end
    end
end

if( holdState==false )
    hold(axesHandle, 'off');
end
% this last part is so that it doesn't show up on legends
set(lineHandles,'tag','vline','handlevisibility','off')
