function [averages_plot, single_subject_plots] = bemobil_plot_ERSPs_channels(average_time_frequency_data,  all_subjects_time_frequency_data, scale)


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
    this_scale = [-this_scale/2 this_scale/2] ;
elseif length(scale)==1
    this_scale = [-this_scale*scale this_scale*scale];
elseif length(scale)==2
    this_scale = scale;
else
    error('Scale has wrong length. Must be 0 (autoscale to 1/2 of max), 1 (autoscale to scale*max), or 2 (completely manual scale of min and max bounds).')
end

if isempty(average_time_frequency_data.info.latencyMeans)
    disp('No means in latency present, no vertical lines for timewarp will be plotted.')
    do_vlines = false;
else
    do_vlines = true;
end

    do_dipoles = false;
    


% start the plot
averages_plot = figure;

% plot dipoles if possible
if do_dipoles
    subplot(n_subplots,1,1);
    std_dipplot(STUDY,ALLEEG,'clusters',average_time_frequency_data.info.cluster,'figure','off');
    set(averages_plot,'Color','white')
end
plot_alpha=false;
% plot average ERSPs
for ERSP = 1:n_ERSPs+1

        
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
        
        subplot(n_subplots,1,ERSP);
        
%         newfig
%         surf(average_time_frequency_data.times,...
%             average_time_frequency_data.freqs,...
%             data_to_plot_this_condition);
%         shading interp
%         red_blue_colormap
%         h =gca
%         this_scale
%         set(gca,'CLim',this_scale)
%         set(gca, 'YScale', 'log')
%         set(gca, 'YLim', [average_time_frequency_data.freqs(1) average_time_frequency_data.freqs(end)])
%         set(gca, 'YTick', [2 4 8 13 22 35 50 80])
        
        imagesclogy(average_time_frequency_data.times,...
            average_time_frequency_data.freqs,...
            data_to_plot_this_condition,...
            this_scale);%,...
        %             [average_time_frequency_data.times(1):200:average_time_frequency_data.times(end)]);
        set(gca,'XMinorTick','on')
        axis xy;
        
        %         xticks(0:200:xlim);
        
        % plot vlines of baseline in grey
        vertical_line(average_time_frequency_data.info.baseline_start_end,'white',2);
        % lines =average_time_frequency_data.info.baseline_start_end;
        % yLimits=get(gca,'ylim');             % Row vector
        % Ylimits=repmat(yLimits', 1, length(lines));
        % hold on
        %  handle = gca;
        %  plot(handle,[lines';lines'], Ylimits,'LineWidth',2);
        
        if do_vlines
            % plot vlines of timewarp in black
            vertical_line(average_time_frequency_data.info.latencyMeans,'black',2);
        end
        
        % @Marius: I had to change the order of arguments for strjoin, why
        % now, I dont know...
        % from: this_cond_title = strjoin(' ', strsplit(this_condition.condition_title,'_'));
        % to:
        this_cond_title = strjoin(strsplit(this_condition.condition_title,'_'), ' ');
        
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
        
    
end

ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');

set(get(ax,'Title'),'Visible','on')

plot_title = {['Channels: ' num2str(average_time_frequency_data.info.ICs_used) ', ' num2str(length(average_time_frequency_data.info.subjects_used))...
    ' subjects']};

if do_dipoles
    plot_title{end+1} = ['Mean residual variance: ' num2str(round(STUDY.cluster(average_time_frequency_data.info.cluster).mean_rv,3)*100) '%'];
end

if plot_alpha
    plot_title{end+1} = ['Significance level : ' num2str(average_time_frequency_data.alpha) ', ' num2str(average_time_frequency_data.n_permutes) ' permutes'];
end

plot_title{end+1} = ['Trial normalization : ' num2str(average_time_frequency_data.info.trial_normalization) ', auto epoch cleaning: ' num2str(average_time_frequency_data.info.do_auto_epoch_cleaning)];

plot_title{end+1} = ' ';
title(plot_title);

disp('done')


% plot single subject ERSPs
if ~isempty(all_subjects_time_frequency_data)
    
    disp('Plotting single subject ERSPs...')
    
    for subject = 1:length(all_subjects_time_frequency_data)
        
        single_subject_plots{subject} = figure;
        this_time_frequency_data = all_subjects_time_frequency_data(subject);
        
        % no difference plots here, but both conditions and the grand average
        for ERSP = 1:n_ERSPs
            try
                
                if ERSP == 1
                    
                    % if grand average plot
                    %                 fprintf('Plotting grand average data...\n')
                    this_condition = this_time_frequency_data.grand_average;
                    
                    plot_title = ['Subject: ' num2str(subject) ', Trial normalization : '...
                        num2str(average_time_frequency_data.info.trial_normalization) ', auto epoch cleaning: '...
                        num2str(average_time_frequency_data.info.do_auto_epoch_cleaning)...
                        ', ' this_condition.condition_title];
                    
                    
                else
                    
                    % if a condition plot
                    %                 fprintf('Plotting condition_%d data...\n', ERSP-1)
                    this_condition = this_time_frequency_data.(['condition_' num2str(ERSP-1)]);
                    
                    plot_title = [this_condition.condition_title];
                    
                end
                
                data_to_plot_this_condition = this_condition.ersp;
                
                
                subplot(n_ERSPs,1,ERSP);
                
                imagesclogy(this_time_frequency_data.times,...
                    this_time_frequency_data.freqs,...
                    data_to_plot_this_condition,...
                    this_scale);
                
                axis xy;
                
                if do_vlines
                    % plot vlines of timewarp in black
                    vline(average_time_frequency_data.info.latencyMeans,'black');
                end
                
                % plot vlines of baseline in grey
                vline(this_time_frequency_data.info.baseline_start_end,'white');
                
                title(plot_title);
                
                hold on
                
                % plot color scaling bar
                cbar;
                
            catch
                warning('Potato error message: U dun fucked up in plot ERSPs.')
            end
            
        end
        
        ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
        
        set(get(ax,'Title'),'Visible','on')
        
    end
    
    disp('...done')
else
    single_subject_plots = [];
end


function [lineHandles] = vertical_line(x_positions, lineType, lineWidth)
%
% vertical_line() - Draws vectical lines into the current plot. Linetype and linewidth can be specified.
%
% Usage:
%   >>  [lineHandles] = vertical_line(x_positions, lineType, lineWidth)
%
% Inputs:
%   x_positions     - positions on the x-axis where the lines should be drawn
%   lineType        - same as the ones used in plot()
%   lineWidth       - width of the line (best to use 1 or 2 at most)
%
% Outputs:
%   lineHandles     - handle to the drawn line
%
% Authors: Marius Klug, 2018
%
% Acknowledgement:  It was based on vline() written by Hoi Wong (2008).

if( ~isvector(x_positions) )
    error('x must be a vector');
else
    x_positions=x_positions(:);
end

if( ~exist('label', 'var') )        label = [];         end
if( ~exist('lineType', 'var') )     lineType = [];      end
if( ~exist('lineWidth', 'var') )    lineWidth = 1;      end

axesHandle = gca;

initial_ishold = ishold(axesHandle);
hold(axesHandle, 'on');

yLimits=get(axesHandle,'ylim');
Ylimits=repmat(yLimits', 1, length(x_positions));

if( isempty(lineType) )
    lineHandles = plot(axesHandle, [x_positions';x_positions'], Ylimits,'LineWidth',lineWidth);
else
    lineHandles = plot(axesHandle, [x_positions';x_positions'], Ylimits, lineType, 'LineWidth',lineWidth);
end

if ~initial_ishold
    hold(axesHandle, 'off');
end

% make it not show up in legends
set(lineHandles,'tag','vline','handlevisibility','off')
