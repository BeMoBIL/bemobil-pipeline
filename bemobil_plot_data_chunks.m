function bemobil_plot_data_chunks(EEG, filepath, filename)
    plotfigure = figure('color','w');
    set(plotfigure, 'Position', get(0,'screensize'))
    ax1 = subplot(231);
    ax2 = subplot(232);
    ax3 = subplot(233);
    ax4 = subplot(234);
    ax5 = subplot(235);
    ax6 = subplot(236);

    starttime = EEG.times(end)/7*1;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax1,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 1 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax1);
    close(fighandle)

    starttime = EEG.times(end)/7*2;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax2,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 2 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax2);
    close(fighandle)

    starttime = EEG.times(end)/7*3;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax3,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 3 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax3);
    close(fighandle)

    starttime = EEG.times(end)/7*4;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax4,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 4 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax4);
    close(fighandle)

    starttime = EEG.times(end)/7*5;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax5,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 5 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax5);
    close(fighandle)

    starttime = EEG.times(end)/7*6;
    vis_artifacts(EEG,EEG,'show_events',1,'time_subset',...
        round([starttime starttime+10000]/1000)); % plot 10s at the first quarter
    axeshandle = gca;
    fighandle = gcf;
    axcp = copyobj(axeshandle, plotfigure);
    set(axcp,'Position',get(ax6,'position'));
    axcp.XTickLabel = [0:10]+round(starttime/1000);
    axcp.YTick=[];
    axcp.Title.String = ['Interpolated channels data section 6 of ' num2str(round(EEG.times(end)/1000)) 's'];
    axcp.XLabel.String = 'seconds';
    delete(ax6);
    close(fighandle)

    savefig(plotfigure,fullfile(filepath, [filename '.fig']))
    print(plotfigure,fullfile(filepath, [filename '.png'),'-dpng')
    close
end

