function bemobil_plot_amica_autoreject(EEG, filepath, filename)
    % plot autorejection
    data2plot = EEG.data(1:round(EEG.nbchan/10):EEG.nbchan,:)';
    figure;
    set(gcf,'color','w','Position', get(0,'screensize'));
    plot(data2plot,'g');
    data2plot(~EEG.etc.bad_samples,:) = NaN;
    hold on
    plot(data2plot,'r');
    xlim([-10000 EEG.pnts+10000])
    ylim([-1000 1000])
    title(['AMICA autorejection, removed ' num2str(round(EEG.etc.bad_samples_percent,2)) '% of the samples'])
    xlabel('Samples')
    ylabel('\muV')
    clear data2plot
    % save figure to disk
    savefig(gcf,fullfile(filepath,[filename '.fig']))
    print(gcf,fullfile(filepath,[filename '.png']),'-dpng')
    close
end

