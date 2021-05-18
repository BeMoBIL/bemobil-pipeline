% CLEAN_DATA_WITH_ZAPLINE - Removial of frequency artifacts using ZapLine to remove line noise from EEG/MEG data. Adds
% automatic detection of the number of components to remove, and chunks the data into segments to account for
% nonstationarities. Based on: de Cheveigne, A. (2020) ZapLine: a simple and effective method to remove power line
% artifacts. Neuroimage, 1, 1-13.
% 
% Requires noisetools to be installed: http://audition.ens.fr/adc/NoiseTools/
%
% Usage:
%   >>  [cleanData, resNremoveFinal, resScores, plothandles] = clean_data_with_zapline(data, srate, linefreqs, varargin);
%   
% Required Inputs:
%   data                    - EEG data matrix
%   srate                   - sampling rate in Hz
%   linefreqs               - vector with one or more line frequencies to be removed
% 
% Optional Parameters:
%   adaptiveNremove         - bool. if automatic adaptation of removal should be used. (default = 1)
%   fixedNremove            - numerical vector. fixed number of removed components. if adaptive removal is used, this 
%                               will be the minimum. can be either a scalar (then it is used for all line freqs) or a 
%                               vector of the same length as the linefreqs (then individual fixed n remove will be
%                               used). (default = 0)
%   chunkLength             - numerical. length of chunks to be cleaned in seconds. if set to 0, no chunks will be used. 
%                               (default = 30)
%   plotResults             - bool. if plot should be created. takes time to compute the spectrum. (default = 1)
%   figBase                 - integer. figure number to be created and plotted in. each iteration of linefreqs increases 
%                               this number by 1. (default = 100)
%   nfft                    - numerical. fft window size for computing the spectrum. (default = 512)
%   nkeep                   - integer. PCA reduction of components before removal. (default = round(20+size(data,2)/4))
%   initialSigma            - numerical. initial iterative outlier detection sigma threshold. (default = 3)
%   sigmaIncrease           - numerical. iterative outlier detection sigma threshold increase per iteration (to ensure 
%                               convergence). (default = 0.1)
% 
% Outputs:
%   cleanData               - clean EEG data matrix
%   zaplineNremoveFinal     - matrix of number of removed components per linefreq and chunk
%   scores                  - matrix of artifact component scores per linefreq and chunk
%   zaplineConfig           - config struct with all used parameters
%   plothandles             - vector of handles to the created figures
% 
% Example:
%   pop_eegplot( EEG, 1, 1, 1);
%   [cleanData, resNremoveFinal, resScores, plothandles] = clean_data_with_zapline(EEG.data,EEG.srate,[50 60]);
%   EEG.data = cleanData;
%   pop_eegplot( EEG, 1, 1, 1);
%
% See also:
%   nt_zapline_plus, iterative_outlier_removal
%
% Author: Marius Klug, 2021

function [cleanData, resNremoveFinal, resScores, zaplineConfig, plothandles] = clean_data_with_zapline(data, srate, linefreqs, varargin)

if nargin == 0
    help clean_data_with_zapline
    return
end

disp('Removing frequency artifacts using ZapLine with adaptations for automatic component selection and chunked data.')
disp('---------------- PLEASE CITE ------------------')
disp('de Cheveigne, A. (2020) ZapLine: a simple and effective method to remove power line artifacts. Neuroimage, 1, 1-13.')
disp('---------------- PLEASE CITE ------------------')

% input parsing settings
p = inputParser;
p.CaseSensitive = false;

addRequired(p, 'data', @(x) validateattributes(x,{'numeric'},{'matrix'},'clean_EEG_with_zapline','data'))
addRequired(p, 'srate', @(x) validateattributes(x,{'numeric'},{'positive','scalar','integer'},'clean_EEG_with_zapline','srate'))
addRequired(p, 'linefreqs', @(x) validateattributes(x,{'numeric'},{'positive','vector'},'clean_EEG_with_zapline','linefreqs'))
addOptional(p, 'adaptiveNremove', 1, @(x) validateattributes(x,{'numeric','logical'},{'scalar','binary'},'clean_EEG_with_zapline','adaptiveNremove'));
addOptional(p, 'fixedNremove', 0, @(x) validateattributes(x,{'numeric'},{'integer'},'clean_EEG_with_zapline','fixedNremove'));
addOptional(p, 'nfft', 512, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'},'clean_EEG_with_zapline','nfft'));
addOptional(p, 'initialSigma', 3, @(x) validateattributes(x,{'numeric'},{'scalar','positive'},'clean_EEG_with_zapline','initialSigma'));
addOptional(p, 'sigmaIncrease', 0.1, @(x) validateattributes(x,{'numeric'},{'scalar'},'clean_EEG_with_zapline','sigmaIncrease'));
addOptional(p, 'chunkLength', 30, @(x) validateattributes(x,{'numeric'},{'scalar','integer'},'clean_EEG_with_zapline','chunkLength'));
addOptional(p, 'nkeep', 0, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'},'clean_EEG_with_zapline','nkeep'));
addOptional(p, 'plotResults', 1, @(x) validateattributes(x,{'numeric','logical'},{'scalar','binary'},'clean_EEG_with_zapline','plotResults'));
addOptional(p, 'figBase', 100, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'},'clean_EEG_with_zapline','figBase'));

% parse the input
parse(p,data,srate,linefreqs,varargin{:});

data = p.Results.data;
srate = p.Results.srate;
linefreqs = p.Results.linefreqs;
adaptiveNremove = p.Results.adaptiveNremove;
fixedNremove = p.Results.fixedNremove;
nfft = p.Results.nfft;
initialSigma = p.Results.initialSigma;
sigmaIncrease = p.Results.sigmaIncrease;
chunkLength = p.Results.chunkLength;
nkeep = p.Results.nkeep;
plotResults = p.Results.plotResults;
figBase = p.Results.figBase;

% finalize
transposeData = size(data,2)>size(data,1);
if transposeData
    data = data';
end
if nkeep == 0
    nkeep = round(20+size(data,2)/4);
    disp(['Reducing the number of components to ' num2str(nkeep) ', set the ''nkeep'' flag to decide otherwise.'])
end
if chunkLength == 0
    chunkLength = size(data,1)/srate;
end

assert(isscalar(fixedNremove)||length(fixedNremove)==length(linefreqs),'''fixedNremove'' has to be either a scalar or the same length as linefreqs!')


zaplineConfig.adaptiveNremove = adaptiveNremove;
zaplineConfig.nkeep = nkeep;
zaplineConfig.linefreqs = linefreqs;
zaplineConfig.adaptiveNremove = adaptiveNremove;
zaplineConfig.fixedNremove = fixedNremove;
zaplineConfig.nfft = nfft;
zaplineConfig.initialSigma = initialSigma;
zaplineConfig.sigmaIncrease = sigmaIncrease;
zaplineConfig.chunkLength = chunkLength;
zaplineConfig.nkeep = nkeep;

%% Clean each frequency one after another
for iLinefreq = 1:length(linefreqs)
    
    linefreq = linefreqs(iLinefreq);
    
    if isscalar(fixedNremove)
        thisFixedNremove = fixedNremove;
    else
        thisFixedNremove = fixedNremove(iLinefreq);
    end
        
    
    % needs to be normalized for zapline
    fline = linefreq/srate;
    
    fprintf('Removing noise at %gHz... \n',linefreq);
    
    figThis = figBase+iLinefreq;
    
    % result data matrix
    cleanData = NaN(size(data));
    % last chunk must be larger than the others, to ensure fft works, but at least 1 chunk must be used
    nChunks = max(floor(size(data,1)/srate/chunkLength),1); 
    scores = NaN(nChunks,nkeep);
    NremoveFinal = NaN(nChunks,1);
    
    for iChunk = 1:nChunks
        if mod(iChunk,round(nChunks/10))==0
            disp(['Chunk ' num2str(iChunk) ' of ' num2str(nChunks)])
        end
        
        if iChunk ~= nChunks
            chunkIndices = 1+chunkLength*srate*(iChunk-1):chunkLength*srate*(iChunk);
        else
            chunkIndices = 1+chunkLength*srate*(iChunk-1):size(data,1);
        end
        
        chunk = data(chunkIndices,:);
        [cleanData(chunkIndices,:),~,NremoveFinal(iChunk),thisScores] =...
            nt_zapline_plus(chunk,fline,thisFixedNremove,zaplineConfig,0);
        scores(iChunk,1:length(thisScores)) = thisScores;
    end
    disp('Done.')
    
    if plotResults
        disp('Computing spectra for plotting.')
        
        % plot nremoved and scores        
        
        plothandles(iLinefreq) = figure(figThis);
        clf; set(gcf,'color','w','Position',[1 1 1920 1080])
        
        subplot(3,10,[1:8]);
        
        plot(NremoveFinal)
        xlabel('Chunk')
        ylabel('Removed Comps')
        title(['Removed ' num2str(linefreq) 'Hz components per ' num2str(chunkLength)...
            's chunk (automatic detection), M = ' num2str(mean(NremoveFinal))])
        
        figure(figThis);
        subplot(3,10,[9:10]);
        
        plot(nanmean(scores,1))
        hold on
        plot([mean(NremoveFinal) mean(NremoveFinal)],ylim,'r')
        title(['Mean artifact scores'])
        xlabel('Component')
        drawnow
        
        proportion_removed = nt_wpwr(data-cleanData)/nt_wpwr(nt_demean(data));
        disp(['proportion of non-DC power removed: ' num2str(proportion_removed)]);
        
        % plot starting spectrum
        
        figure(figThis);
        subplot(3,10,[11:15 21:25]);
        
        [pxx,f]=nt_spect_plot(data/sqrt(mean(data(:).^2)),nfft,[],[],1/fline);
        divisor=sum(pxx);
        
        figure(figThis);
        subplot(3,10,[11:15 21:25]);
        
        semilogy(f,abs(pxx)/divisor);
        legend('original'); legend boxoff
        set(gca,'ygrid','on','xgrid','on');
        xlabel('frequency (relative to line)');
        ylabel('relative power');
        yl1=get(gca,'ylim');
        hh=get(gca,'children');
        set(hh(1),'color','k')
        title(['Line frequency: ' num2str(linefreq) 'Hz'])
        drawnow
        
        % plot clean spectrum
        
        figure(figThis);
        subplot(3,10,[16:20 26:30]);
        
        [pxx,f]=nt_spect_plot(cleanData/sqrt(mean(data(:).^2)),nfft,[],[],1/fline);
        
        figure(figThis);
        subplot(3,10,[16:20 26:30]);
        
        semilogy(f,abs(pxx)/divisor);
        drawnow
        
        % plot removed noise spectrum
        hold on
        [pxx,f]=nt_spect_plot((data-cleanData)/sqrt(mean(data(:).^2)),nfft,[],[],1/fline);
        
        figure(figThis);
        subplot(3,10,[16:20 26:30]);
        hold on
        
        semilogy(f,abs(pxx)/divisor);
        legend('clean', 'removed'); legend boxoff
        set(gca,'ygrid','on','xgrid','on');
        set(gca,'yticklabel',[]); ylabel([]);
        xlabel('frequency (relative to line)');
        yl2=get(gca,'ylim');
        hh=get(gca,'children');
        set(hh(1),'color',[1 .5 .5]); set(hh(2), 'color', [ 0 .7 0]);
        set(hh(2),'linewidth', 2);
        yl(1)=min(yl1(1),yl2(1)); yl(2)=max(yl1(2),yl2(2));
        title(['proportion of non-DC power removed: ' num2str(proportion_removed)])
        
        % adjust scales
        ylim(yl);
        
        figure(figThis);
        subplot(3,10,[11:15 21:25]);
        
        ylim(yl);
        
        drawnow
        
        disp('...done.')
    end
    
    data = cleanData;
    resScores(iLinefreq,1:size(scores,1),1:size(scores,2)) = scores;
    resNremoveFinal(iLinefreq,1:size(NremoveFinal,1),1:size(NremoveFinal,2)) = NremoveFinal;
    
    % store in EEG file
%     EEG.data = y';
%     EEG.etc.zapline.nRemovedFinal(iLinefreq,:) = zaplineNremoveFinal;
%     EEG.etc.zapline.scores(iLinefreq,:,:) = scores;
end

if transposeData
    cleanData = cleanData';
end

if ~exist('plothandles','var')
    plothandles = [];
end