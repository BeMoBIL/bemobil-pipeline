function [y,yy,nremove,scores]=nt_zapline_bemobil(x,fline,nremove,p,plotflag)
%[y,yy]=nt_zapline(x,fline,nremove,p,plotflag) - remove power line artifact
%
%  y: denoised data
%  yy: artifact
%  nremove: adaptive nremove result
%  score: score of the components on which the adaptive algorithm is based
%
%  x: data
%  fline: line frequency (normalized to sr)
%  nremove: number of components to remove [default: 1]
%  p: additional parameters:
%    p.nfft: size of FFT [default:1024]
%    p.nkeep: number of components to keep in DSS [default: all]
%    p.niterations: number of iterations for smoothing filter
%    p.fig1: figure to use for DSS score [default: 100]
%    p.fig2: figure to use for results [default: 101]
%    p.adaptiveNremove: use adaptive detection method of artifact scores for
%		removal instead of predefined nremove
%  plotflag: plot
%
%Examples:
%  nt_zapline(x,60/1000) 
%    apply to x, assuming line frequency=60Hz and sampling rate=1000Hz, plot results
%  nt_zapline(x,60/1000,4)
%    same, removing 4 line-dominated components 
%  p=[];p.nkeep=30; nt_zapline(x,60/1000,4,p);
%    same, truncating PCs beyond the 30th to avoid overfitting
%  [y,yy]=nt_zapline(x,60/1000)
%    return cleaned data in y, noise in yy, don't plot
%
% Marius Klug: Added support for adaptive detection of nremove using outlier detection (2020)

% NoiseTools
try, nt_greetings; catch, disp('You must download NoiseNools from http://audition.ens.fr/adc/NoiseTools/'); return; end

assert(nargin>=2, '!'); 
if nargin<3||isempty(nremove); nremove=1; end
if nargin<4; p=[]; end
if ~isfield(p,'nfft'); p.nfft=1024; end
if ~isfield(p,'nkeep'); p.nkeep=[]; end
if ~isfield(p,'niterations'); p.niterations=1; end
if ~isfield(p,'fig1'); p.fig1=100; end
if ~isfield(p, 'fig2'); p.fig2=101; end
if ~isfield(p, 'adaptiveNremove'); p.adaptiveNremove=1; end
if nargin<5||isempty(plotflag); plotflag=0; end

if isempty(x); error('!'); end
if nremove>=size(x,1); error('!'); end
if fline>1/2; error('fline should be less than Nyquist'); end
if size(x,1)<p.nfft; warning(['reducing nfft to ',num2str(size(x,1))]); p.nfft=size(x,1); end

if ~nargout || plotflag
    % print result and display spectra
    [y,yy,nremove,scores]=nt_zapline_bemobil(x,fline,nremove,p); %%% MK added outputs
    disp('proportion of non-DC power removed:');
    disp(nt_wpwr(x-y)/nt_wpwr(nt_demean(x)));
    
    figure(p.fig2); clf;    
    subplot 121
    [pxx,f]=nt_spect_plot(x/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
    divisor=sum(pxx);
    semilogy(f,abs(pxx)/divisor);
    legend('original'); legend boxoff
    set(gca,'ygrid','on','xgrid','on');
    xlabel('frequency (relative to line)');
    ylabel('relative power');
    yl1=get(gca,'ylim');
    hh=get(gca,'children');
    set(hh(1),'color','k')
    subplot 122
    [pxx,f]=nt_spect_plot(y/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
    semilogy(f,abs(pxx)/divisor);
    if nremove~=0 %%% MK
    hold on
    [pxx,f]=nt_spect_plot((x-y)/sqrt(mean(x(:).^2)),p.nfft,[],[],1/fline);
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
    subplot 121; ylim(yl); subplot 122; ylim(yl);
    end
      
    drawnow;
    return
end
if ~nargout
    clear y yy
    return
end

xx=nt_smooth(x,1/fline,p.niterations); % cancels line_frequency and harmonics, light lowpass
if isempty(p.nkeep); p.nkeep=size(x,2); end
xxxx=nt_pca(x-xx,[],p.nkeep); % reduce dimensionality to avoid overfitting

% DSS to isolate line components from residual:
nHarmonics=floor((1/2)/fline);
[c0,c1]=nt_bias_fft(xxxx,fline*(1:nHarmonics), p.nfft);

[todss,pwr0,pwr1]=nt_dss0(c0,c1);
scores = pwr1./pwr0; %%% MK
assert(size(todss,2)>0, '!'); 
if ~nargout;
    figure(p.fig1); clf;
    plot(pwr1./pwr0, '.-'); xlabel('component'); ylabel('score'); title('DSS to enhance line frequencies');
end

%% MK
if p.adaptiveNremove == 1
% 	nremove = triangle_threshold(pwr1./pwr0,'R',1); 
%   elbow detection does not work very well because the elbow is more
% 	significant for noisier datasets which menas that clean datasets get more components removed which defeats the
% 	purpose
    
    [adaptiveNremove, ~] = bemobil_iterative_threshold_detection(scores);
    
	
    if adaptiveNremove >= nremove
        fprintf('Adaptive score outlier detection found %d components, predefined nremove is %d. Removing %d components. This does not reduce the data rank!\n',adaptiveNremove,nremove,adaptiveNremove);
        nremove = adaptiveNremove;
    else
        fprintf('Adaptive score outlier detection found %d components, which is LESS than predefined nremove (%d). Using predefined removal to remove %d components. This does not reduce the data rank!\n',adaptiveNremove,nremove,nremove);
    end
    
end
%% MK end

xxxx=nt_mmat(xxxx,todss(:,1:nremove)); % line-dominated components
xxx=nt_tsr(x-xx,xxxx); % project them out
clear xxxx

% reconstruct clean signal
y=xx+xxx; clear xx xxx
yy=x-y;

% test code
if 0
    sr=400;
    nsamples=100000; nchans=100;
    signal=randn(nsamples,nchans);
    artifact=sin((1:nsamples)'/sr*2*pi*50);
    artifact=max(artifact,0).^3; % introduce harmonics
    artifact=3*nt_demean(artifact*randn(1,nchans));
    disp(nt_wpwr(artifact)/nt_wpwr(signal+artifact));
    nt_zapline_bemobil(signal+artifact,50/sr);
end

