% bemobil_dipoles() - Searches and stores dipole xyz information for all clusters and their 
% centroids in the STUDY.
%   
% THIS FUNCTION RELIES HEAVILY ON A EEGLAB FUNCTION STD_CENTROID, ALL CREDITS FOR THE PREVIOUS WORK
% GO TO Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005!
%
% Usage:
%   >>  STUDY = bemobil_dipoles(STUDY,ALLEEG)
%
% Inputs:
%   STUDY                   - STUDY data set
%   ALLEEG                  - complete EEGLAB data set structure
%
% Outputs:
%   STUDY                   - STUDY data set (with computed cluster locations)
%
% See also: 
%   EEGLAB, std_centroid
% 
% Authors: Marius Klug, 2017, Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005

function STUDY = bemobil_dipoles(STUDY,ALLEEG)

    % compute for all clusters except parentclust and outlier clust
    clsind = 1:length(STUDY.cluster)-2; 
    
    for clust = 1:length(clsind)
        clear all_diplocs
        max_r = 0;
        len = length(STUDY.cluster(clsind(clust)).comps);
        tmppos = [ 0 0 0 ];
        tmpmom = [ 0 0 0 ];
        tmprv = 0;
        ndip = 0;
        for k = 1:len 
            fprintf('.');
            comp  = STUDY.cluster(clsind(clust)).comps(k);
            abset = STUDY.cluster(clsind(clust)).sets(1,k);
            if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ], 'Aborting compute centroid dipole');
               return;
            end
            if ~isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
                ndip   = ndip +1;
                posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
                momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
                if size(posxyz,1) == 2
                    if all(posxyz(2,:) == [ 0 0 0 ])
                        posxyz(2,:) = [];
                        momxyz(2,:) = [];
                    end;
                end;
                tmppos = tmppos + mean(posxyz,1);
                tmpmom = tmpmom + mean(momxyz,1);
                tmprv = tmprv + ALLEEG(abset).dipfit.model(comp).rv;
                if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       load('-mat', ALLEEG(abset).dipfit.hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
                   end
                end
                all_diplocs(k,:) = posxyz;
            end
        end
        centroid{clust}.dipole.posxyz =  tmppos/ndip;
        centroid{clust}.dipole.momxyz =  tmpmom/ndip;
        centroid{clust}.dipole.rv =  tmprv/ndip;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical') & (~isfield(ALLEEG(abset).dipfit, 'hdmfile')) %old dipfit
            centroid{clust}.dipole.maxr = max_r;
        end
        STUDY.cluster(clsind(clust)).dipole = centroid{clust}.dipole;
        STUDY.cluster(clsind(clust)).all_diplocs = all_diplocs;
    end
    fprintf('\n');
end

