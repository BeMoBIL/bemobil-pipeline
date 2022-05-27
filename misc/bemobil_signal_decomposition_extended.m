% bemobil_signal_decomposition_extended() - Computes a spatial filter for the EEG data set, which decomposes
% the data into components after optionally first using SSD to extract principle components (e.g. statistically independent components using AMICA - by default).
%
% SPoC and SSD combinations have been used but mSPoC and CCA are in
% experimental phase 
%
% Usage:
%   >>  [ALLEEG EEG CURRENTSET] = [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition_extended(ALLEEG, EEG, CURRENTSET, 'CCA','Z_target', epoch_zs, 'out_filename',out_filename,'out_filepath',out_filepath);
%   >>  [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition_extended(ALLEEG, EEG, CURRENTSET, 'SPOC_R2','Z_target', epoch_zs,'n_spoc_components',1, 'out_filename',out_filename,'out_filepath',out_filepath);
%
% Inputs:
%   ALLEEG                  - complete EEGLAB data set structure
%   EEG                     - current EEGLAB EEG structure
%   CURRENTSET              - index of current EEGLAB EEG structure within ALLEEG
%   algorithm               - Algorithm to use for extraction filters (no default, "AMICA", "CCA", "SSD", "MSPOC", "SPOC_LAMBDA", "SPOC_R2")
% 
% Optional Inputs: (Note not all inputs are optional for their respective algorithm. See algorithm header for details)
% Note: some algorithms take additional parameters not provided through the
% wrapper function. Use them alone for fine-grained param tuning. 
% 
%   preprocess_ssd          - Use SSD to reduce dimensionality of data (default false)
%   denoise_ssd             - Use SSD to denoise the data (default false)
%   ssd_n_components        - SSD: Num components to keep with dimn. reduction (default number channels)
%   num_models              - AMICA: number of models to learn (default 1)
%   max_threads             - AMICA: maximum of CPU threads to be used (default 1)
%   ssd_freq                - SSD: 3 x 2 matrix with the cut-off frequencies. (default [])
%                               First row: cut-off frequencies for band-pass of the to be extracted 
%                               oscillations.
%                               Second row: cut-off frequencies for the lowest and highest 
%                               frequencies defining flanking intervals.
%                               Third row: cut-off frequencies for the band-stop filtering of 
%                               the central frequency process.
%
%   ssd_sampling_freq       - SSD: the sampling frequency (in Hz) of X (default 250)
%   ssd_filter_order        - SSD: filter order used for butterworth bandpass and 
%                               bandstop filtering. If unsure about it, use [] then
%                               the default value of order = 2 will be used
%   ssd_epoch_indices      - SSD:  a matrix of size N x 2, where N is the number of 
%                               good (i.e. artifact free) data segments. Each row of 
%                               this matrix contains the first and the last sample index of
%                               a data segment, i.e. epoch_indices(n,:) = [1000, 5000]
%                               means that the n'th segment starts at sample 1000 and
%                               ends at sample 5000. If all data is useable or if unsure,
%                               use [], then the default of [1, size(X,1)] will be used.   
%    Z_target              - CCA/SPOC/MSPOC: The target function to be ran
%                               against. Check the respective algorithm for sizing requirments.
%
%    n_bootstrapping_iterations - SPOC: The number of repetitions to bootstrap the
%                                   distribution of lambda values with shuffled z.
%                                   Default value is 0, i.e. no bootstrapping and
%                                   therefore no p-values. 
%
%    pca_X_var_explained   - Optional dimensioality reduction with PCA, using
%                             the number of PCA components is determined by how
%                             much cumulative variance they should explain. This
%                             number must be between 0 and 1. Default is 1, i.e.
%                             no dimensionality reduction. (Not reccomended with SSD reduction)
%   n_spoc_components      - SPOC/MSPOC: The number of SPOC components to
%                             return, or in the case of MSPOC, the number of component pairs
% 
%   out_filename           - output filename 
%   out_filepath           - output filepath (File will only be saved on disk
%                               if both a name and a path are provided)
%
% Outputs:
%   ALLEEG                 - complete EEGLAB data set structure
%   EEG                    - current EEGLAB EEG structure
%   Currentset             - index of current EEGLAB EEG structure within ALLEEG
%
%   .set data file of current EEGLAB EEG structure stored on disk (OPTIONALLY)
%
% See also:
%    EEGLAB, runamica15
%
% Authors: Lukas Gehrke, Friederike Hohlefeld, Marius Klug, Luke Guerdan 2018

%%
function [ALLEEG EEG CURRENTSET] = bemobil_signal_decomposition_extended(ALLEEG, EEG, CURRENTSET, algorithm, varargin)


% Paramater checks
opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
                  'pca_X_var_explained', 1, ...
                  'numb_models', 1, ...
                  'maxx_threads', 1, ...
                  'preprocess_ssd', 0,...
                  'denoise_ssd', 0,...
                  'ssd_n_components', EEG.nbchan,...
                  'ssd_freq', [],...
                  'ssd_sampling_freq', 250,...
                  'ssd_filter_order', 2,...
                  'ssd_epoch_indices', [],...                 
                  'Z_target',[],...
                  'n_bootstrapping_iterations',0,...
                  'pca_X_var_explained',1,...
                  'n_spoc_components',EEG.nbchan,...
                  'out_filename', '',...
                  'out_filepath', '');
              
                                                                               
                                                                                                                          
 %%                                                                                                                      
% only save a file on disk if both a name and a path are provided
save_file_on_disk = ((~strcmp(opt.out_filename, '')) && (~strcmp(opt.out_filepath, '')));

% check if file already exist and show warning if it does
if save_file_on_disk
    mkdir(opt.out_filepath); % make sure that folder exists, nothing happens if so
    dir_files = dir(opt.out_filepath);
    if ismember(opt.out_filename, {dir_files.name})
        warning([opt.out_filename ' file already exists in: ' opt.out_filepath '. File will be overwritten...']);
    end
else
    % necessary for AMICA write-out, there won't be other data savings than the AMICA folder!
    opt.out_filepath = EEG.filepath;
    opt.out_filename = EEG.filename;
end

%UPDATE
if isfield(EEG,'datfile') && ~isempty(EEG.datfile)
     disp('Found datfile.');
     data = [EEG.filepath '\' EEG.datfile];        
else
     disp('No datfile field found in EEG structure. Will write temp file in current directory.');
     data = EEG.data(:,:);
end


% delete potentially preexistent folder since it will interfere in case algorithm crashes
if ~isempty(dir([opt.out_filepath '\' opt.out_filename '_' algorithm])); rmdir([opt.out_filepath '\' opt.out_filename '_' algorithm],'s'); end
EEG_data = EEG.data;


if((opt.preprocess_ssd == true && strcmp(algorithm, 'SSD')) ||...
   (opt.denoise_ssd == true && strcmp(algorithm, 'SSD')) ||...
   (opt.denoise_ssd == true && opt.preprocess_ssd))
   
    error('SSD settings incompatible')
end


if opt.preprocess_ssd
    EEG.etc.spoc.preprocess_ssd = true;
    [EEG_data, W_ssd, C_ssd] = ssd_preprocess(EEG_data, 'dim_reduction', opt.ssd_n_components, opt.ssd_freq, opt.ssd_sampling_freq, opt.ssd_filter_order, opt.ssd_epoch_indices);      
elseif opt.denoise_ssd
    EEG.etc.spoc.denoise_ssd = true;
    [EEG_data, W_ssd, C_ssd] = ssd_preprocess(EEG_data, 'deoising', opt.ssd_n_components, opt.ssd_freq, opt.ssd_sampling_freq, opt.ssd_filter_order, opt.ssd_epoch_indices);       
end
if opt.preprocess_ssd || opt.denoise_ssd
    EEG.etc.ssd.ssd_n_components = opt.ssd_n_components;
    EEG.etc.ssd.ssd_freq = opt.ssd_freq; 
    EEG.etc.ssd.W_ssd = W_ssd;
    EEG.etc.ssd.C_ssd = C_ssd;
end

%Check that the default or fed-in num components is correct (may not be after SSD preprocessing)
if strcmp(algorithm,'SPOC_LAMBDA') || strcmp(algorithm,'SPOC_R2') 
    
    if(ndims(EEG_data) < 3)
        error('SPoC requires epoched data.');
    end 
    data_rank = rank(EEG_data(:,:,1));
    if (data_rank < opt.n_spoc_components)
        disp(['Low data rank detected. Calculating ' num2str(data_rank) ' components'])
        opt.n_spoc_components =data_rank;
    end
end

disp(['Starting ' algorithm '...'])

switch algorithm
    
    case 'AMICA'
        %%Amica has not been tested in this function but has been included
        %%for completeness
        
        
        %%AMICA output interpretation currently not consistant with SPOC
        %%aka ICA will not work in this signal decomposition func yet. I
        %%need to fix implementation first
        while opt.maxx_threads > 0
            % try/catch loop because AMICA can crash dependent on the data set and the number of threads
            try
                [w, s, mods] = runamica15(data,...
                    'num_models', opt.numb_models,...
                    'max_threads', opt.maxx_threads,...
                    'outdir', [opt.out_filepath '\' opt.out_filename '_AMICA'],...
                    'num_chans', EEG.nbchan,...
                    'writestep', 2000);
                disp('AMICA successfull, storing weights and sphere.');
                EEG.etc.spatial_filter.AMICAmods = mods; %if this is used we need to replace this on all algos. Stored just for records or used? 
                EEG.icaweights = w;
                EEG.icasphere = s;
                % if successful, get out of the loop
                break 

            catch
                % if error, reduce threads by one
                opt.maxx_threads = opt.maxx_threads - 1;
                warning(['AMICA crashed. Reducing maximum threads to ' num2str(opt.maxx_threads)]);
            end
        end

        if opt.maxx_threads == 0
            warning('AMICA crashed with all possible maximum thread options. Try increasing the maximum usable threads of your CPU. If the maximum number of threads has already been tried, you''re pretty much fucked. Ask Jason Palmer, the creator of AMICA.');
            disp('Continuing with default EEGLAB runica() ...');
            
            [w, s] = runica(EEG.data);
            EEG.icaweights = w;
            EEG.icasphere = s;
            disp('runica successfull, storing weights and sphere.');
            algorithm = 'RUNICA';
        end
    
    case 'SSD' 
        
        if(ndims(EEG_data) > 2)
            error('SSD only accepts continuous datasets.');
        end 
        
        [W, A, lambda, C_s, X_ssd]  = ssd(EEG_data, opt.ssd_freq, opt.ssd_sampling_freq, opt.ssd_filter_order, opt.ssd_epoch_indices);       
        EEG.etc.spatial_filter.ssd.W = W;
        EEG.etc.spatial_filter.ssd.A = A;
        EEG.etc.spatial_filter.ssd.lamda = lambda;
        EEG.etc.spatial_filter.ssd.C_s = C_s;
        EEG.etc.spatial_filter.ssd.X_ssd = X_ssd;
        EEG.etc.spatial_filter.ssd.ssd_freq = opt.ssd_freq;
    
    case 'SPOC_LAMBDA'
        
        [W, A, lambda_values, p_values_lambda, Cxx, Cxxz, Cxxe] = spoc(EEG_data, opt.Z_target, 'n_spoc_components', opt.n_spoc_components, 'n_bootstrapping_iterations', opt.n_bootstrapping_iterations, 'pca_X_var_explained', opt.pca_X_var_explained);
        EEG.etc.spatial_filter.spoc.type = 'SPOC_LAMBDA';
        EEG.etc.spatial_filter.spoc.W = W;
        EEG.etc.spatial_filter.spoc.A = A;
        EEG.etc.spatial_filter.spoc.lambda_values = lambda_values;
        EEG.etc.spatial_filter.spoc.p_values_lamda = p_values_lambda;
        EEG.etc.spatial_filter.spoc.Cxx = Cxx;
        EEG.etc.spatial_filter.spoc.Cxxz = Cxxz;
        EEG.etc.spatial_filter.spoc.Cxxe = Cxxe;
        
    case 'SPOC_R2' 
        
        [W, A, Cxx, Cxxz, Cxxe, r2_values, r_values, lambda_values, ~] = spoc_r2(EEG_data, opt.Z_target, 'n_spoc_components', opt.n_spoc_components, 'pca_X_var_explained', opt.pca_X_var_explained);
        EEG.etc.spatial_filter.spoc.type = 'SPOC_R2';
        EEG.etc.spatial_filter.spoc.W = W;
        EEG.etc.spatial_filter.spoc.A = A;
        EEG.etc.spatial_filter.spoc.lambda_values = lambda_values;
        EEG.etc.spatial_filter.spoc.r2_values = r2_values;
        EEG.etc.spatial_filter.spoc.r_values = r_values;
        EEG.etc.spatial_filter.spoc.Cxx = Cxx;
        EEG.etc.spatial_filter.spoc.Cxxz = Cxxz;
        EEG.etc.spatial_filter.spoc.Cxxe = Cxxe;
        
    case 'MSPOC'
         opt.n_spoc_components = min(rank(EEG_data()))
        
         [W, Wy, Wtau, A, Ay, ~] = mspoc(EEG_data, opt.Z_target,'n_component_sets', opt.n_spoc_components);
         EEG.etc.spatial_filter.spoc.type = 'MSPOC';
         EEG.etc.spatial_filter.spoc.W = W;
         EEG.etc.spatial_filter.spoc.Wy = Wy;
         EEG.etc.spatial_filter.spoc.Wtau = Wtau;
         EEG.etc.spatial_filter.spoc.A = A;
         EEG.etc.spatial_filter.spoc.Ay = Ay;
         
    case 'CCA'
        
        %temporary, assumes epoched data, should run check later
        if(ndims(EEG_data) > 2)
            error('SSD only accepts continuous datasets.');
        end
             
        [W,~] = canoncorr(EEG_data,opt.Z_target); 
        %W are conncor weights maping EEG to optimially correlated subspace
        A = pinv(W');   
        EEG.etc.spatial_filter.cca.W = W;
        EEG.etc.spatial_filter.cca.A = A;
end


if not(strcmp(algorithm, 'AMICA'))

    % If SSD dim reduction is enabled spacial filters found are in SSD space. 
    % We need to map them back to the input data space to be useful.
    % SSD denoiseing doesn't require this back transformation
    if opt.preprocess_ssd
        W = W_ssd * W;
        EEG.icawinv = C_ssd * W;       
    else
        EEG.icawinv = A;
    end
%     EEG.spocweights = W; 
    
%     EEG = compute_filter_activations(EEG);
    %dummy values to pass EEG checkset
%     EEG.icasphere = zeros(size(EEG.icawinv.')); 
%     EEG.icaweights = zeros(size(EEG.icawinv,2));

    EEG.icachansind = 1:size(EEG.icawinv,1);    
	
	EEG.icaweights = eye(size(W,2));
	EEG.icasphere = W';
end


EEG.etc.spatial_filter.algorithm = algorithm;
EEG.etc.spatial_filter.original_data_path = [opt.out_filepath '\' opt.out_filename];

%Passing checkset criteria
if EEG.trials > 1 && ~isfield(EEG.event, 'epoch')
    EEG.trials = 1;
    EEG.pnts = size(EEG.data, 2)* size(EEG.data, 3);
    EEG.data = reshape(EEG.data, [size(EEG.data, 1),size(EEG.data, 2)* size(EEG.data, 3)]);
    EEG.icaact = reshape(EEG.icaact, [size(EEG.icaact,1), size(EEG.icaact,2)*size(EEG.icaact,3) ]);   
end

EEG = eeg_checkset(EEG); 

%[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
%EEG = eeg_checkset( EEG );

% save on disk
if save_file_on_disk
    EEG = pop_saveset( EEG, 'filename',opt.out_filename,'filepath', opt.out_filepath);
    disp('...done');
end

%[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);