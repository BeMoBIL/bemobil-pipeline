function elec = bemobil_create_elec(EEG)
% creates an elec matlab struct from a loaded EEG data set that can be stored as .mat file and used during import in the
% bemobil_xdf2bids function.
%
% example:
%   elec_struct_loaded            = load(['P:\situation_awareness\data\0_raw-data\s' num2str(subject) '\s' num2str(subject) '_elec_struct.mat']);
%   config.eeg.elec_struct        = elec_struct_loaded.elec;
%
% see also:
% bemobil_bids2set

if nargin == 0
    help bemobil_create_elec
    return
end

chanpos = [];
chanpos(:,1) = [EEG.chanlocs.X];
chanpos(:,2) = [EEG.chanlocs.Y];
chanpos(:,3) = [EEG.chanlocs.Z];

chantype = cell(size(chanpos,1),1);
chantype(1:size(chanpos,1)) = deal({'eeg'});

chanunit = cell(size(chanpos,1),1);
chanunit(1:size(chanpos,1)) = deal({'uV'});

elecpos = chanpos;

label = {EEG.chanlocs.labels}';

unit = 'mm';

elec.chanpos = chanpos;
elec.chantype = chantype;
elec.chanunit = chanunit;
elec.elecpos = elecpos;
elec.label = label;
elec.unit = unit;