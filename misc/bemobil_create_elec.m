function elec = bemobil_create_elec(EEG)

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