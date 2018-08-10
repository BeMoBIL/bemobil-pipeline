function eeglab_basepath = bemobil_find_eeglab_basepath()

eeglab_path = which('eeglab');
eeglab_path_split = strsplit(eeglab_path,'\eeglab.m');
eeglab_basepath = eeglab_path_split{1};