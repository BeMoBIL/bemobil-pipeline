% evaluating the decomposition quality with SASICA

rejects = EEG.reject.gcompreject;
keeps = 1 - rejects;
numberBrain = sum(keeps);
autoCorr = sum(EEG.reject.SASICA.icaautocorr(logical(keeps)));
focal = sum(EEG.reject.SASICA.icafocalcomp(logical(keeps)));
scaledAutoCorr = autoCorr / numberBrain;
scaledFocal = focal / numberBrain;
final = scaledAutoCorr * scaledFocal

