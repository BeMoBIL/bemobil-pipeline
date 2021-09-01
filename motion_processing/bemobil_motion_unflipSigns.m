function EEG_motion_out = bemobil_motion_unflipSigns(EEG_motion_in)
% Heuristic for unflipping the sign of quaternion values to anable filtering. Quaternions can represent the
% same value two ways, whereas only the sign changes. Sometimes the representation flips in the time series.
% If this gets filtered, it creates an artifact, so this is why we want to unflip it first. It won't make a
% difference later when transforming the values to eul angles.
% The principle idea is to check if the difference between consecutive values becomes smaller if we flip them.
%
% Input arguments:
%       EEG dataset containing motion channels with quaternion data (must have exactly 4 channels with labels
%       containing "quat_x/y/z/w")
%
% Output argument:
%       EEG dataset with unflipped quaternion data
%
% Usage:


for channel = 1:EEG_motion_in.nbchan
    
    % checking for already present euls
    if any(~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'eul')))
        error('You can only unflip Quaternions, this dataset contains eul angles, try it with the original data set.')
    end
    
end


data = EEG_motion_in.data';

% find correct channelindex for the quaternion values of
% this RB
quaternionX = ~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat_x'));
quaternionY = ~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat_y'));
quaternionZ = ~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat_z'));
quaternionW = ~cellfun(@isempty,strfind(lower({EEG_motion_in.chanlocs.labels}),'quat_w'));

% take the values
X = data(:,quaternionX);
Y = data(:,quaternionY);
Z = data(:,quaternionZ);
W = data(:,quaternionW);

% unflip loop
for dataPoint = 5:size(data,1)
    
    epsilon = 0.5;
    
    Adiff = abs(X(dataPoint-1) - X(dataPoint));
    Adiff2 = abs(X(dataPoint-2) - X(dataPoint-1));
    Adiff3 = abs(X(dataPoint-3) - X(dataPoint-2));
    Adiff4 = abs(X(dataPoint-4) - X(dataPoint-3));
    AsumOfPreviousDiffs = Adiff2 + Adiff3 + Adiff4;
    
    AflippedDataPoint = -X(dataPoint);
    
    AdiffFlipped = abs(X(dataPoint-1) - AflippedDataPoint);
    
    AconditionMet = AdiffFlipped<Adiff;% & Adiff > AsumOfPreviousDiffs;
    AbigJump = AdiffFlipped> epsilon;
    
    Bdiff = abs(Y(dataPoint-1) - Y(dataPoint));
    Bdiff2 = abs(Y(dataPoint-2) - Y(dataPoint-1));
    Bdiff3 = abs(Y(dataPoint-3) - Y(dataPoint-2));
    Bdiff4 = abs(Y(dataPoint-4) - Y(dataPoint-3));
    BsumOfPreviousDiffs = Bdiff2 + Bdiff3 + Bdiff4;
    
    BflippedDataPoint = -Y(dataPoint);
    
    BdiffFlipped = abs(Y(dataPoint-1) - BflippedDataPoint);
    
    BconditionMet = BdiffFlipped<Bdiff;% & Bdiff > BsumOfPreviousDiffs;
    BbigJump = BdiffFlipped> epsilon;
    
    Cdiff = abs(Z(dataPoint-1) - Z(dataPoint));
    Cdiff2 = abs(Z(dataPoint-2) - Z(dataPoint-1));
    Cdiff3 = abs(Z(dataPoint-3) - Z(dataPoint-2));
    Cdiff4 = abs(Z(dataPoint-4) - Z(dataPoint-3));
    CsumOfPreviousDiffs = Cdiff2 + Cdiff3 + Cdiff4;
    
    CflippedDataPoint = -Z(dataPoint);
    
    CdiffFlipped = abs(Z(dataPoint-1) - CflippedDataPoint);
    
    CconditionMet = CdiffFlipped<Cdiff;% & Cdiff > CsumOfPreviousDiffs;
    CbigJump = CdiffFlipped> epsilon;
    
    Ddiff = abs(W(dataPoint-1) - W(dataPoint));
    Ddiff2 = abs(W(dataPoint-2) - W(dataPoint-1));
    Ddiff3 = abs(W(dataPoint-3) - W(dataPoint-2));
    Ddiff4 = abs(W(dataPoint-4) - W(dataPoint-3));
    DsumOfPreviousDiffs = Ddiff2 + Ddiff3 + Ddiff4;
    
    DflippedDataPoint = -W(dataPoint);
    
    DdiffFlipped = abs(W(dataPoint-1) - DflippedDataPoint);
    
    DconditionMet = DdiffFlipped<Ddiff;% & Ddiff > DsumOfPreviousDiffs;
    DbigJump = DdiffFlipped> epsilon;
    
    
    
    if (AconditionMet+BconditionMet+CconditionMet+DconditionMet >= 2 && AbigJump+BbigJump+CbigJump+DbigJump == 0)
        
        X(dataPoint) = -X(dataPoint);
        Y(dataPoint) = -Y(dataPoint);
        Z(dataPoint) = -Z(dataPoint);
        W(dataPoint) = -W(dataPoint);
        
    end
    
end

data(:,quaternionX) = X;
data(:,quaternionY) = Y;
data(:,quaternionZ) = Z;
data(:,quaternionW) = W;

%%
EEG_motion_out = EEG_motion_in;

EEG_motion_out.data = data';

end