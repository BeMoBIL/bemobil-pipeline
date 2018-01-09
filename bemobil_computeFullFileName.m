% compute full file names
% -----------------------
function res = bemobil_computeFullFileName(filePaths, fileNames)
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end;
end