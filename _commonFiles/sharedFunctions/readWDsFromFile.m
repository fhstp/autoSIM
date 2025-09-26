function [workingDirectories, staticC3dFiles] = readWDsFromFile(rootDirectory, startIdx, stopIdx)

% This file looks for a 'workindDirectories.xlsx' file located in the
% rootDirectory and reads it. With the specified index the user can set an
% index for which data to use.
% -------------------------------------------------------------------------

% Read the *.xlsx file
file2read = fullfile(rootDirectory, 'workingDirectories.xlsx');
data = readtable(file2read);

% Get info.
if isnan(stopIdx)
    workingDirectories = data.WorkingDirectories(startIdx:end)';
    staticC3dFiles = data.Statics(startIdx:end)';
else
    workingDirectories = data.WorkingDirectories(startIdx:stopIdx)';
    staticC3dFiles = data.Statics(startIdx:stopIdx)';
end