function renameFolders(workingDirectory, strs)
% This function checks all paths in the working directory and renames them
% if they contain strings specified in strs (e.g. ',', ' ') or German 
% umlauts (ä, ö, ü).

% IN:
%       workingDirectory = 'D:\TestFolder\'
%       strs = {',', ' ', '#'}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the folder contents
tmp = dir(fullfile(workingDirectory,'**'));

% Remove all files (isdir property is 0)
folders = tmp([tmp(:).isdir]);

% Remove '.' and '..'
folders = folders(~ismember({folders(:).name},{'.','..'}));

% Find all folders containing the strs or umlauts
umlauts = {'ä', 'ö', 'ü'};
allStrs = [strs, umlauts];
idx = contains({folders(:).name}, allStrs);
relevantFolders = {folders(idx).name};
relevantPaths = {folders(idx).folder};

% Sort from longest to shortest path so that renaming starts with
% the deepest folders. Otherwise, paths will be incorrect in the second loop.
[~, sortIdx]= sort(cellfun('length', cellstr(relevantPaths)),'descend');
relevantFolders = relevantFolders(sortIdx);
relevantPaths = relevantPaths(sortIdx);

% Define paths with wrong strings
paths = fullfile(relevantPaths, relevantFolders);

% Define replacement for umlauts
umlautReplacements = {'ae', 'oe', 'ue'};

for i = 1 : numel(paths)
    % Replace unwanted characters
    newPath = paths{i};
    for j = 1 : length(strs)
        newPath = strrep(newPath, strs{j}, ''); % delete unwanted strings
    end
    
    % Replace umlauts
    for j = 1 : length(umlauts)
        newPath = strrep(newPath, umlauts{j}, umlautReplacements{j});
    end
    
    % Rename path
    path = paths{i};
    movefile(path, newPath);
    
    % Write message
    disp(char(strcat('>>>>> Path changed to:', {' '}, newPath)));
end
end