function [folderNames_out, staticNames_out] = getTopLvlFoldersStatics(rootDirectory, method, staticNamePattern, staticNamePattern2)
% First, this function looks in the specifid folder for *.c3d files and
% gets a unique list of folder paths populating *.c3d files. It
% then looks into these folders and creates a list of stand/static files
% (based on the pattern specified with 'staticName'. When more than one
% file is availabe, it assumes that the trial with the highest number in
% its filename is the appropriate one. You also can specify an alternative
% pattern using staticNamePattern2 in case the first one was not found.

% It will also save folderNames_out, staticNames_out to a spreadsheet and
% store it in the rootDirectory as "files2process".

% -------------------------------------------------------------------------

disp('****************************************  Getting Folders and statics ...  *****************************************');

% Get all folder paths that populate *.c3d files.
files = dir(fullfile(strcat(rootDirectory,'/**/*.c3d')));

%% Get the workingDirectories.

% Get a logical vector that tells which is a directory.
dirFlags = ~[files.isdir];

% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.

% Get only the unique folder names into a cell array.
tmpFolders = fullfile(unique({subFolders.folder})); % Start at 3 to skip . and ..

% Add trailing '\' to the paths.
for j = 1:length(tmpFolders)
    FolderNames{j} = strcat(tmpFolders{j}, '\');
end

%% Get the static files.

for i = 1:length(FolderNames)

    switch method
        case 'byC3dFilePatternName'
            % Get the files.
            tmpFiles = dir(fullfile(strcat(FolderNames{i},'\*.c3d')));

            % Get a logical vector that tells which is a directory.
            fileFlags = ~[tmpFiles.isdir];
            tmpFiles = tmpFiles(fileFlags);

            % Get stand or static names.
            tmpOut = tmpFiles(startsWith({tmpFiles.name}, staticNamePattern, 'IgnoreCase', true));
            tmpOut = cellfun(@(x) x(1:end-4), {tmpOut.name}, 'un', 0); % Remove *.c3d

            % If empty try different approach.
            if isempty(tmpOut)
                tmpOut = tmpFiles(contains({tmpFiles.name}, staticNamePattern, 'IgnoreCase', true));
                tmpOut = cellfun(@(x) x(1:end-4), {tmpOut.name}, 'un', 0); % Remove *.c3d
            end

            % If still empty raise error.
            if isempty(tmpOut)
                warning off backtrace
                warning(strcat('No static file found while creating the static input file lists for <', FolderNames{i}, '>.'));
                warning on backtrace
                staticNames_out{i} = 'NoFileFound';
                folderNames_out{i} = FolderNames{i};

            elseif length(tmpOut) == 1
                staticNames_out{i} = char(strcat(tmpOut,'.c3d')); % add file extension back
                folderNames_out{i} = FolderNames{i}; % only use folder name when static was found

            else
                % Assuming that the trial number ist last and that the trial with
                % the highest number is the relevant one.
                static_Idx = [];
                for j = 1:length(tmpOut)
                    static_Idx(j) = str2double(regexp(tmpOut{j},'\d.','match'));
                end
                [~, idx] = max(static_Idx); %get the latest file

                staticNames_out{i} = char(strcat(tmpOut{idx},'.c3d')); % add file extension back
                folderNames_out{i} = FolderNames{i}; % only use folder name when static was found
            end

        case 'byEnfDescription'

            % Get the files
            cd(FolderNames{i});

            % Find all *.enf files
            tmp_enf = strtrim(string(ls('*.*Trial*.enf')));

            % Check if tmp_enf is empty
            if strcmp(tmp_enf, "")

                % Nothing found.
                staticNames_out{i} = 'NoFileFound';
                folderNames_out{i} = FolderNames{i};
                warning off backtrace
                warning(strcat('Still no "static trial" found in <', FolderNames{i}, '>!" ... '));
                warning on backtrace

            else

                % Make file list
                files = struct();
                files.enf = tmp_enf;

                % Now find the file with the Pattern
                potFiles = {}; % Init ialize
                for k = 1 : length(files.enf)
                    % Get condition
                    fileID = fopen(files.enf{k});
                    tline = fgetl(fileID);
                    while ischar(tline)
                        %disp(tline)
                        tline = fgetl(fileID);
                        if ischar(tline) && contains(tline, 'DESCRIPTION=')
                            potFiles{k} = tline(13:end);
                            break
                        end
                    end
                    fclose(fileID);
                end

                % Get the file
                staticsOut = {};
                potFiles = potFiles(~cellfun('isempty',potFiles)); % remove empty cells
                staticsOut = potFiles(find(contains(potFiles, staticNamePattern,"IgnoreCase",true)));

                % If no files were found try different patterns: "stand" and "static"
                if isempty(staticsOut)
                    warning off backtrace
                    warning(strcat('No "static trial" found with pattern <', staticNamePattern, '> in <', FolderNames{i}, '>!" ... I am trying Stand and Static ...'));
                    warning on backtrace
                    staticsOut = potFiles(find(contains(potFiles, {'Static', 'Stand'},"IgnoreCase",true)));

                    if isempty(staticsOut)
                        % Do not raise error here, error handling is introduced later
                        warning off backtrace
                        warning(strcat('Still no "static trial" found in <', FolderNames{i}, '>!" ... '));
                        warning on backtrace
                    end
                end

                % Now make selection
                if length(staticsOut) == 1 % if one was found
                    staticNames_out{i} = char(strcat(staticsOut,'.c3d')); % add file extension back
                    folderNames_out{i} = FolderNames{i}; % only use folder name when static was found
                    disp(['<<<<< Found static trial ... ', char(staticNames_out{i})]);

                elseif length(staticsOut) > 1 % if more were found

                    if ~isempty(staticNamePattern2)
                        % Try to find static based on the condition Pattern Walk"A"
                        idx = find(endsWith(staticsOut, staticNamePattern2, 'IgnoreCase', true));

                        if ~isempty(idx) % if found use it
                            staticNames_out{i} = char(strcat(staticsOut(idx),'.c3d'));
                            folderNames_out{i} = FolderNames{i}; % only use folder name when static was found
                            disp(['<<<<< Found static trial ... ', char(staticNames_out{i})]);
                        else
                            % If walk "A" condition is not reflected by static trial names

                            % Sort the files and get the first one
                            [~,idx] = sort(lower(staticsOut)); % make sure that the case does not affect sorting
                            staticsOut = staticsOut(idx);

                            staticNames_out{i} = char(strcat(staticsOut(end),'.c3d')); %take first one
                            folderNames_out{i} = FolderNames{i}; % only use folder name when static was found
                            % Inform user
                            warning off backtrace
                            warning(strcat('"Static" trial name convention does not follow the "WalkA" convention in <', FolderNames{i}, '>. I will use the last one: ... ',char(staticNames_out{i}),'!'));
                            warning on backtrace
                            disp(['<<<<< Selected static trial ... ', char(staticNames_out{i})]);
                        end

                    else
                        % Inform user if more than one file was found.

                        % Sort the files and get the first one.
                        [~,idx] = sort(lower(staticsOut)); % make sure that the case does not affect sorting
                        staticsOut = staticsOut(idx);

                        staticNames_out{i} = char(strcat(staticsOut(1),'.c3d')); %take first one
                        folderNames_out{i} = FolderNames{i}; % only use folder name when static was found
                        warning off backtrace
                        warning(strcat('More than one "static" trial found in <', FolderNames{i}, '>. I will use the first one: ... ',char(staticNames_out{i}),'!'));
                        warning on backtrace
                    end
                else
                    % Still nothing found
                    % Do not raise error here, error handling is introduced later
                    staticNames_out{i} = 'NoFileFound';
                    folderNames_out{i} = FolderNames{i};
                    warning off backtrace
                    warning(strcat('Still no "static trial" found in <', FolderNames{i}, '>!" ... '));
                    warning on backtrace
                end
            end
    end
end

% Now store folderNames_out, staticNames_out to a spreadsheet and
% store it in the rootDirectory as a files.
wdTable = cell2table([num2cell(1:length(staticNames_out))', folderNames_out', staticNames_out', cell(length(staticNames_out),1), cell(length(staticNames_out),1)], 'VariableNames',{'Index', 'WorkingDirectories', 'Statics', 'Simulation', 'Processed'});
destination = fullfile(rootDirectory, 'workingDirectories.xlsx');
writetable(wdTable, destination);

%% Clear variables except output to prevet memory leak.
clearvars -except folderNames_out staticNames_out
end