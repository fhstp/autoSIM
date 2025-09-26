function renameSEBTc3D(workingDirectory)
% This files finds all mSEBT *.c3d & *.enf files in a working directory and
% renames the files so that they have the directions AT - PM - PL in their
% filenames.

% Find all *.c3d files
cd(workingDirectory);
tmp_c3d = strtrim(string(ls('*.c3d')));

% Make sure that only *.enf files are listed which have a corresponding *.c3d file
tmp_enf = strings();
for k = 1:length(tmp_c3d)
    %tmp_enf(k,1) = strtrim(string(ls(strcat('*',tmp_c3d{k}(1:7),'*.enf'))));
    tmp_enf(k,:) = strcat(tmp_c3d{k}(1:end-4),'.Trial.enf');
end

% Exclude static, etc. ...
files.c3d = tmp_c3d(contains(tmp_c3d, 'mSEBT'));
files.enf = tmp_enf(contains(tmp_enf, 'mSEBT'));

% Now rename
for i = 1 : length(files.c3d)
    % Get condition
    fileID = fopen(files.enf{i});
    dat = fscanf(fileID,'%s');
    idx = strfind(dat,'DESCRIPTION='); % find position of description
    condition = dat(idx + 12 : idx + 13); % get text after 'DESCRIPTION='
    fclose(fileID);
    
    % Rename
    % c3d
    if ~contains(files.c3d{i},condition) % check if files were already renamed
        movefile(files.c3d{i}, strcat(files.c3d{i}(1:5), '_', condition,'_',files.c3d{i}(end-5:end)), 'f'); % don´t use other characters than '_', as this will resuts later in an error where these names are used as fieldnames in structs
    else
        disp('>>>>> Renamed *.c3d file detected - no changes made.')
    end
    % enf
    if ~contains(files.enf{i},condition) % check if files were already renamed
        movefile(files.enf{i}, strcat(files.enf{i}(1:5), '_', condition,'_',files.enf{i}(end-11:end)), 'f');
    else
        disp('>>>>> Renamed *.enf file detected - no changes made.')
    end
end

disp('***************** *.c3d & *.enf files renamed *****************')

%% Clear variables except output to prevet memory leak.
clearvars
end

