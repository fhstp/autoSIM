function renameC3D2enfDescription(workingDirectory)
% This files finds all *.c3d & *.enf files in a working directory and
% renames the c3d files so  match the discriptions in the *.enf files

% Find all *.c3d files
cd(workingDirectory);
tmp_c3d = strtrim(string(ls('*.c3d')));
tmp_enf = strtrim(string(ls('*.enf')));

% Make sure that only *.enf files are listed which have a corresponding *.c3d file
tmpNamesc3d = erase(tmp_c3d,'.c3d'); % remove file ending
tmpNamesenf = erase(tmp_enf,'.Trial.enf'); % remove file ending
commonFileNames = intersect(tmpNamesc3d, tmpNamesenf);

% Make file list
files.c3d = strcat(commonFileNames,'.c3d');
files.enf = strcat(commonFileNames,'.Trial.enf');

% Now rename all files
description = {}; % initialize
for i = 1 : length(files.enf)
    
    % Get condition
    fileID = fopen(files.enf{i});
    tline = fgetl(fileID);
    while ischar(tline)
        %disp(tline)
        tline = fgetl(fileID);
        if ischar(tline) && contains(tline, 'DESCRIPTION=')
            description = char(tline(13:end));
            break           
        end
    end
    fclose(fileID);

    % In case there is no description still rename so there is no confusion
    % with other files
    if isempty(description)
        description = strcat('NoEnfDescriptionFound_', num2str(i));
    end    
    
    % Rename
    % c3d
    if ~strcmp(files.c3d{i}(1:end-4),description) % check if files were already renamed
        movefile(files.c3d{i}, strcat(description, '.c3d'), 'f'); % don´t use other characters than '_', as this will resuts later in an error where these names are used as fieldnames in structs
        disp('>>>>> Renamed one *.c3d file.')
    else
        disp('>>>>> Renamed *.c3d file detected - no changes made.')
    end
    
    % enf
    if ~strcmp(files.enf{i}(1:end-10),description) % check if files were already renamed
        movefile(files.enf{i}, strcat(description, '.Trial.enf'), 'f');
        disp('>>>>> Renamed one *.c3d file.')
    else
        disp('>>>>> Renamed *.enf file detected - no changes made.')
    end
end

disp('***************** *.c3d & *.enf files renamed *****************')

end

