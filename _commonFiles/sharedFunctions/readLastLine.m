function lastLine = readLastLine(filename)
% This function reads the last line of a text file.

    % Open the file for reading
    fid = fopen(filename, 'rt');
    
    % Check if the file opened successfully
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Initialize variables
    lastLine = '';
    currentLine = fgetl(fid);
    
    % Loop through the file until the end
    while ischar(currentLine)
        lastLine = currentLine;
        currentLine = fgetl(fid);
    end
    
    % Close the file
    fclose(fid);

%% Clear variables except output to prevet memory leak.
clearvars -except lastLine
end
