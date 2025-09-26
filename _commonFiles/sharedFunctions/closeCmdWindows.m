function closeCmdWindows()
% This function finds all open cmd windows with the name "OpenConsole.exe"
% and forces them to close.

% Search for open CMD windows
[~, cmdOutput] = system('tasklist /fi "imagename eq OpenConsole.exe" /v /fo csv');

% Split the output into lines
cmdOutputLines = strsplit(cmdOutput, '\n');
cmdOutputLines = strrep(cmdOutputLines,'"','');

% Process each line to check if it's a CMD window and extract the PID
for i = 2:length(cmdOutputLines)-1

    line = strsplit(cmdOutputLines{i}, ',');

    % Check if it's a CMD window
    if ~isempty(line) && strcmp(line{1}, 'OpenConsole.exe')
        % Extract the PID
        pid = str2double(line{2});

        % Close the CMD window using taskkill command
        system(['taskkill /f /pid ' num2str(pid)]);
    end
end
end