function nCMD = getNCmdWindows()
% This function finds all open cmd windows with the name 
%"OpenConsole.exe" or cmd.exe".

% Search for open CMD windows

% Try one version
[~, cmdOutput] = system('tasklist /fi "imagename eq OpenConsole.exe" /v /fo csv');
searchStr = 'OpenConsole.exe';

% If nothing was found try another version
if contains(cmdOutput, 'INFO: No tasks are running which match the specified criteria.')

    [~, cmdOutput] = system('tasklist /fi "imagename eq cmd.exe" /v /fo csv');
    searchStr = 'cmd.exe';
	
end


% Split the output into lines
cmdOutputLines = strsplit(cmdOutput, '\n');
cmdOutputLines = strrep(cmdOutputLines,'"','');
nCMD = 0;

% Process each line to check if it's a CMD window and extract the PID
for i = 2:length(cmdOutputLines)-1

    line = strsplit(cmdOutputLines{i}, ',');

    % Check if it's a CMD window
    if ~isempty(line) && strcmp(line{1}, searchStr)
        nCMD = nCMD + 1;
    end
end

end