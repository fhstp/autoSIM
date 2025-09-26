function [flagExists, skipSideInBatch] = checkIfSecConstrFileExists(workingDirectory, file2find, side, checkInterval)
% This function checks every <checkInterval> time period if a <file2find> exists. Note you can
% use wildcards such as '*' in the file2find string, e.g. 'IK_second_coord_constraint_functions_right*.xml'.
% It also generates a flag for the side for which an sec. constrain
% simulation failed so that I can skipp upcoming trials of that side for
% the current subject.

pause('on');

%Start timer
ticStart = tic;
while true
    
    % Stop if file exists
    file = dir(fullfile(workingDirectory,'\JAM\', file2find));
    if ~isempty(file)
        disp('>>>>> Sec. constrain *.xml file found. Proceeding with calculations ... <<<<');
        flagExists = true;
        skipSideInBatch = 'none';
        break % break out of while loop
    end
 
    % Check how long it was since starting the while loop and stop if it is
    % longer than one hour
    tEnd = (toc(ticStart)) / 60; % in minutes
    if tEnd > 60
        flagExists = false;
        skipSideInBatch = side;
        break
    end

    % Stop e.g. 10 Minutes (= 60*10)
    pause(checkInterval)
    disp('>>>>> Checking if sec. constrain functions *.xml file exists ... <<<<');

end
