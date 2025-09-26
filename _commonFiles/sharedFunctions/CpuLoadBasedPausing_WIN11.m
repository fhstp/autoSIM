function CpuLoadBasedPausing_WIN11(thresholdCpuLoad, waitInterval)
% This function estimates the CPU Load in a while loop. If
% <thresholdCpuLoad> is reached, it will break out of the loop.
% This function can be used to let Matlab wait until the CPU load is low
% before running a script or processing the next batch of files.
%
% Input:
% - thresholdCpuLoad: CPU load threshold to reach to break out of the loop
% - waitInterval: the interval between CPU load estimations (in seconds)
%
% Written by: Brian Horsak
% Modified: 03/2025
% -------------------------------------------------------------------------

pause('on');

while true
    % Start timer
    tic

    % Execute typeperf command and capture output
    [status, cmdout] = system('wmic cpu get loadpercentage'); %system('typeperf "\Processor Information(_Total)\% Processor Utility" -sc 6');

    if status ~= 0
        warning('Failed to retrieve CPU load. Make sure typeperf is available for your PC. I will wait for 5 minutes and then proceed.');
        pause(60*5) % wait five minutes.
        break
    end

    % Get CPU load for 60 seconds.
    cpuValues = [];
    for i = 1:6
        % Extract CPU load values from command output
        %cpuValues = regexp(cmdout, '([0-9]+\.[0-9]+)', 'match');
        %cpuValues = str2double(cpuValues);
        pattern = 'LoadPercentage\s+(\d+)';
        matches = regexp(cmdout, pattern, 'tokens');
        cpuValues(i) = str2double(matches{1}{1});
        pause(10)
    end

    % Clean potential nans.
    cpuValues = cpuValues(~isnan(cpuValues));

    if isempty(cpuValues)
        warning('No CPU values retrieved. Retrying...');
        pause(waitInterval - toc);
        continue;
    end

    % Compute median CPU load
    mdnCPULoad = median(cpuValues);

    % Stop pausing if CPU load is below threshold
    if mdnCPULoad < thresholdCpuLoad
        disp('>>>>> CPU Load threshold reached! Continuing ... <<<<');
        break;
    end

    % Display status
    disp(strcat('>>>>> Average CPU Load:', {' '}, string(round(mdnCPULoad)), '% ... waiting to fall below', {' '}, string(thresholdCpuLoad), '% to continue!'));

    % Pause for the remaining interval time
    pause(waitInterval - toc);
end

end
