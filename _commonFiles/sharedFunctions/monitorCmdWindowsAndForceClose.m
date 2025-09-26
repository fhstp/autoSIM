% MONITOR_CMD_WINDOWS_AND_FORCE_CLOSE Monitors and manages cmd windows.
%   This function monitors the number of open cmd windows over a specified
%   duration. It checks at regular intervals if the number of cmd windows
%   has decreased. If the number of cmd windows
%   does not decrease within the entire duration, it forces
%   all cmd windows to close. If the number of cmd windows decreases,
%   it restarts the monitoring.
%
%   Usage:
%       monitorCmdWindowsAndForceClose(num_checks, wait_interval, n_cmd, task)
%
%   Parameters:
%       num_checks - Number of checks to perform during the monitoring period.
%       wait_interval - Time interval (in seconds) between each check.
%
%   Example:
%       monitorCmdWindowsAndForceClose(5, 60); 5 x 60 sec. == 5 minutes
%
%   See also:
%       getNCmdWindows, closeCmdWindows
%
%   Author: Brian Horsak
%   Date: 03/2025

function monitorCmdWindowsAndForceClose(num_checks, wait_interval)

while true
    % Initialize array to store the number of open cmd windows
    cmd_counts = zeros(num_checks, 1);

    % Monitor the number of open cmd windows
    for i = 1:num_checks
        cmd_counts(i) = getNCmdWindows;
        pause(wait_interval);
    end

    total_time = wait_interval * num_checks;
    total_N_cmds = getNCmdWindows();

    % Final check if the number of open cmd windows decreased.
    if total_N_cmds == 0 || all(cmd_counts(1) == cmd_counts(2:end)) % make sure there was no change for all iterations.
        fprintf('Number of open cmd-windows did not decrease within the last %d seconds. I will force all cmd windows to close now ...\n', total_time);
        closeCmdWindows();
        break;
    else
        fprintf('Number of open cmd-windows is still decreasing. I will wait...\n');
    end

end
end