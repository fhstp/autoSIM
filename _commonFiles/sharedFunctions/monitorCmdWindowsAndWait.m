% MONITOR_CMD_WINDOWS_AND_WAIT Monitors cmd windows and waits for changes.
%   This function monitors the number of open cmd windows over a specified
%   duration. It checks at regular intervals if the number of cmd windows
%   has increased or decreased by a specified threshold. If the number of
%   cmd windows does not change as expected, it continues waiting.
%
%   Usage:
%       monitorCmdWindowsAndWait(num_checks, wait_interval, N_cmd)
%
%   Parameters:
%       num_checks - Number of checks to perform during the monitoring period.
%       wait_interval - Time interval (in seconds) between each check.
%       N_cmd - Threshold for the number of cmd windows.
%
%   Example:
%       monitorCmdWindowsAndWait(5, 60, 50)
%
%   See also:
%       getNCmdWindows, closeCmdWindows
%
%   Author: Brian Horsak
%   Date: 03/2025

function monitorCmdWindowsAndWait(num_checks, wait_interval, max_cmd)

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

    % Check if the number of open cmd windows is above threshold and increasing.
    if  total_N_cmds > max_cmd && cmd_counts(1) <= cmd_counts(end)
        fprintf('\nThe number of open cmd-windows (N=%d) is above the max. allowed threshold (N=%d) and was not decreasing within the last %d seconds. I will keep waiting ...\n\n', total_N_cmds, max_cmd, total_time);
    else
		if total_N_cmds < max_cmd
			fprintf('\nI will proceed ... The number of open cmd-windows (N=%d) is below the max. allowed threshold (N=%d).\n\n', total_N_cmds, max_cmd);
		else
			fprintf('\nI will proceed ... The number of open cmd-windows (N=%d) is still above the max. allowed threshold (N=%d) but cmd-windows are decreasing over the last %d seconds.\n\n', total_N_cmds, max_cmd, total_time);
		end
        
		% Jump out of while loop.
        break
    end
end
end