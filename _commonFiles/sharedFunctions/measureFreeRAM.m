% Function to measure and display free RAM percentage
function freeRAMPercentage =  measureFreeRAM
    % Get memory usage information
    [~, sys] = memory;
    
    % Calculate the percentage of free RAM
    freeRAMPercentage = (sys.PhysicalMemory.Available / sys.PhysicalMemory.Total) * 100;
    
    % Display the percentage of free RAM
    %fprintf('Free RAM: %.2f%%\n', freeRAMPercentage);
end

% Timer setup to run the measureFreeRAM function every 60 seconds
%t = timer('ExecutionMode', 'fixedRate', 'Period', 60, 'TimerFcn', @measureFreeRAM);

% Start the timer
%start(t);

% To stop the timer, you can use the following command:
% stop(t);

% Optionally, if you want to clean up the timer after stopping it, use:
% delete(t);