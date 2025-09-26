function CpuLoadBasedPausing(thresholdCpuLoad, waitInterval)
% This function estimates the CPU Load in a while loop. If
% <thresholdCpuLoad> is reached it will break out of the loop. This
% function can be used, to let Matlab wait until the CPU load is low, to
% run a script or, in case of COMAK, run the next batch of files.

% Input: 
% - thresholdCpuLoad: CPU load threhsold to reach to break out of the loop
% - waitInterval: the interval between wthe CPU load is estimated (in
% seconds!);

% Writte by: Brian Horsak
% Last changed. 03/2023
% -------------------------------------------------------------------------

pause('on');

% Initialize CPU load monitoring
cpuIdleProcess = 'Idle';
mon.NumOfCPU = double(System.Environment.ProcessorCount);
mon.ProcPerfCounter.cpuIdle = System.Diagnostics.PerformanceCounter('Process', '% Processor Time', cpuIdleProcess);
% mon.ProcPerfCounter.Matlab  = System.Diagnostics.PerformanceCounter('Process', '% Processor Time', MatlabProcess.ProcessName);
% MatlabProcess = System.Diagnostics.Process.GetCurrentProcess(); %// "Matlab" process

while true
    % Start timer
    tic
    % Calculate the cpu usage over one minute
    cpu_total =  [];

    % Start to get the CPU load
    mon.ProcPerfCounter.cpuIdle.NextValue;
    for i = 1 : 6
        cpu_total(i) = 100 - mon.ProcPerfCounter.cpuIdle.NextValue / mon.NumOfCPU;
        %cpu.matlab = mon.ProcPerfCounter.Matlab.NextValue / mon.NumOfCPU;
        %disp(cpu_total(i))
        pause(10)
    end
    mdnCPULoad = (median(cpu_total));
    disp(strcat('>>>>> Average CPU-Load:',{' '}, string(round(mdnCPULoad)),'% ... waiting to fall below', {' '}, string(thresholdCpuLoad),'% to continue!'));

    % Stop pausing if CPU load is below threhsold
    if mdnCPULoad < thresholdCpuLoad
        disp('>>>>> CPU Load threshold reached! Continuing ... <<<<');
        break % break out of while loop
    end

    % Stop e.g. 10 Minutes (= 60*10)
    pause(waitInterval-toc)
end
end