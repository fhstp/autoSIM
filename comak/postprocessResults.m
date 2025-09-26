clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing of results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs a few pipelines to postprocess the output in each
% working directoy. Please be patient, this can take up to a couple of
% minutes depending on the amount of files, trials, etc..
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         04/2022

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardcoded paths - nothing to do here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For debugging only
% dbstop if error  % turn on
% dbclear if error % turn off

% Clean Matlab Path (optional)
restoredefaultpath

% Get path where this file is located on disc and set that as the repoPath
% Note that the file needs to to located in its original location, 
% within the repo top lvl folder.
tmpPath = mfilename('fullpath');
[repoPath, ~, ~] = fileparts(tmpPath);
repoPath = strcat(repoPath,'\');

path2setupFiles = fullfile(repoPath, '\setupFiles\');
addpath(genpath(fullfile(repoPath)));

% Also add common Files
addpath(genpath(fullfile(strcat(repoPath,'\..\_commonFiles'))));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Settings - make your choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ===== Set the data to process ==========================================

% The root folder populating all of your working directories
rootDirectory = 'E:\LocDat\GitHub\autoSIM\_commonFiles\dataExamples\OSS'; 

% Select Option 1 out of 3:
Option = 1;

switch Option
    case 1
        %----- Option #1 ----------------------------------------------------------
        % Set manual list of workingDirectories and staticC3dFiles:
        % You can only specify one single folder. These folders contain all *.c3d files to analyze for a single subject
        % OpenSim automatically looks for that folder here. If not found in workingDirectory\Geometry, it will look at the standard OpenSim Paths.
        % Make sure that all paths have a '\' at the end!
        %---
        workingDirectories =   {'E:\LocDat\GitHub\autoSIM\_commonFiles\dataExamples\OSS\'};    % {'D:\...\', 'C:\...\', ...} or {'D:\....\'}

    case 2
        %----- Option #2 ----------------------------------------------------------
        % Create workingDirectories and staticC3dFiles automatically and store a local file in the rootDirectory for later usage.
        % If the following line is uncommented it will create the lists of <workingDirectories> and <staticC3dFiles>
        % automatically based on the specified string pattern, e.g., 'Static' or 'Stand', etc.
        % You can choose between two methods: 'byC3dFilePatternName' == looking for the appropriate file by the c3d file name pattern provided by the last
        % inputvar (e.g., 'static') OR by 'byEnfDescription' which will look for the file based on the description in the *.enf file. The latter is the standard for the OSS.
        % Note that the third input is a string pattern for which the script would look at the end of potential static trials, i.e. "A" for "StandA" trials.
        %---
        [workingDirectories, ~] = getTopLvlFoldersStatics(rootDirectory, 'byEnfDescription', 'Stand', ''); % default (OSS only!) = 'byEnfDescription', 'Stand', '';  default (general) = 'byC3dFilePatternName', 'Static', '';

    case 3
        %----- Option #3 ----------------------------------------------------------
        % Read [workingDirectories, staticC3dFiles] from local file which was
        % created either manually or by Option #2 in an earlier stage. You can
        % specify here from which index it should start. This is the best option for large-scale simulation tasks.
        %---
        startIdx = 1; % default = 1
        stopIdx = NaN; % default = NaN (== to last wd); you can also specify a stop index here.
        [workingDirectories, ~] = readWDsFromFile(rootDirectory, startIdx, stopIdx);
end 

%% ===== More settings ====================================================

% Define several conditions to loop through. You also can have only one condition
conditions =  {'WalkA'}; % e.g.: {'Walk'} or {'mSEBT_AT','mSEBT_PM','mSEBT_PL'}

% Set trial type for plotting: walking or not "something else" (=notWalking)
% This is only necessray to have the plots nicely formatted.
trialType = 'walking'; % default = 'walking'; 'walking' or 'notWalking';

% You can use the prefix (if you used it during processing) to split up the data report in different "groups"
prefix = {'willi-v3'};  % default = {'standard'} or {'standard-CE00', 'standard-CE0MW', 'standard-CE150'}

% Force to time normalize all data to 100% activity time OR use settings from JAM
% processing. Note that the latter migth return partly time-normalized data (for
% joint contact force and non-normalized ones for e.g. muscle activity.
timeNormFlag = true; % default = true; true or false;

% Disable this for debugging when developing the code. Enable for running
% big file batches so that errors are caught and Matlab won`t stop on errors.
catchErrors = false; % default = true; true or false

% Shall the reporter process the *.vtp files. This is necessary if you want
% to get pressure maps for the TF/PF JCF separated for the medial and
% lateral components.
processVtp = false; % default = true; true or false

% Use Matlab's parallel toolbox? Note that in this case errors are catched
% but not tracked since the parallel workers do not allow for this.
useParallelToolBox = true; % default = true; true or false
M = 20; % Set number of workers; for a 32-core server 20 seems ok.

%% Run postprocessing pipelines in loops
loops4PostProcessing(rootDirectory, workingDirectories, conditions, path2setupFiles, prefix, timeNormFlag, catchErrors, trialType, 'comak', processVtp, useParallelToolBox, M)