clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing of results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file runs a few pipelines to postprocess the output in each
% working directoy. Please be patient, this can take up to a couple of
% minutes depending on the amount of files, trials, etc..
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         08/2022

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardcoded paths - nothing to do here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Which model do you want to use?
Model2Use = 'lernergopal'; % 'lernergopal', 'rajagopal', RajagopalLaiUhlrich2023 ....

% The root folder populating all of your working directories
rootDirectory = 'E:\LocDat\GitHub\autoSIM\_commonFiles\dataExamples\OSS\';

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
conditions =  {'WalkA'}; % e.g.: {'Dynamic'} or {'mSEBT_AT','mSEBT_PM','mSEBT_PL'}

% Set trial type for plotting: walking or not "something else" (=notWalking)
% This is only necessray to have the plots nicely formatted.
trialType = 'walking'; % default = 'walking'; 'walking' or 'notWalking';

% You can use the prefix (if you used it during processing) to split up the data report in different "groups"
prefix = {'standard'};  % default = {'standard'} or {'standard-CE00', 'standard-CE0MW', 'standard-CE150'}

% Force to time normalize all data to 100% activity time OR use settings from JAM
% processing. Note that the latter migth return partly time-normalized data (for
% joint contact force and non-normalized ones for e.g. muscle activity.
timeNormFlag = true; % default = true; true or false;

% Disable this for debugging when developing the code. Enable for running
% big file batches so that errors are caught and Matlab won`t stop on errors.
catchErrors = true; % default = true; true or false

% Use Matlab's parallel toolbox? Note that in this case errors are catched
% but not tracked since the parallel workers do not allow for this.
useParallelToolBox = true; % default = true; true or false
M = 20; % Set number of workers; for a 32-core server 20 seems ok.

%% Run postprocessing pipelines in loops

% Start timer
tStart = tic;

% Loop through prefix
for pref = 1 : length(prefix)

    % Set current prefix
    cPrefix = prefix{pref};

    % Set maximum N of WDs.
    Nmax = length(workingDirectories);

    % Define the number of splits and compute batch size
    if Nmax < 30 % If Nmax is below e.g., 30 we do not need to create separated batches, since this only takes longer than using a single batch.
        split = 1; % ==1; If Nmax is below e.g., 30 we do not need to create separated batches.
        batchSize = ceil(Nmax / split);  % Use ceil to ensure no directories are missed.
    else
        split = 20; % default = 5-20; this highly depends on your use case. E.g.: we processed 1.500 directorires, had 128GB RAM and needed 5 splits to prevent running out of RAM due to memory leak.
        batchSize = ceil(Nmax / split);  % Use ceil to ensure no directories are missed.
    end

    % Loop over each batch - since there is often a memory leak, this might help
    % to prevent that Matlab runs out of RAM, since after each batch the
    % parpool gets closed and associated RAM released.
    for batchIdx = 1:split
        % Compute the start and end index for the current batch
        startIdx = (batchIdx - 1) * batchSize + 1;
        endIdx = min(batchIdx * batchSize, Nmax);  % Ensure not to exceed Nmax

        % Run the parallel loop for the current batch
        for i_parfor = startIdx:endIdx; warning("parfor not activated!"); %#> for development only
        %parfor (i_parfor = startIdx:endIdx, M)
            loops4PostProcessingVarModels(rootDirectory, workingDirectories, conditions, path2setupFiles, cPrefix, timeNormFlag, catchErrors, trialType, Model2Use, i_parfor)
        end

        % Shut down parallel pool to release allocated RAM by the parpool.
        poolobj = gcp('nocreate');
        delete(poolobj);

        % Inform user about finished batch.
        disp(strcat('>>>>>  Batch #', string(batchIdx), {' '}, 'out of #', string(split), {' '}, 'finished. <<<<<'));
    end

    % Clean up WD_parfooLoop files to single WD file.
    cleanUpParforWDFiles4PostPro(fullfile(rootDirectory, strcat('_', Model2Use, '-WD-tmpLogFiles')), length(workingDirectories),  char(fullfile(rootDirectory, strcat('workingDirectories-', Model2Use, '_', cPrefix, '.xlsx'))));

end
%% Final  message
% End timer
tEnd = seconds(toc(tStart));
tEnd.Format = 'hh:mm';
disp(''); % new line
disp('************************************************************************************************************');
disp(strcat('>>>>>  Total processing duration:',{' '}, string(tEnd),'(hh:mm). <<<<<'));
disp('>>>>> Post processing of all results finished. <<<<<');
