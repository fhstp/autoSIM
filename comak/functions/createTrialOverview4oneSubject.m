function createTrialOverview4oneSubject(workingDirectory, rootDirectory, path2setupFiles, condition, prefix, timeNormFlag, trialType, processVtp)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create reports for one single subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file creates for all plots in one JAM folder (of a single subject) a
% data report. It only reports the data which belog to a specific
% condition. The prefix is optional, if left empty 'prefix' will not be
% used to filter the data used for plotting. If for example prefix =
% 'test', only files will be selected which hold the substring 'test' in
% their names. 
%
% NOTE: 
%          # You will have to define the path of the folder manually below.
%
%
% Inout:    # workingDirectory: e.g. 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\data\AK-data\ContactEnergy\' 
%           # rootDirectory: e.g.  top level dir containing alls working dirs         
%           # path2setupFiles:  Path to the repository stored somewhere on your hard disc (necessary for norm data display)
%             e.g. 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\code\my-comak-workflow\osimJAM\osimProcessing-GitLab\setupFiles\';
%           # prefix = Define type of files to analyze, leave empty '' or use a string, e.g. 'test'
%           # condition = 'Dynamic';   % Mandatory (e.g. Dynamic, Walking, .... a substring in the filenames)
%           # indicates if all data will be time-normalized to 100 %
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
%
% Last changed:         03/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start timer
tStart = tic;

% Prepare input
workingDirectory = fullfile(workingDirectory);
path2setupFiles = fullfile(path2setupFiles);

% Find and add all JAM output from the working directory for one condition to one output file and save plots.
subResults = analyzeSubjectFiles(workingDirectory, rootDirectory, path2setupFiles, condition, prefix, timeNormFlag, trialType, processVtp);

%% FINAL %%

% End timer
tEnd = toc(tStart);
disp('********************************************');
fprintf('Data Report processing duration: %d minutes and %0.f seconds.\n', floor(tEnd/60), rem(tEnd,60));
disp('Data report was created ...');
end