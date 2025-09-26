function createErrorReportsPerTrial(workingDirectory, rootDirectory, path2setupFiles, condition, prefix, trialType, processVtp)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create error reports for each single trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file creates for each folder in the JAM results folder an error
% report.
%
% NOTE:
%       # You will have to define the path of the folder manually below.
%         Before you run this skript you will need to run the
%         <postprocessVtpFiles.m> script to ahve all vtp files converted to
%         ascii version so matlab can read them.
%
%       # If you want to process files with a prefix you have to add the
%         prefix in the Seetings section below!
%
%
% Input:    - workingDirectory = e.g. 'C:\LokaleDaten\CloudStation\MyPapers\2020_GuidedGrowth\data\AK-data\NocontactMScaling\';
%           - prefix = empty or e.g. 'contactEnergy100'; If empty all files will be used
%           - condition = e.g. 'Dynamic';  mandatory (e.g. Dynamic, Walking, .... a substring in the filenames)
%
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
%
% Last changed:         01/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start timer
tStart = tic;

% Prepare input
workingDirectory = fullfile(workingDirectory);

%% Find all relevant error report files

% Find all IK kinematics marker locations files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, ~] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesIKMrkLocationsList = strcat(folderPaths,'\inverse-kinematics\','_ik_model_marker_locations.sto');
FilesIKMrkLocationsList = FilesIKMrkLocationsList(isfile(FilesIKMrkLocationsList)); % get only files that exist

% Find all IK kinematics error files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, ~] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesIKErrorList = strcat(folderPaths,'\inverse-kinematics\','_ik_marker_errors.sto');
FilesIKErrorList = FilesIKErrorList(isfile(FilesIKErrorList)); % get only files that exist

% Find all muscle activation files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, folderNames] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesActivationErrorList = strcat(folderPaths,'\comak\',folderNames,'__activation.sto');
FilesActivationErrorList = FilesActivationErrorList(isfile(FilesActivationErrorList)); % get only files that exist

% Find all muscle forces files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, folderNames] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesForceList = strcat(folderPaths,'\comak\',folderNames,'__force.sto');
FilesForceList = FilesForceList(isfile(FilesForceList)); % get only files that exist

% Find all inverse-dynamics files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, folderNames] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesIDList = strcat(folderPaths,'\inverse-dynamics\',folderNames,'_inverse-dynamics.sto');
FilesIDList = FilesIDList(isfile(FilesIDList)); % get only files that exist

% Find all comak kinematics files
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, folderNames] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesKinemList = strcat(folderPaths,'\comak\',folderNames,'__values.sto');
FilesKinemList = FilesKinemList(isfile(FilesKinemList)); % get only files that exist

% Find all force reporter files - I want to get everything except the muscle in here. The muscle forces in the ForceReport are wrong!
Folder = strcat(workingDirectory,'JAM\');
[folderPaths, folderNames] = getAllFoldersToAnalyze(Folder, prefix, condition);
FilesCFList = strcat(folderPaths,'\jam\',folderNames,'__ForceReporter_forces.sto');
FilesCFList = FilesCFList(isfile(FilesCFList)); % get only files that exist

% Check if File lists are empty
if isempty(FilesKinemList)
    warning(strcat('No file(s) found in <',workingDirectory , '> for condition <', condition,'> and prefix <', prefix,'>.'));
else
    
    %% Deal unequal file lists
    tmp = strrep(FilesKinemList, strcat(workingDirectory,'JAM\'), ''); % start with the results which has highest chances to have results. This will help that I do not miss data during dealing!

    cnt = 1;
    for k = 1 : length(tmp)
        idxSlash = strfind(tmp{k},'\');
        idxSlash = idxSlash(1) -1;
        currentTrialName = tmp{k}(1:idxSlash);

        if any(contains(FilesIKMrkLocationsList, currentTrialName)) && any(contains(FilesIKErrorList, currentTrialName)) && ...
                any(contains(FilesActivationErrorList, currentTrialName)) && any(contains(FilesForceList, currentTrialName)) && ...
                any(contains(FilesIDList, currentTrialName)) && any(contains(FilesKinemList, currentTrialName)) && any(contains(FilesCFList, currentTrialName))

            % Deale lists
            FilesIKMrkLocationsList_dealed{cnt,1} = FilesIKMrkLocationsList{find(contains(FilesIKMrkLocationsList, currentTrialName))};
            FilesIKErrorList_dealed{cnt,1} = FilesIKErrorList{find(contains(FilesIKErrorList, currentTrialName))};
            FilesActivationErrorList_dealed{cnt,1} = FilesActivationErrorList{find(contains(FilesActivationErrorList, currentTrialName))};
            FilesForceList_dealed{cnt,1} = FilesForceList{find(contains(FilesForceList, currentTrialName))};
            FilesIDList_dealed{cnt,1} = FilesIDList{find(contains(FilesIDList, currentTrialName))};
            FilesKinemList_dealed{cnt,1} = FilesKinemList{find(contains(FilesKinemList, currentTrialName))};
            FilesCFList_dealed{cnt,1} = FilesCFList{find(contains(FilesCFList, currentTrialName))};

            % Increase cnt
            cnt = cnt + 1;
        end
    end

    % Check to make sure
    lengths2check = [length(FilesIKMrkLocationsList_dealed), length(FilesIKErrorList_dealed), length(FilesActivationErrorList_dealed), length(FilesForceList_dealed), length(FilesIDList_dealed), length(FilesKinemList_dealed), length(FilesCFList_dealed)];
    if all(diff(lengths2check) == 0)
        disp('>>>>> File lists look ok.')
    else
        warning('Number of files is uneuqal for results categories! Please check!')
        return
    end


    %% Now do the rest
    for i = 1 : size(FilesCFList_dealed,1)
        try
            %% Load saved workspace variable of current working directory
            disp('>>>>> Starting to plot current batch.')
            % Get root folder of current subject
            id_tmp = strfind(FilesKinemList_dealed{i},'\');
            % Load saved work space
            ws = load(strcat(FilesKinemList_dealed{i}(1:id_tmp(end-1)),'workspace')); % this file is created during the end of <osimjam_workflow>.

            %% Load InputData File
            InputFile = dir(strcat(workingDirectory,'JAM\','*', condition,'*InputData*.*'));
            load(horzcat(InputFile.folder,'\', InputFile.name));

            % Get file names
            FileNames_InputData = fieldnames(InputData);
            subjectName = char(strrep(InputData.(FileNames_InputData{1}).subjectName,' ','_'));

            %% Plot IK marker location vs. experimental marker
            path2MrkLocations = FilesIKMrkLocationsList_dealed{i};
            [hfig0, IKerrorLocs] = reportIKmarkerErrors(path2MrkLocations, ws.path2trc, workingDirectory, ws.statesFileName);

            %% Plot errors to check quality
            % Prepare input data
            path2IKerrorReport = FilesIKErrorList_dealed{i};
            path2MuscleActivationReport = FilesActivationErrorList_dealed{i};
            path2MuscleForceReport = FilesForceList_dealed{i};
            path2IDReport = FilesIDList_dealed{i};

            [errors, hfig1] = generate_error_plots_comak(path2IKerrorReport,path2MuscleActivationReport,path2MuscleForceReport,path2IDReport, ws.statesFileName, ws.cTO, ws.TO);

            %% Plot results to check quality
            % Prepare input data
            path2KneeKinemReport = FilesKinemList_dealed{i};
            path2ContactFReport = FilesCFList_dealed{i};
            path2setupFilesReport = ws.path.setupFiles;
            path2normData = fullfile(path2setupFiles,'..\..\_commonFiles\normData\');
            path2origGeometry = fullfile(path2setupFiles,'Models\_Geometry\');

            [results, hfig2, hfig3, hfig4, hfig5] = generate_results_plots_comak(trialType, workingDirectory, path2KneeKinemReport, path2ContactFReport, path2MuscleActivationReport, path2setupFilesReport, path2IDReport, path2MuscleForceReport, path2normData, path2origGeometry, ws.statesFileName, ws.side, ws.cTO, ws.TO, ws.IC, ws.cIC, ws.BW, processVtp);

            %% Add Meta data to results
            results.metadata.BW = ws.BW;
            results.metadata.IC = ws.IC;
            results.metadata.ICi = ws.ICi;
            results.metadata.TO = ws.TO;
            results.metadata.cTO = ws.cTO;
            results.metadata.side = ws.side;

            %% Save figures and results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Create output folder
            outputFolder = char(strcat(workingDirectory,'JAM\allErrorReportFigures\'));
            if ~logical(exist(outputFolder, 'dir'))
                mkdir(outputFolder)
            end

            % Create output folder at top level
            outputFolderGroupData = char(fullfile(strcat(rootDirectory,'_comak-groupData\IKerrors\')));
            if ~logical(exist(outputFolderGroupData, 'dir'))
                mkdir(outputFolderGroupData)
            end

            % Make sure that no other data from same individual are already
            % stored in the folder, if so add custome label to "subjectName"
            ls = dir(fullfile(outputFolderGroupData, strcat(subjectName,'*', prefix, condition)));
            cfiles = {ls.name};
            nFiles = length(cfiles);

            if nFiles == 0
                subjectName2save = subjectName;
            else
                suffix = 'a':'z';
                subjectName2save = strcat(subjectName, suffix(nFiles));
            end

            % Combine error struct
            errors.IKerrorLocs = IKerrorLocs;

            % Save stuff
            saveas(hfig0,char(strcat(outputFolderGroupData, subjectName2save, '_', ws.statesFileName,'_IKModelvsExperimentalMrkLocationErrorReport.png'))); % also copy to top lvl folder for convenience

            saveas(hfig0,char(strcat(outputFolder, ws.statesFileName,'_IKModelvsExperimentalMrkLocationErrorReport.png')));
            saveas(hfig1,char(strcat(outputFolder, ws.statesFileName,'_ErrorReport.png')));
            saveas(hfig2,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportKinematicsForce.png')));
            saveas(hfig3,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportMuscleActivation.png')));
            if processVtp
                saveas(hfig4,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportContactPressureMap_TF.png')));
                saveas(hfig5,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportContactPressureMap_PF.png')));
            end
            save(char(strcat(ws.path.COMAKresults, ws.statesFileName,'FinalResultsStruct.mat')),'results');
            save(char(strcat(ws.path.COMAKresults, ws.statesFileName,'ErrorStruct.mat')),'errors');
            close all
            clear ws
        catch
            %% Clear some vars and close figs
            close all
            clear ws
        end

        %% FINAL %%
        % End timer
        tEnd = toc(tStart);
        disp('********************************************');
        fprintf('>>>>> Processing duration: %d minutes and %0.f seconds.\n', floor(tEnd/60), rem(tEnd,60));
        disp('>>>>> Error reports were created ...');
    end
end