function createErrorReportsPerTrial(workingDirectory, rootDirectory, path2setupFiles, condition, prefix, trialType, Model2Use)
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
% Last changed:         08/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start timer
tStart = tic;

% Prepare input
workingDirectory = fullfile(workingDirectory);

%% Find all relevant report files
cd(strcat(workingDirectory,'Simulation\'));
Folder = cd;
if strcmp(prefix, '') % if prefix is empty
    FilesIKMrkLocations =   dir(fullfile(Folder, '**', strcat('*',condition,'*','\**\inverse-kinematics\'), '*_ik_model_marker_locations.sto')); % these files don´t have the prefix in their names
    FilesIKError =          dir(fullfile(Folder, '**', strcat('*',condition,'*','\**\inverse-kinematics\'), '*_ik_marker_errors.sto')); % these files don´t have the prefix in their names
    FilesActivationError =  dir(fullfile(Folder, '**', strcat('*',condition,'*','\**'), '*_activation.sto'));
    FilesForce =            dir(fullfile(Folder, '**', strcat('*',condition,'*','\**'), '*_force.sto'));
    FilesID =               dir(fullfile(Folder, '**', strcat('*',condition,'*_inverse-dynamics.sto')));
    FilesKinem =            dir(fullfile(Folder, '**', strcat('*',condition,'*_IK_motion_file.mot')));
    FilesCF =               dir(fullfile(Folder, '**', strcat('*',condition,'*','\**'), '*_ReactionLoads.sto'));
    FilesGRFs =             dir(fullfile(Folder, '**', strcat('*',condition,'*','\**'), '*_ForceReporter_forces.sto'));
else
%     FilesIKMrkLocations =   dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**\inverse-kinematics\'), '*_ik_model_marker_locations.sto')); % these files don´t have the prefix in their names
%     FilesIKError =          dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**\inverse-kinematics\'), '*_ik_marker_errors.sto')); % these files don´t have the prefix in their names
%     FilesActivationError =  dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**'), '*_activation.sto'));
%     FilesForce =            dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**'), '*_force.sto'));
%     FilesID =               dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*_inverse-dynamics.sto')));
%     FilesKinem =            dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*_IK_motion_file.mot')));
%     FilesCF =               dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**'), '*_ReactionLoads.sto'));
%     FilesGRFs =             dir(fullfile(Folder, '**', strcat('*',prefix,'*', condition,'*','\**'), '*_ForceReporter_forces.sto'));
    FilesIKMrkLocations =   dir(fullfile(Folder, '**',strcat('*',prefix,'*', condition,'*','\**\inverse-kinematics\'), '*_ik_model_marker_locations.sto')); % these files don´t have the prefix in their names
    FilesIKError =          dir(fullfile(Folder, '**',strcat('*',prefix,'*', condition,'*','\**\inverse-kinematics\'), '*_ik_marker_errors.sto')); % these files don´t have the prefix in their names
    FilesActivationError =  dir(fullfile(Folder, '**\static-optimization\', strcat('*',prefix,'*',condition,'*_activation.sto')));
    FilesForce =            dir(fullfile(Folder, '**\static-optimization\', strcat('*',prefix,'*',condition,'*_force.sto')));
    FilesID =               dir(fullfile(Folder, '**\inverse-dynamics\', strcat('*',prefix,'*',condition,'*_inverse-dynamics.sto')));
    FilesKinem =            dir(fullfile(Folder, '**',strcat('*',prefix,'*', condition,'*','\**\inverse-kinematics\'), '*_IK_motion_file.mot'));
    FilesCF =               dir(fullfile(Folder, '**\analyze\', strcat('*',prefix,'*',condition,'*_ReactionLoads.sto')));
    FilesGRFs =             dir(fullfile(Folder, '**\analyze\', strcat('*',prefix,'*', condition,'*_ForceReporter_forces.sto')));
end

%% Create File List

FilesIKMrkLocationsList = {}; % initialize cell
for i = 1 : length(FilesIKMrkLocations)
    FilesIKMrkLocationsList{i,1} = [FilesIKMrkLocations(i).folder,'\', FilesIKMrkLocations(i).name];
end

FilesIKErrorList = {}; % initialize cell
for i = 1 : length(FilesIKError)
    FilesIKErrorList{i,1} = [FilesIKError(i).folder,'\', FilesIKError(i).name];
end

FilesActivationErrorList = {}; % initialize cell
for i = 1 : length(FilesActivationError)
    FilesActivationErrorList{i,1} = [FilesActivationError(i).folder,'\', FilesActivationError(i).name];
end

FilesForceList = {}; % initialize cell
for i = 1 : length(FilesForce)
    FilesForceList{i,1} = [FilesForce(i).folder,'\', FilesForce(i).name];
end

FilesIDList = {}; % initialize cell
for i = 1 : length(FilesID)
    FilesIDList{i,1} = [FilesID(i).folder,'\', FilesID(i).name];
end

FilesKinemList = {}; % initialize cell
for i = 1 : length(FilesKinem)
    FilesKinemList{i,1} = [FilesKinem(i).folder,'\', FilesKinem(i).name];
end

FilesCFList = {}; % initialize cell
for i = 1 : length(FilesCF)
    FilesCFList{i,1} = [FilesCF(i).folder,'\', FilesCF(i).name];
end

FilesGRFList = {}; % initialize cell
for i = 1 : length(FilesGRFs)
    FilesGRFList{i,1} = [FilesGRFs(i).folder,'\', FilesGRFs(i).name];
end

%% Check if all file lists have the same size

% Deal unequal file lists
tmp = strrep(FilesKinemList, strcat(workingDirectory,'Simulation\'), ''); % start with the results which has highest chances to have results. This will help that I do not miss data during dealing!

cnt = 1;
for k = 1 : length(tmp)
    idxSlash = strfind(tmp{k},'\');
    idxSlash = idxSlash(1) -1;
    currentTrialName = tmp{k}(1:idxSlash);

    if any(contains(FilesIKMrkLocationsList, currentTrialName)) && any(contains(FilesIKErrorList, currentTrialName)) && ...
            any(contains(FilesActivationErrorList, currentTrialName)) && any(contains(FilesForceList, currentTrialName)) && ...
            any(contains(FilesIDList, currentTrialName)) && any(contains(FilesKinemList, currentTrialName)) && ...
            any(contains(FilesCFList, currentTrialName)) && any(contains(FilesGRFList, currentTrialName))

        % Deale lists
        FilesIKMrkLocationsList_dealed{cnt,1} = FilesIKMrkLocationsList{find(contains(FilesIKMrkLocationsList, currentTrialName))};
        FilesIKErrorList_dealed{cnt,1} = FilesIKErrorList{find(contains(FilesIKErrorList, currentTrialName))};
        FilesActivationErrorList_dealed{cnt,1} = FilesActivationErrorList{find(contains(FilesActivationErrorList, currentTrialName))};
        FilesForceList_dealed{cnt,1} = FilesForceList{find(contains(FilesForceList, currentTrialName))};
        FilesIDList_dealed{cnt,1} = FilesIDList{find(contains(FilesIDList, currentTrialName))};
        FilesKinemList_dealed{cnt,1} = FilesKinemList{find(contains(FilesKinemList, currentTrialName))};
        FilesCFList_dealed{cnt,1} = FilesCFList{find(contains(FilesCFList, currentTrialName))};
        FilesGRFList_dealed{cnt,1} = FilesGRFList{find(contains(FilesGRFList, currentTrialName))};

        % Increase cnt
        cnt = cnt + 1;
    end
end

% Check to make sure
lengths2check = [length(FilesIKMrkLocationsList_dealed), length(FilesIKErrorList_dealed), length(FilesActivationErrorList_dealed), ...
                 length(FilesForceList_dealed), length(FilesIDList_dealed), length(FilesKinemList_dealed), length(FilesCFList_dealed), ...
                 length(FilesGRFList_dealed)];

if all(diff(lengths2check) == 0)
    disp('>>>>> File lists look ok.')
else
    warning('Number of files is uneuqal for results categories! Please check!')
    return
end


%% Now do the rest
for i = 1 : size(FilesKinemList_dealed,1)
     try
        %% Load saved workspace variable of current working directory
        disp('>>>>> Starting to plot current batch.')
        % Get root folder of current subject
        id_tmp = strfind(FilesKinemList_dealed{i},'\');
        % Load saved work space
        ws = load(strcat(FilesKinemList_dealed{i}(1:id_tmp(end-1)),'workspace')); % this file is created during the end of <osimjam_workflow>.

        %% Load InputData File
        InputFile = dir(strcat(workingDirectory,'Simulation\','*', condition,'*InputData*.*'));
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

        [errors, hfig1] = generate_error_plots(path2IKerrorReport,path2MuscleActivationReport,path2MuscleForceReport,path2IDReport, ws.statesFileName, ws.cTO, ws.TO);


        %% Plot results to check quality
        % Prepare input data
        path2KneeKinemReport = FilesKinemList_dealed{i};
        path2ContactFReport = FilesCFList_dealed{i};
        path2GRFReport = FilesGRFList_dealed{i};

        path2setupFilesReport = ws.path.setupFiles;
        path2normData = fullfile(path2setupFiles,'..\..\_commonFiles\normData\');

        [results, hfig2, hfig3, hfig4] = generate_results_plots(trialType, workingDirectory, path2KneeKinemReport, path2ContactFReport, path2MuscleActivationReport, path2setupFilesReport, path2IDReport, path2MuscleForceReport, path2GRFReport, path2normData, ws.statesFileName, ws.side, ws.cTO, ws.TO, ws.IC, ws.cIC, ws.BW, Model2Use);

        %% Add Meta data to results
        results.metadata.BW = ws.BW;
        results.metadata.IC = ws.IC;
        results.metadata.ICi = ws.ICi;
        results.metadata.TO = ws.TO;
        results.metadata.cTO = ws.cTO;
        results.metadata.side = ws.side;

        %% Save figures and results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Handle empty prefix
        if isempty(prefix); prefix2save = strcat(prefix,'-');
		else
			prefix2save = prefix;
		end

        % Create output folder
        outputFolder = char(strcat(workingDirectory,'Simulation\allErrorReportFigures\'));
        if ~logical(exist(outputFolder, 'dir'))
            mkdir(outputFolder)
        end

        % Create output folder at top level
        outputFolderGroupData = char(fullfile(strcat(rootDirectory, '_', Model2Use, '-', prefix2save, '-groupData\IKerrors\')));
        if ~logical(exist(outputFolderGroupData, 'dir'))
            mkdir(outputFolderGroupData)
        end

        % Make sure that no other data from same individual are already
        % stored in the folder, if so add custome label to "subjectName"
        ls = dir(fullfile(outputFolderGroupData, strcat(subjectName, '*', ws.statesFileName,'*')));
        cfiles = {ls.name};
        nFiles = length(cfiles);

        if nFiles == 0
            subjectName2save = strcat(subjectName, '_', ws.statesFileName);
        else
            suffix = {'_a', '_b', '_c', '_d', '_e', '_f', '_g'};
            subjectName2save = strcat(subjectName, suffix{nFiles}, '_', ws.statesFileName);
        end

        % combine error struct
        errors.IKerrorLocs = IKerrorLocs;

        % Save stuff
        if exist('hfig0','var'); saveas(hfig0,char(strcat(outputFolderGroupData, subjectName2save, '_', ws.statesFileName,'_IKModelvsExperimentalMrkLocationErrorReport.png'))); end% also copy to top lvl folder for convenience


        if exist('hfig0','var'); saveas(hfig0,char(strcat(outputFolder, ws.statesFileName,'_IKModelvsExperimentalMrkLocationErrorReport.png'))); end
        if exist('hfig1','var'); saveas(hfig1,char(strcat(outputFolder, ws.statesFileName,'_ErrorReport.png'))); end
        if exist('hfig2','var'); saveas(hfig2,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportKinematicsForce.png'))); end
        if exist('hfig3','var'); saveas(hfig3,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportMuscleActivation.png'))); end
        if exist('hfig4','var'); saveas(hfig4,char(strcat(outputFolder, ws.statesFileName,'_ResultsReportGRFs.png'))); end
        if exist('results','var'); save(char(strcat(ws.path.results, ws.statesFileName,'FinalResultsStruct.mat')),'results'); end
        if exist('errors','var'); save(char(strcat(ws.path.results, ws.statesFileName,'ErrorStruct.mat')),'errors'); end

        %% Clear some vars
        close all
        clear ws
    catch
        %% Clear some vars and close figs
        close all
        clear ws
    end
end

%% FINAL %%
% End timer
tEnd = toc(tStart);
disp('********************************************');
fprintf('>>>>> Processing duration: %d minutes and %0.f seconds.\n', floor(tEnd/60), rem(tEnd,60));
disp('>>>>> Error reports were created ...');
end