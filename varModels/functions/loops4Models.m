function loops4Models(rootDirectory, workingDirectories, staticC3dFiles, conditions, labFlag, path2GenericModels, path2bin, path2opensim, path2setupFiles, tf_angle_r, tf_angle_l, ...
    firstContact_L, firstContact_R, prefix, timeNorm, maxCmd, thresholdCpuLoad, catchErrors, lockSubtalar4Scaling, ...
    scaleMuscleStrength, manualMusScaleF, markerSet, bodyheightGenericModel, addPelvisHelperMarker, pelvisMarker4nonUniformScaling, tf_angle_fromSource, torsiontool, useStatic4TibRotEstimationAsFallback, ...
    tib_torsion_LeftMarkers, tib_torsion_RightMarkers, forceTrcMotCreation, dataAugmentation, ForceModelCreation, performPostProcessing, trialType, timeNormFlag, renameC3DFiles2enfDescription, ...
    Model2Use, tasks, thresholdFreeRAM, allowAutoRestart, checkAndAdaptMomArms, useASTool, repoPaths, useStatic4FrontAlignmentAsFallback, useC3Devents, useCPUThreshold)

% In case user selected the option that matlab can restart automatically in
% case of low memory.

% The follwoing <currentPath> needs to be created automatically since for a restart we do not have
% any inputs to the function before loading the temp workspace.
tmpPath = mfilename('fullpath');
[currentPath, ~, ~] = fileparts(fileparts(fileparts(tmpPath))); % jump up three folder lvls.
currentPath = fullfile(currentPath);

% Check if a workspace file exists to determine if we are resuming.
checkPointFilePath = fullfile(currentPath, 'workspaceCheckpoint4AutoRestart.mat');
if exist(checkPointFilePath, 'file')
    load(checkPointFilePath);
    delete(checkPointFilePath); % Delete the file after loading to prevent further restarts

    % Set the paths for Matlab again since they got lost by the restart.
    addpath(genpath(repoPaths.repo));
    addpath(genpath(repoPaths.commonFiles));

    % Make user aware that workspace was loaded from file instead of
    % created to prevent an unintended loading.
    warning("The workspace was loaded from an earlier *.tmp file. If this was not intended check if an old <workspaceCheckpoint4AutoRestart.mat is located in the repo folder.>")
    pause(60*5)
else
    i_wd = 1; % Initialize i_wd if starting fresh
    kickOffRestart = false; % Set to false initially
end


% Now loop through folders
errN = 0; % Number of errors occured and catched
failures = struct('cond', {}, 'wd', {}, 'err', {}); % collect catched errors
tStart = tic; % Start timer
numFiles = 1; %count number of files processed
overFlowThreshold = 1000; % default = 1000 (500-1000); Used to pause the skript every eg 500 trials and then force to close all open cmd windows.

% Batchcount for comak for each trial in the working directory
batchCount = 0; % Used to limit the number of open cmd windows

% Loop through wd
for i_wd = i_wd : length(workingDirectories)
    try
        % Set current workingDirectory & static trial
        % Check if comak subfolder in working directory exists, otherwise
        % create folder
        rootWorkingDirectory = fullfile(workingDirectories{i_wd});

        % Set current prefix
        cPrefix = prefix{1}; % I either get one prefix from the current workflow OR a set of prefixes from the dataaugmentation mode

        % Create prefix-specific WD
        workingDirectory = fullfile(strcat(workingDirectories{i_wd}, Model2Use, '-', cPrefix,'\'));

        if ~logical(exist(workingDirectory, 'dir'))
            mkdir(workingDirectory)
        end

        % Set static file
        staticC3d = staticC3dFiles{i_wd};

        % Check if it is a SEBT trial and then rename files
        if sum(contains(conditions, 'SEBT')) > 0
            % Prepare mSEBT data and rename *.c3d/*.enf files in the working directory so that each file contains the mSEBT direction
            renameSEBTc3D(workingDirectory);
        end

        %Rename c3d files in accordacne with the *.enf file
        %description, this is an OSS specific workflow and necessary to
        %get the correct naming for the conditions!
        if renameC3DFiles2enfDescription
            disp('>>>>> All *.c3d files will be renamed to the *.enf description!')
            renameC3D2enfDescription(rootWorkingDirectory);
        end

        % Loop for condition
        for cond = 1 : size(conditions,2)
            try
                % Initialize var to track the number of processed files in a WD for specific cond.
                proFilesCount = 0;

                % Set mSEBT condition: 'mSEBT-AT','mSEBT-PM','mSEBT-PL' and run comak separately
                condition = conditions{cond};

                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % RUN Pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Start message
                disp(['****************************************  ', upper(Model2Use), ' pipeline started  **********************************************']);

                % Calculate Input data for workflow for all files specified by 'condition' and creates all *.trc and *.mot files
                % Note that the last input defines how the grf data are rotated, e.g. for the 'FHSTP-BIZ' lab or the 'OSS' lab.
                [InputData] = prepareInputData(rootWorkingDirectory, workingDirectory, staticC3d, condition, labFlag, markerSet, addPelvisHelperMarker, pelvisMarker4nonUniformScaling, forceTrcMotCreation, Model2Use, useC3Devents);

                % Get trial names for all condition files
                trials = fieldnames(InputData);

                % Copy Geometry folder once to working directory
                if ~logical(exist(fullfile(workingDirectory,'Geometry'),'dir'))
                    copyfile(fullfile(path2GenericModels,'_Geometry'), fullfile(workingDirectory,'Geometry'));
                end

                % Prepare (scale) the models for both body sides and also scale contact geometry manually, because these are not scaled by opensim (seems to be a bug)
                bodymass = InputData.(trials{1}).Bodymass;      % get BW from first trial
                BodyHeight = InputData.(trials{1}).BodyHeight;      % get body mass from first trial
                static = InputData.(trials{1}).static_trcPath;  % static trc path from first trial
                sub_name = InputData.(trials{1}).subjectName;   % subject name

                % Create one sample *.trc file for R or L. This is used for the function <adaptWrappingObjects...> that checks the muscle moment
                % arms in the prepareScaledModel function.
                if checkAndAdaptMomArms

                    % Get one trial for L/R for IK.
                    idxTrialsR = find(contains(cellfun(@(x) x(end), trials, 'UniformOutput', false),'r'),1);
                    idxTrialsL = find(contains(cellfun(@(x) x(end), trials, 'UniformOutput', false),'l'),1);

                    % Check if a file was found for R.
                    if ~isempty(idxTrialsR)
                        sampleInputR = InputData.(trials{idxTrialsR});
                        checkAndAdaptMomArmsR = true;
                    else
                        sampleInputR = 'NA';
                        checkAndAdaptMomArmsR = false;
                    end

                    % Check if a file was found for L.
                    if ~isempty(idxTrialsL)
                        sampleInputL = InputData.(trials{idxTrialsL});
                        checkAndAdaptMomArmsL = true;
                    else
                        sampleInputL = 'NA';
                        checkAndAdaptMomArmsL = false;
                    end
                else
                    sampleInputR = 'NA';
                    sampleInputL = 'NA';
                    checkAndAdaptMomArmsR = false;
                    checkAndAdaptMomArmsL = false;
                end

                % Select L or R trial randomly for MomArm CheckUp to make
                % sure not always the right side is used.
                
                if checkAndAdaptMomArmsR && checkAndAdaptMomArmsL % both are available.
                    rndVal = randi([0, 1]); % random selection.
                    if rndVal == 0 % Right
                        sampleInput = sampleInputR;
                        checkAndAdaptMomArms = checkAndAdaptMomArmsR;
                    else
                        sampleInput = sampleInputL; % Left
                        checkAndAdaptMomArms = checkAndAdaptMomArmsL;
                    end
                elseif checkAndAdaptMomArmsR
                    sampleInput = sampleInputR; % Right
                    checkAndAdaptMomArms = checkAndAdaptMomArmsR;
                else
                    sampleInput = sampleInputL; % Left
                    checkAndAdaptMomArms = checkAndAdaptMomArmsL;
                end

                % Start to collect some data for subject info file.
                dataFile = fullfile(rootWorkingDirectory, 'data.xml');
                if isfile(dataFile)
                    persInfo = readstruct(dataFile);

                    % Make sure there are no field types that cannot be read in python.
                    persInfo = convertFieldsToChar(persInfo);

                    SubjInfo.MetaData = persInfo;
                    SubjInfo.MetaData.bodymassFromC3D = bodymass;
                    SubjInfo.MetaData.bodyheightFromC3D = BodyHeight;
                    SubjInfo.MetaData.description = 'This node contains information stored in a separate <data> file in the working directory. Note that if a variable such as AVR exists here, it does not necessarily mean it was used for personalization.';
                else
                    SubjInfo.MetaData.bodymassFromC3D = bodymass;
                    SubjInfo.MetaData.bodyheightFromC3D = BodyHeight;
                    SubjInfo.MetaData.description = 'There was no additional "data-file" available.';
                end

                % Add version number so one knows with which version
                % data were simulated.
                try
                    path2Version = path2setupFiles;

                    % Get the parent directory three levels above.
                    for i = 1:3
                        path2Version = fileparts(path2Version);
                    end

                    curVersion = readLastLine(fullfile(path2Version, 'version.md'));
                    SubjInfo.VersionProcessing = curVersion;
                end

                % Set flags to specify that yet not models were scaled.
                flag_modelScaled = false;

                for i = 1 : length(trials)

                    IC = InputData.(trials{i}).IC;
                    ICi = InputData.(trials{i}).ICi;
                    cIC = InputData.(trials{i}).cIC;
                    cTO = InputData.(trials{i}).cTO;
                    TO = InputData.(trials{i}).TO;
                    BW = InputData.(trials{i}).Bodymass * 9.81;
                    side = InputData.(trials{i}).Side;
                    path2c3d = (InputData.(trials{i}).c3dPath);
                    path2trc = (InputData.(trials{i}).trcPath);
                    path2mot = (InputData.(trials{i}).motPath);
                    path2enf = (InputData.(trials{i}).enfPath);
                    name = (InputData.(trials{i}).name);


                    % Prepare (scale) the models for both body sides.
                    if ~flag_modelScaled
                        [path2Model, torsionTool, tf_angle_fromSource_right, ...
                            tf_angle_fromSource_left, tf_angle_right, tf_angle_left, ...
                            MomArmsResolved, persInfoFromMarkers] = prepareScaledModel(rootWorkingDirectory, workingDirectory, path2GenericModels, path2opensim, bodymass, static, ...
                                                                    labFlag, lockSubtalar4Scaling, scaleMuscleStrength, manualMusScaleF, cPrefix, bodyheightGenericModel, BodyHeight, ...
                                                                    tf_angle_fromSource, tf_angle_r, tf_angle_l, strcat(rootWorkingDirectory, staticC3d), Model2Use, ...
                                                                    torsiontool, tib_torsion_RightMarkers, tib_torsion_LeftMarkers, ForceModelCreation, checkAndAdaptMomArms, sampleInput, useASTool, useStatic4TibRotEstimationAsFallback, useStatic4FrontAlignmentAsFallback);
    
                        % Add personalization info for fallback solutions
                        % just for documentation here.
                        SubjInfo.persInfoFromMarkers = persInfoFromMarkers;
                        
                        % Right info.
                        % Add info about model personalization to InputData
                        % logicals
                        SubjInfo.UsedData4Personalization.r.Logicals.TF_adaption = char(tf_angle_fromSource_right);
                        SubjInfo.UsedData4Personalization.r.Logicals.AV_adaption = char(torsionTool.femurAntetorsionAdaption);
                        SubjInfo.UsedData4Personalization.r.Logicals.TT_adaption = char(torsionTool.tibTorsionAdaption);
                        SubjInfo.UsedData4Personalization.r.Logicals.NS_adaption = char(torsionTool.neckShaftAdaption);
                        SubjInfo.UsedData4Personalization.r.Logicals.MomArmsResolved = MomArmsResolved;

                        % Values
                        SubjInfo.UsedData4Personalization.r.Values.TF_angle = tf_angle_right;
                        SubjInfo.UsedData4Personalization.r.Values.TT_angle = torsionTool.TTR;
                        SubjInfo.UsedData4Personalization.r.Values.NS_angle = torsionTool.NSAR;
                        SubjInfo.UsedData4Personalization.r.Values.AV_angle = torsionTool.AVR;

                        % Left info.
                        % Add info about model personalization to InputData
                        % logicals
                        SubjInfo.UsedData4Personalization.l.Logicals.TF_adaption = char(tf_angle_fromSource_left);
                        SubjInfo.UsedData4Personalization.l.Logicals.AV_adaption = char(torsionTool.femurAntetorsionAdaption);
                        SubjInfo.UsedData4Personalization.l.Logicals.TT_adaption = char(torsionTool.tibTorsionAdaption);
                        SubjInfo.UsedData4Personalization.l.Logicals.NS_adaption = char(torsionTool.neckShaftAdaption);
                        SubjInfo.UsedData4Personalization.l.Logicals.MomArmsResolved = MomArmsResolved;

                        % Values
                        SubjInfo.UsedData4Personalization.l.Values.TF_angle = tf_angle_left;
                        SubjInfo.UsedData4Personalization.l.Values.TT_angle = torsionTool.TTL;
                        SubjInfo.UsedData4Personalization.l.Values.NS_angle = torsionTool.NSAL;
                        SubjInfo.UsedData4Personalization.l.Values.AV_angle = torsionTool.AVL;

                        % Save Subject Info.
                        save(strcat(workingDirectory,'Simulation\', char(strrep(sub_name,' ','_')),'-',condition,'-SubjInfo.mat'),'SubjInfo');

                        % Set flag to know for next round that model was scaled already.
                        flag_modelScaled = true;
                    end

                    % Create external loads file
                    path2extLoad = createExtLoadsFile(path2enf, path2mot, name, firstContact_L, firstContact_R);

                    % Run LERNER pipeline for left or right side or full trial.
                    if strcmp(side, 'left')

                        % Run the pipeline, each trial in one cmd window
                        disp(char(strcat('>>>>> Started trial:',{' '},name, '_', side)));

                        % Run simulation if everything is ready
                        Model_workflow(workingDirectory, path2opensim, path2setupFiles, side, name, path2trc, path2mot, path2extLoad, path2Model, IC, cTO, cIC, TO, ICi, BW, cPrefix, tf_angle_r, tf_angle_l, timeNorm, labFlag, Model2Use,tasks);

                    elseif strcmp(side, 'right')

                        % Run the pipeline, each trial in one cmd window
                        disp(char(strcat('>>>>> Started trial:',{' '}, name, '_', side)));

                        % Run simulation if everything is ready
                        Model_workflow(workingDirectory, path2opensim, path2setupFiles, side, name, path2trc, path2mot, path2extLoad, path2Model, IC, cTO, cIC, TO, ICi, BW, cPrefix, tf_angle_r, tf_angle_l, timeNorm, labFlag, Model2Use,tasks);

                    else % Full trial processing

                        % Run the pipeline, each trial in one cmd window
                        disp(char(strcat('>>>>> Started trial:',{' '}, name, '_', side)));

                        % Run simulation if everything is ready
                        Model_workflow(workingDirectory, path2opensim, path2setupFiles, side, name, path2trc, path2mot, path2extLoad, path2Model, IC, cTO, cIC, TO, ICi, BW, cPrefix, tf_angle_r, tf_angle_l, timeNorm, labFlag, Model2Use,tasks);
                    end

                    %% Pause after maximum number of cmd windows are open
                    batchCount = batchCount + 1;
                    proFilesCount = proFilesCount + 1;

                    % Get processing time update
                    tEnd = seconds(toc(tStart));
                    tEnd.Format = 'hh:mm';

                    % If max. number of files reached or all files from condition forwarded display message and wait for next batch.
                    if batchCount == maxCmd
                        % Set the current number of processed files
                        % Display status
                        disp(' ');
                        disp(' ');
                        disp('*** Status Update ******************************************************************************************');
                        disp(char(strcat('Current working directory:', {' '}, workingDirectory)));
                        disp(char(strcat('Working directory',{' #'}, string(i_wd), {' '}, 'out of', {' #'}, string(length(workingDirectories)), '.')));
                        disp(char(strcat('Current condition:', {' '}, condition, {' out of { '}, strjoin(conditions, ', '), ' }.')));
                        if dataAugmentation; disp(char(strcat('Current data augmentation loop',{' #'}, string(i_dataAug), {' '}, 'out of', {' #'}, string(augN), '.'))); end
                        disp(strcat('Total files to process in current WD:..................', num2str(length(trials))));
                        disp(strcat('Processed files in current WD:.........................', num2str(proFilesCount)));
                        disp(strcat('Files waiting to be processed in current WD:...........', num2str(length(trials)-i)));
                        disp(strcat('Processed files in current batch:......................', num2str(batchCount)));
                        disp(strcat('Number of files forwarded to cmd window since start:...', num2str(numFiles)));
                        disp(strcat('Processing duration since <comak> start:',{' '}, string(tEnd),'(hh:mm).'));
                        disp('************************************************************************************************************');
                        disp(' ');
                        % To potentially prevent buffer overflow, make sure to close all fids in case fids are still open
                        fclose('all');
						
						if ~useCPUThreshold
							warning('Using CPU-based thresholding, since other solutions are not implemented yet');
						end
						
                        % Now check every X minutes if CPU load is below the threshold (e.g. 35%), before next batch of files.
                        if batchCount == maxCmd
                            disp(strcat('-> Matlab paused! When CPU-load drops below the threshold of', {' '}, string(thresholdCpuLoad) ,{'% '},'the next batch of files will kick off ...'));
                            CpuLoadBasedPausing_WIN11(thresholdCpuLoad, 60*1.5);

                            % Set batch cnt back to zero.
                            batchCount = 0; % Used to limit the number of open cmd windows
                        end

                        % Check how much RAM is left and take action if too low
                        freeRAMinPerc = measureFreeRAM;
                        if  freeRAMinPerc < thresholdFreeRAM % if below threshold then take action.
                            warning('Only %.0f%% free RAM left! I will wait for 1h and then will close all open CMD windows', freeRAMinPerc)
                            pause (60*60) % wait for 1h to finish open cmd windows
                            closeCmdWindows(); % now force to close open cmd windows
                            pause(120);
                        end

                        % Now measure again to make sure action was succesfull.
                        freeRAMinPerc = measureFreeRAM;
                        if  freeRAMinPerc < thresholdFreeRAM % if still below threshold then take action.
                            warning('Still only %.0f%% free RAM left after all CMD were closed! User action required! Free up RAM, e.g. by restarting Matlab.', freeRAMinPerc)

                            % In case user allowed for autorestart, set flag to kickoff autorestart at end of loop.
                            if allowAutoRestart
                                kickOffRestart = true;
                            else
                                pause();
                            end
                        end

                        % Make sure that cmd windows that do not converge get closed after a while. This hopefully prevents any "buffer overflow" errors.
                        if numFiles > overFlowThreshold
                            disp(['>>>>>  After ', num2str(overFlowThreshold) ,' files I will take a 2h break. Afterwards all remaining cmd windows will be terminated.']);
                            pause(60*60*2);
                            closeCmdWindows(); % now force to close open cmd windows
                            pause(60);

                            % Adjust overFlowThreshold
                            overFlowThreshold = overFlowThreshold + overFlowThreshold;
                        end
                    end
                    numFiles = numFiles + 1;
                end

            catch ME
                if catchErrors
                    disp(char(strcat('>>>>> An error occured during processing of condition', {' "'}, condition, {'" '},'in wd:', {' '},  workingDirectory,{'. '}, 'Skipped condition in current wd!')));
                    fprintf('>>>>> Error message: %s\n', ME.message);
                    errN = errN + 1;
                    failures(errN).cond = condition;
                    failures(errN).wd = workingDirectory;
                    failures(errN).err = getReport(ME);
                    continue % Jump to next iteration of: for i
                else
                    rethrow(ME)
                end
            end
        end

        %% Now run postprocessing if user selected this option
        if performPostProcessing

            % Wait to make sure all batch files are finished
            disp('>>>>> Waiting to postprocess the current batch of files ...');
            thresholdCpuLoadPostProcessing = 10; % Define a threshold the CPU load has to fall below, default ~ 10%
            CpuLoadBasedPausing(thresholdCpuLoadPostProcessing, 60*3) % check every X minutes if CPU load is below the threshold (e.g. 35%). Wait extra long here.
            %pause

            % Make sure the condition/prefix/workingDirectory are cell arrays
            Condition4PostPro = {condition};
            Prefix4PostP = {cPrefix};
            workingDirectory4PostPro = {workingDirectories{i_wd}}; % use the wd without the added pipeline here; this will be done later.

            % Run the post processing loop
            loops4PostProcessingVarModels(rootDirectory, workingDirectory4PostPro, Condition4PostPro, path2setupFiles, Prefix4PostP, timeNormFlag, catchErrors, trialType, Model2Use);

            % Final postprocessing message
            disp(strcat('*********************************  MODEL postprocessing finished for <',workingDirectory, 'current trial  *********************************'));
        end

    catch ME
        if catchErrors
            disp(char(strcat('>>>>> An error occured during processing of wd:', {' '},  workingDirectory,{'. '}, 'Skipped wd!')));
            fprintf('>>>>> Error message: %s\n', ME.message);
            errN = errN + 1;
            failures(errN).cond = condition;
            failures(errN).wd = workingDirectory;
            failures(errN).err = getReport(ME);
            continue % Jump to next iteration of: for i
        else
            rethrow(ME)
        end
    end
end

% Update status in local WorkingDirectories file located in the
% rootDirectory. Note that it is important to update and save the file
% during each WD loop. In case of a Matlab crash otherwise all info
% would be lost.

maxAttempts = 5; % Set the maximum number of attempts
currentAttempt = 1; % Initialize attempt counter

while currentAttempt <= maxAttempts
    try
        wdFile = fullfile(rootDirectory, 'workingDirectories.xlsx');
        data = readtable(wdFile);
        c_idx = find(strcmp(data.WorkingDirectories, rootWorkingDirectory));
        data.Simulation(c_idx) = 1;
        data.ModelWOadaption(c_idx) = {MomArmsResolved};
        writetable(data, wdFile);
        break;
    catch
        % Wait for some time before the next attempt
        pause(5);

        % Increment the attempt counter
        currentAttempt = currentAttempt + 1;
    end
end

% In case user selected the option that matlab can restart automatically in
% case of low memory.
if allowAutoRestart

    if kickOffRestart

        % Increase i_wd by one so when script reloads it willstart with next loop.
        i_wd = i_wd + 1;

        % Set restartFlag back to false to ensure it doesn't restart again.
        kickOffRestart = false; % Update the restart flag.

        % Save the current workspace to a file
        pause (60)
        warning('Saving workspace and restarting MATLAB!');
        save(checkPointFilePath);

        % Construct the command to restart MATLAB and run the same script
        matlab_command = ['matlab -r "run(''' mfilename('fullpath') ''');" &'];
        system(matlab_command);

        % Exit the current MATLAB session
        exit;
    end
end

% Save errors log as *.txt in root dir.
writetable(struct2table(failures), fullfile(rootDirectory,'MatlabErrorLog.xlsx'));

% End timer
tEnd = seconds(toc(tStart));
tEnd.Format = 'hh:mm';
disp(''); % new line
disp('************************************************************************************************************');
disp(strcat('>>>>>  Total number of matlab errors catched:',{' '}, string(errN), '. <<<<<'));
disp(strcat('>>>>>  Total processing duration:',{' '}, string(tEnd),'(hh:mm). <<<<<'));

