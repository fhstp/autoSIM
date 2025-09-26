function loops4PostProcessingVarModels(rootDirectory, workingDirectories, conditions, path2setupFiles, cPrefix, timeNormFlag, catchErrors, trialType, Model2Use, i_wd)

%% Run pipelines
disp('>>>>> Starting to postprocess the data ...');

% Init error log
errN = 0; % Number of errors occured and catched

% Add loops here:
% Loop through conditions
for cond = 1 : length(conditions)

    try
        % Set current condition
        condition = conditions{cond};

        % Set current workingDirectory
        workingDirectory = char(fullfile(strcat(workingDirectories{i_wd}, Model2Use, '-', cPrefix,'\')));

        % Create error reports for each trial and save output in working directory
        disp('>>>>> Started creating error and per trial reports <<');
        createErrorReportsPerTrial(workingDirectory, rootDirectory, path2setupFiles, condition, cPrefix, trialType, Model2Use);

        % Create a data report of all (selected) trials within one working
        % directory and create summarized output files and plots.
        disp('>>>>> Started creating summarized trials reports <<');
        createTrialOverview4oneSubject(workingDirectory, rootDirectory, path2setupFiles, condition, cPrefix, timeNormFlag, trialType, Model2Use);

    catch ME
        if catchErrors
            disp(char(strcat('>>>>> An error occured during processing of condition', {' "'}, condition, {'" '},'with prefix:', {' "'}, cPrefix, {'"'}, ', in wd:', {' '},  workingDirectory,{'. '}, 'Skipped prefix/condition in current wd!')));
            fprintf('>>>>> Error message: %s\n', ME.message);
            errN = errN + 1;
        else
            rethrow(ME)
        end
    end
end
% Update status in local WorkingDirectories file located in the
% rootDirectory. Note that it is important to update and save the file
% during each loop. In case of a Matlab crash otherwise all info
% would be lost.

% For the parfor version of this script we need to save each time a
% separate file, since otherwise Matlab might crash. The files will be
% concatenated after the parfor loop is finished.

% Create a directory to store the individual files
outputDir = fullfile(rootDirectory, strcat('_', Model2Use, '-WD-tmpLogFiles')); % Note this path name is also hardcoded in startVarModels.m
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Define the header of the excel file.
headers = {'Index', 'WorkingDirectories', 'Statics', 'Simulation', 'Processed', 'ModelWOadaption'};

% Create a table with the data
data = table(i_wd, string(workingDirectories{i_wd}), NaN, NaN, 1, NaN, 'VariableNames', headers);

% Use parfor to create individual files
wdFile = fullfile(outputDir, sprintf('postPro_parforLoop_%d.xlsx', i_wd)); % Note this file name is also hardcoded in cleanUpParforWDFiles4PostPro.m
writetable(data, wdFile);

% ---- deprecated

% maxAttempts = 5; % Set the maximum number of attempts
% currentAttempt = 1; % Initialize attempt counter

% while currentAttempt <= maxAttempts
%     try
%         % Try to write the information in the prefix specific WD file.
%         % If not available write data to "standard" WD file.
%         wdFile_v1 = fullfile(rootDirectory, strcat('workingDirectories-', Model2Use, '_', prefix, '.xlsx'));
%         wdFile_v2 = fullfile(rootDirectory, 'workingDirectories.xlsx');
%         if isfile(wdFile_v1)
%             wdFile = char(wdFile_v1);
%         else
%             wdFile = char(wdFile_v2);
%         end
%
%         % Now write the data.
%         data = readtable(wdFile);
%         c_idx = find(strcmp(data.WorkingDirectories, workingDirectories{i_wd}));
%         data.Processed(c_idx) = 1;
%         writetable(data, wdFile);
%         break;
%     catch
%         % Wait for some time before the next attempt
%         pause(5);
%
%         % Increment the attempt counter
%         currentAttempt = currentAttempt + 1;
%     end
% end