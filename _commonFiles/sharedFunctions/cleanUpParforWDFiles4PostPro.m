function cleanUpParforWDFiles4PostPro(folder, numFiles, outputFileName)
% This function processes Excel files postPro_parforLoop_1.xlsx to postPro_parforLoop_N.xlsx.
% If outputFileName exists, it updates the "Processed" column in-place using values from the input files.
% If it doesn't exist, it concatenates and saves all input files into a new output file.
%
% Parameters:
% folder (string): Directory with input Excel files
% numFiles (integer): Number of files to process
% outputFileName (string): Final combined or updated Excel file name

combinedData = [];

% Check if the output file already exists
if isfile(outputFileName)
    existingData = readtable(outputFileName);
    outputExists = true;
else
    outputExists = false;
end

for i = 1:numFiles
    % Generate the filename for the current file
    fileName = fullfile(folder, sprintf('postPro_parforLoop_%d.xlsx', i));

    if isfile(fileName)
        % Read the current file
        data = readtable(fileName);

        if outputExists
            % Verify required columns exist
            if all(ismember({'WorkingDirectories', 'Processed'}, existingData.Properties.VariableNames)) && ...
                    all(ismember({'WorkingDirectories', 'Processed'}, data.Properties.VariableNames))

                wd = data.WorkingDirectories{1};  % Assuming cell array of strings
                procValue = data.Processed(1);

                % Match in existing data
                matchIdx = find(strcmp(existingData.WorkingDirectories, wd));

                if ~isempty(matchIdx)
                    existingData.Processed(matchIdx) = procValue;
                else
                    warning('No matching WorkingDirectory found in existing file for: %s', wd);
                end
            else
                warning('Missing "WorkingDirectories" or "Processed" column in either file.');
            end
        else
            % Accumulate new data for concatenation
            combinedData = [combinedData; data]; %#ok<AGROW>
        end

        % Optionally delete file after processing
        delete(fileName);
    else
        warning('File %s does not exist', fileName);
    end
end

% Save the results
if outputExists
    writetable(existingData, outputFileName);
    disp('Updated existing output file with new "Processed" values.');
else
    writetable(combinedData, outputFileName);
    disp('Created new combined output file from individual WD files.');
end
end
