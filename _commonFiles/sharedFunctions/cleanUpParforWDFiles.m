function cleanUpParforWDFiles(folder, numFiles, outputFileName)
    % This function concatenates Excel files from parforLoop_1.xlsx to parforLoop_N.xlsx
    % into a single Excel file and deletes the original files after processing.
    %
    % Parameters:
    % folder (string): The directory where the parforLoop files are located
    % numFiles (integer): The number of parforLoop files to process (e.g., 100)
    % outputFileName (string): The name of the output file (e.g., 'CombinedFile.xlsx')

    % Initialize a table to store all data.
    combinedData = [];

    % Loop through the files from parforLoop_1 to parforLoop_numFiles.
    for i = 1:numFiles
        
        % Create the filename.
        fileName = fullfile(folder, sprintf('parforLoop_%d.xlsx', i));
        
        % Check if the file exists.
        if isfile(fileName)
            % Read the Excel file.
            data = readtable(fileName);
            
            % Concatenate the data into combinedData.
            combinedData = [combinedData; data]; %#ok<AGROW>
            
            % Delete the file after reading.
            delete(fileName);
        else
            warning('File %s does not exist', fileName);
        end
    end

    % Write the combined data to an Excel file.
    writetable(combinedData, outputFileName);

    disp('Concatenation of WD log files completed and individual files deleted.');
end
