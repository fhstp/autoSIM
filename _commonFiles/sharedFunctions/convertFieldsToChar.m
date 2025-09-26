function structWithChar = convertFieldsToChar(inputStruct)
% Function to convert string fields or datetime to char

    structWithChar = inputStruct; % Create a copy of the original structure

    % Iterate through the field names
    fieldNames = fieldnames(structWithChar);
    for i = 1:length(fieldNames)
        % Check if the field contains a string, or datetime and convert to char if true
        if isstring(structWithChar.(fieldNames{i})) || isdatetime(structWithChar.(fieldNames{i}))
            structWithChar.(fieldNames{i}) = {char(structWithChar.(fieldNames{i}))};
        end
    end
end