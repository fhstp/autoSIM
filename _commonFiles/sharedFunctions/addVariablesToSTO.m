function addVariablesToSTO(stoFile, newVarNames, newVarData)
% addVariablesToSTO Reads an OpenSim .sto file, adds new variables,
% and writes the modified file back to disk.
%
% INPUTS:
%   stoFile   - path to the .sto file
%   newVarNames - cell array of variable names, e.g. {'Var1','Var2','Var3','Var4'}
%   newVarData  - matrix with same number of rows as the STO file,
%                 columns correspond to newVarNames
%
% Example:
%   data = rand(100,4); % 4 new variables, 100 rows
%   addVariablesToSTO('input.sto','output.sto',...
%       {'NewVar1','NewVar2','NewVar3','NewVar4'}, data)
%
% Author: B. Horsak, 09/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For dev only.
% stoFile_out = [stoFile(1:end-4), '_changed_.sto'];

%--- Step 1: Read the .sto file ---
raw = fileread(stoFile);
lines = strsplit(raw,'\n');

% Find header end (line with "endheader")
headerEndIdx = find(contains(lines,'endheader'),1);
header = lines(1:headerEndIdx);

% Read numeric data into a table
data = readtable(stoFile,'FileType','text','Delimiter','\t','HeaderLines',headerEndIdx);

%--- Step 2: Validate input sizes ---
if height(data) ~= size(newVarData,1)
    error('Row count mismatch: newVarData must have same number of rows as STO file.');
end
if length(newVarNames) ~= size(newVarData,2)
    error('Number of variable names must equal number of columns in newVarData.');
end

%--- Step 3: Add new variables ---
for k = 1:length(newVarNames)
    data.(newVarNames{k}) = newVarData(:,k);
end

%--- Step 4: Update header nColumns ---
colLineIdx = find(startsWith(strtrim(header),'nColumns'),1);
if ~isempty(colLineIdx)
    newCols = width(data);
    header{colLineIdx} = sprintf('nColumns=%d', newCols);
end

%--- Step 5: Write back to .sto file ---
fid = fopen(stoFile,'w');
if fid == -1
    error('Cannot open output file: %s', outputFile);
end

% Write header
for i = 1:numel(header)
    fprintf(fid,'%s\n',header{i});
end

% Write column names
varNames = data.Properties.VariableNames;
fprintf(fid,'%s\t',varNames{1:end-1});
fprintf(fid,'%s\n',varNames{end});

% Write data row by row
fmt = [repmat('%g\t',1,width(data)-1), '%g\n'];
for r = 1:height(data)
    fprintf(fid, fmt, table2array(data(r,:)));
end

fclose(fid);
end
