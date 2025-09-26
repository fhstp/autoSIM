function [out1, out2] = getHeaderAlignedDataSet(var1, var2)
% This function expects to have the sub ID before the first '-' in the col
% header string. It then alignes both data to have only the shared ID in
% the output, and both outputs in the same order.

%% Get var names
var1Names = var1.Properties.VariableNames;
var2Names = var2.Properties.VariableNames;

%% Reduce to ID information for var1 & var2
% Var 1
for i = 1 : length(var1Names)
    varTmp = var1Names{i};
    idx = strfind(varTmp,'-');
    if ~isempty(idx)
    idx = idx(1)-1; % get rid of '-'
    var1Names{i} = varTmp(1:idx);
    end
end

% Var2
for i = 1 : length(var2Names)
    varTmp = var2Names{i};
    idx = strfind(varTmp,'-');
    if ~isempty(idx)
    idx = idx(1)-1; % get rid of '-'
    var2Names{i} = varTmp(1:idx);
    end
end

%% Now deal both vars
[C, ~, ~] = intersect(var1Names,var2Names, 'stable');
idx_1 = ~ismember(var1Names, C);
idx_2 = ~ismember(var2Names, C);
out1 = removevars(var1, idx_1);
out2 = removevars(var2, idx_2);

%% Now check that varnames are aligned in both variables
ID1_2Check = out1.Properties.VariableNames;
ID2_2Check = out2.Properties.VariableNames;

% Reduce to ID information for var1 & var2
% Var 1
for i = 1 : length(ID1_2Check)
    varTmp = ID1_2Check{i};
    idx = strfind(varTmp,'-');
    if ~isempty(idx)
    idx = idx(1)-1; % get rid of '-'
    ID1_2Check{i} = varTmp(1:idx);
    end
end

% Var2
for i = 1 : length(ID2_2Check)
    varTmp = ID2_2Check{i};
    idx = strfind(varTmp,'-');
    if ~isempty(idx)
    idx = idx(1)-1; % get rid of '-'
    ID2_2Check{i} = varTmp(1:idx);
    end
end

if ~isequal(ID1_2Check, ID2_2Check)
    error('>>>>> WARNING: variable names are not identical!');
end