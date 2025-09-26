function [out, fig] = inspectSingleVariable(tmpList, dataRedFlag, cleanOutlierFlag, node1, node2, variable)
% This file reads the JAM overall mat file, and  creates for a selected
% variable a table populating columns-wise the data of
% the entire sample. For this purpose the data are aggregated subject-wise
% i.e. by calculating a mean avergae curve.
%
% Input: 
% tmpList: file list
% dataRedFlag: 'data reduction Flag', secifies the method used to redcue
%               the data, e.g. mean
% 
% cleanOutlierFlag: specifies if and by which method the data should be
%                   clean variable-wise from outliers before aggregating
% node1, node2, variable: the strcut nodes necessayr to navigate to the
%                         data
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         02/2022
% -------------------------------------------------------------------------
fig = figure;
N = length(tmpList); % Number of subjects
out = table();
% Loop per subject
for i = 1 : N

    % Get file and data
    tmpFile = horzcat(tmpList(i).folder,'\', tmpList(i).name);
    tmpDat = load(tmpFile);
    tmpVar = tmpDat.subResults.(node1).(node2).(variable);
    [~,name,ext] = fileparts(tmpFile);
    trialNameLong = [name,ext];
    trialNameShort = strrep(trialNameLong, '-Dynamic-JAM-Results-All-Trials.mat', '');

    % Get some subject specific infos
    idx = strfind(variable, '_');
    side = variable(idx(end):end); % find end tag defining the side
    tmpTrials = tmpDat.subResults.perSubject.trialName(contains(tmpDat.subResults.perSubject.trialName,side)); % get only trials for specific side

    % Clean or not clean the data, that's the question ...
    switch cleanOutlierFlag
        case 'false'
            % Plot consistency for each subjects
            nCols = 5;
            if nCols > N; nCols = N; end
            nRows = ceil(N/nCols);
            subplot(nRows,nCols,i);
            plot(tmpVar, 'LineWidth', 1.2); % plot everything except outliers
            sgt = sgtitle({'Consistency plots for all subjects:', variable});
            sgt.Interpreter = 'none';
            lgd = legend(tmpTrials);
            lgd.FontSize = 5;
            lgd.Interpreter = 'none';

        case 'multiFun'
            % This function does not like nans
            if any(isnan(tmpVar(:)))
                disp(strcat('>> Attention: NaNs detected and column-wise removed in:  <', trialNameLong, '> for var:','<', variable,'>', 'Please inspect data!'))
                tmpVar = tmpVar(:,any(~isnan(tmpVar)));
                tmpTrials = tmpTrials(any(~isnan(tmpVar))); % make sure to skip also in the trials names
            end
            Array = transpose(tmpVar);
            [D, i_D, Outliers, LRT, F] = MultiFunDepth(Array, 3);
            Outliers = transpose(Outliers);
            tmpVarOrig = tmpVar; % var including the outliers for plotting
            tmpVar = tmpVar(:,~Outliers);
            index = Outliers;

            % Plot consistency for each subjects and highlight potential outliers
            nCols = 5;
            if nCols > N; nCols = N; end
            nRows = ceil(N/nCols);
            subplot(nRows,nCols,i);
            plot(tmpVarOrig(:,~index), 'LineWidth', 1.2); % plot everything except outliers
            hold on
            plot(tmpVarOrig(:,index),'--r', 'LineWidth', 2); % plot outliers in diferent color
            sgt = sgtitle({'Consistency plots for all subjects:', variable});
            sgt.Interpreter = 'none';
            lgd = legend([tmpTrials(~index), tmpTrials(index)]);
            lgd.FontSize = 5;
            lgd.Interpreter = 'none';
    end



    % Aggregate the data for each subject
    switch dataRedFlag
        case 'avg'
            out.(trialNameShort) = mean(tmpVar,2, 'omitnan');
        case 'all'
            % -- Add something here ---
    end
end

