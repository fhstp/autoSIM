function [out, usedFileList] = getAllVariables(tmpList, dataRedFlag, files2skip, thresholdS, maxIKcrit, rmsIKcrit)
% This file reads the JAM overall mat file, creates a similar struct, but
% it creates for each variable a table populating column-wise the data of
% the entire sample. For this purpose the data are aggregated subject-wise
% i.e. by calculating a mean average curve.
%
% Input:
% tmpList: file list
% dataRedFlag: 'data reduction Flag', specifies the method used to redcue
%               the data, e.g. mean
%
%
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         04/2022
% -------------------------------------------------------------------------

% Number of subjects
N = length(tmpList);

% Initialize empty struct
tmpFile = horzcat(tmpList(1).folder,'\', tmpList(1).name);
tmpDat = load(tmpFile);
fnR = fieldnames(tmpDat.subResults.perSubject.r);
fnL = fieldnames(tmpDat.subResults.perSubject.l);

% Initialize empty UsedFileLiest
usedFileList = {};

% Remove cell arrays that can't be averaged
vars2exc = {'IKerrors','trialName', 'cTO','TO','IC','cIC','ICi','Bodymass'};
for i = 1 : length(vars2exc)
    fnR(strcmp(fnR,vars2exc{i})) = [];
    fnL(strcmp(fnL,vars2exc{i})) = [];
end
cL = cell(length(fnL),1);
out.l = cell2struct(cL,fnL);

cR = cell(length(fnR),1);
out.r = cell2struct(cR,fnR);

cntUsedFiles = 1; % to count number of files used
% Loop per subject
for i_subj = 1 : N %

    % Get file and data per Subject
    tmpFile = horzcat(tmpList(i_subj).folder,'\', tmpList(i_subj).name);
    tmpDat = load(tmpFile);
    [~,name,ext] = fileparts(tmpFile);
    trialNameLong = [name,ext];
    trialNameShort = strrep(trialNameLong, '-JAM-Results-All-Trials.mat', '');


    % Loop for side
    side = {'r','l'};
    for i_side = 1 : 2
        sideTmp = side{i_side};
        fnSide = fieldnames(tmpDat.subResults.perSubject.(sideTmp));

        % Check if file should be skipped or not
        if contains(trialNameShort,files2skip) && ~isempty(files2skip)
            disp(char(strcat('>>>>> User selected to skip trial:', {' '}, trialNameShort)))

        elseif ~sum(strcmp(fieldnames(tmpDat.subResults.perSubject.(sideTmp).ContactForces), strcat('time_',sideTmp))) == 1 % check if comak data exist
            disp(char(strcat('>>>>> ', trialNameLong, {' '}, 'Attention: no data found for comak results on <', sideTmp,'> side. Skipped trial!')));

        else % or proceed

            % Create vars 2 check for IK errors
            vars2CheckR = {'RKNE', 'RANK', 'RHEE', 'RS1', 'RS2', 'RS3', 'RT1', 'RT2', 'RT3'};
            vars2CheckL = {'LKNE', 'LANK', 'LHEE', 'LS1', 'LS2', 'LS3', 'LT1', 'LT2', 'LT3'};
            vars2CheckAlways = {'SACR', 'RASI', 'LASI'};

            if strcmp(sideTmp, 'l')
                vars2check = [vars2CheckL, vars2CheckAlways];
            elseif strcmp(sideTmp, 'r')
                vars2check = [vars2CheckR, vars2CheckAlways];
            end

            % Check which markers are within the IK error threshold
            check_rms = [];
            check_max = [];
            for s = 1 : length(vars2check)
                check_rms(s,:) = logical(rms(tmpDat.subResults.perSubject.(sideTmp).IKerrors.(vars2check{s})) <= rmsIKcrit);
                check_max(s,:) = logical(max(tmpDat.subResults.perSubject.(sideTmp).IKerrors.(vars2check{s})) <= maxIKcrit);
            end

            % Calculate knock-out crit
            idxPassedIK = all(check_rms,1) & all(check_max,1); 

            % Check if there are trials which meet IKerror criteria
            if sum(idxPassedIK) == 0 
                disp(char(strcat('>>>>> ', trialNameLong, {' '}, 'Attention: data are out of IKerror bounds on <', sideTmp,'> side. Skipped trial!')));
            
            else % or proceed
                % Find all single trial names
                R = tmpDat.subResults.perSubject.(sideTmp).trialName;
                SingleTrialNames = R(~cellfun('isempty',R));
                SingleTrialNames = SingleTrialNames(idxPassedIK);

                % Create variable list to check if comak results exist
                if strcmp(sideTmp, 'r')
                    comakVars2Check = {...
                        'tf_contact_casting_total_contact_force_x_r', ...
                        'tf_contact_casting_total_contact_force_y_r', ...
                        'tf_contact_casting_total_contact_force_z_r'};
                elseif strcmp(sideTmp, 'l')

                    comakVars2Check = {...
                        'tf_contact_casting_total_contact_force_x_l', ...
                        'tf_contact_casting_total_contact_force_y_l', ...
                        'tf_contact_casting_total_contact_force_z_l'};
                end

                % Create the most representative trial index | Array = nStride x nEpoch x nVariable
                if strcmp(sideTmp, 'r')
                    vars4Sangeux = {...
                        'pelvis_tilt', ...
                        'pelvis_list', ...
                        'pelvis_rot', ...
                        'hip_flex_r_ipsilateral', ...
                        'hip_add_r_ipsilateral', ...
                        'hip_rot_r_ipsilateral', ...
                        'pf_flex_r_ipsilateral', ...
                        'pf_rot_r_ipsilateral', ...
                        'pf_tilt_r_ipsilateral', ...
                        'knee_flex_r_ipsilateral', ...
                        'knee_add_r_ipsilateral', ...
                        'knee_rot_r_ipsilateral', ...
                        'ankle_flex_r_ipsilateral'};
                
                elseif strcmp(sideTmp, 'l')
                
                    vars4Sangeux = {...
                        'pelvis_tilt', ...
                        'pelvis_list', ...
                        'pelvis_rot', ...
                        'hip_flex_l_ipsilateral', ...
                        'hip_add_l_ipsilateral', ...
                        'hip_rot_l_ipsilateral', ...
                        'pf_flex_l_ipsilateral', ...
                        'pf_rot_l_ipsilateral', ...
                        'pf_tilt_l_ipsilateral', ...
                        'knee_flex_l_ipsilateral', ...
                        'knee_add_l_ipsilateral', ...
                        'knee_rot_l_ipsilateral', ...
                        'ankle_flex_l_ipsilateral'};
                end

                % Check if field exists
                checkFn = 0;
                for i_check = 1 : length(comakVars2Check)
                    checkFn(i_check) = isfield(tmpDat.subResults.perSubject.(sideTmp).ContactForces, comakVars2Check{i_check});
                end

                % Check if there are data
                checkData = 0;
                if sum(checkFn) == length(comakVars2Check) % do that only when fields exist
                    for i_check = 1 : length(comakVars2Check)
                        checkData(i_check) = ~isempty(tmpDat.subResults.perSubject.(sideTmp).ContactForces.(comakVars2Check{i_check})(:,idxPassedIK));
                    end
                end

                % Only get data if we have comak results
                if sum(checkFn) ~= length(comakVars2Check) && sum(checkData) ~= length(checkData)
                    disp(char(strcat('>>>>> ', trialNameLong, {' '}, 'Attention: no data found for comak results. Skipped trial!')));
                else
                    
                    % Create most representative trial array
                    arraySangeuxR = [];
                    for j = 1 : length(vars4Sangeux)
                        arraySangeuxR = cat(3, arraySangeuxR, tmpDat.subResults.perSubject.(sideTmp).InverseKinematics.(vars4Sangeux{j})(:,idxPassedIK)');
                    end

                    % Run the method on array
                    [~, i_D, Outliers, ~, ~] = MultiFunDepth(arraySangeuxR, thresholdS);

                    % Store the info about which file was used in second output
                    switch dataRedFlag
                        case 'avg-outlierReduced'
                            usedFileList(cntUsedFiles,1) = {trialNameShort};
                            usedFileList(cntUsedFiles,2:length(SingleTrialNames(:, ~Outliers))+1) = SingleTrialNames(:, ~Outliers);
                        case 'avg-all-omitnan'
                            usedFileList(cntUsedFiles,1) = {trialNameShort};
                            usedFileList(cntUsedFiles,2:length(SingleTrialNames(:))+1) = SingleTrialNames(:)';
                        case 'firstTrialOnly'
                            usedFileList(cntUsedFiles,1) = {trialNameShort};
                            usedFileList(cntUsedFiles,2) = SingleTrialNames(1);
                        case 'mostRepTrial'
                            usedFileList(cntUsedFiles,1) = {trialNameShort};
                            usedFileList(cntUsedFiles,2) = SingleTrialNames(i_D(1));
                    end

                    cntUsedFiles = cntUsedFiles + 1;

                    % Loop for Variable Category (ContacForce, ContactPressure, etc.)
                    for i_vCat = 1 : length(fnSide)
                        % current fieldname
                        fn_vCat = fnSide{i_vCat};

                        if sum(strcmp(fn_vCat, {'trialName', 'IKerrors'})) == 0 % exclude trial name, IKerrors

                            % Loop per variable
                            for i_Var = 1 : length(fieldnames(tmpDat.subResults.perSubject.(sideTmp).(fn_vCat)))

                                % current fieldname
                                fn_vars = fieldnames(tmpDat.subResults.perSubject.(sideTmp).(fn_vCat));
                                fn_cVar = fn_vars{i_Var};

                                % Reduce data to IK passed
                                tmpVar = tmpDat.subResults.perSubject.(sideTmp).(fn_vCat).(fn_cVar)(:,idxPassedIK);

                                % omit time var because avg does not make sense here
                                if ~contains(fn_cVar, 'time')

                                    % initialize struct field if it does not exist
                                    if ~isfield(out.(sideTmp).(fn_vCat), fn_cVar)
                                        out.(sideTmp).(fn_vCat).(fn_cVar) = table();
                                    end

                                    % Aggregate the data for each subject
                                    if ~isempty(tmpVar) % Check if it consists data
                                        switch dataRedFlag
                                            case 'avg-outlierReduced'
                                                Outliers = transpose(Outliers);
                                                tmpVar = tmpVar(:,~Outliers);
                                                out.(sideTmp).(fn_vCat).(fn_cVar).(trialNameShort) = mean(tmpVar,2, 'omitnan');

                                            case 'avg-all-omitnan'
                                                out.(sideTmp).(fn_vCat).(fn_cVar).(trialNameShort) = mean(tmpVar,2, 'omitnan');

                                            case 'firstTrialOnly'
                                                out.(sideTmp).(fn_vCat).(fn_cVar).(trialNameShort) = tmpVar(:,1);

                                            case 'mostRepTrial'
                                                if size(tmpVar,2) < i_D(1)
                                                    disp(char(strcat('>>>>> ', trialNameLong, 'Warning: the index of the most rep. trial does not match with the size of', {' '}, '<', fn_cVar,'>', {' '}, 'Last index used! Please inspect data!')));
                                                    out.(sideTmp).(fn_vCat).(fn_cVar).(trialNameShort) = tmpVar(:, size(tmpVar,2));
                                                else
                                                    out.(sideTmp).(fn_vCat).(fn_cVar).(trialNameShort) = tmpVar(:, i_D(1));
                                                end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end