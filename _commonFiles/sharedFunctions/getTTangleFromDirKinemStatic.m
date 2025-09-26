function [tt_angle_used] = getTTangleFromDirKinemStatic(static, InputData, side, varNameKneeAngle_c3d)
%
%   This function estimates the tibial torsion angle based on knee rotation
%   during gait, extracted from the direct kinematic outputs of the 
% CleveLand Markerset model from the OSS. It computes either:
%     - the median of multiple dynamic trials ('fromDynamic'), or
%     - a single static trial ('fromStatic').
%
%   INPUTS:
%       static     - string path to a static C3D file (used if 'fromStatic' is selected)
%       InputData  - struct containing trial information with fields including 'c3dPath'
%       side       - character ('L' or 'R') indicating which leg to analyze
%
%   OUTPUT:
%       tt_angle_used - estimated tibial torsion angle in degrees (external rotation positive)
%
%   METHOD:
%       - 'fromDynamic' (default): 
%           Extracts knee rotation from the stance phase (foot strike to foot off)
%           of each dynamic trial and computes the median across steps and files.
%       - 'fromStatic': 
%           Extracts the median knee rotation from the static trial directly.
%
%   NOTE:
%       The rotation is multiplied by -1 to match the convention of the TorsionTool
%       where external rotation is considered positive.
%
%   REQUIREMENTS:
%       - BTK Toolbox (for reading C3D files and extracting events/angles)
%
% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Last changed:         06/2025
% -------------------------------------------------------------------------

setMethod = 'fromDynamic'; % default = 'fromDynamic'; OR 'fromStatic'

% Using the dynamic trials is based on the publication from Radler et al. (2010) we will use 'fromDynamic' as default, see Radler, et al., 2010. Torsional profile versus gait analysis: consistency between the anatomic torsion and the resulting gait pattern in patients with rotational malalignment of the lower extremity. Gait Posture 32, 405–410. https://doi.org/10.1016/j.gaitpost.2010.06.019


switch setMethod
    case 'fromDynamic'

        %% Get TT angle median direct from direct kinematic model outputs from the dynamic *.c3d trials.
        % Get all c3d file paths.
        fn = fieldnames(InputData);
        allC3DPaths = {};

        for i_fn = 1:length(fn)
            c_fn = fn{i_fn};
            if isfield(InputData.(c_fn), 'c3dPath')
                allC3DPaths{end+1} = InputData.(c_fn).c3dPath; %#ok<SAGROW>
            else
                warning('No c3dPath field in InputData.%s', c_fn);
            end
        end

        uniqueC3DPaths = unique(allC3DPaths);

        % Initialize output storage.
        kneeRots = NaN(length(uniqueC3DPaths), 1);

        % Determine correct variables and event names
        if lower(side) == 'l'
            kneeVariable = varNameKneeAngle_c3d.L;
            strikeField = "Left_Foot_Strike";
            footOffField = "Left_Foot_Off";
        elseif lower(side) == 'r'
            kneeVariable = varNameKneeAngle_c3d.R;
            strikeField = "Right_Foot_Strike";
            footOffField = "Right_Foot_Off";
        else
            error("The variable 'side' is not defined correctly in <getTTangleFromDirKinemStatic>")
        end

        % Loop through trials
        for i = 1:length(uniqueC3DPaths)
            c_File = uniqueC3DPaths{i};

            if ~isfile(c_File)
                warning('File not found: %s', c_File);
                continue;
            end

            try
                acq = btkReadAcquisition(c_File);
                angles = btkGetAngles(acq);
                events = btkGetEvents(acq);
                metadata = btkGetMetaData(acq);
                framerate = metadata.children.TRIAL.children.CAMERA_RATE.info.values;
                
                % Close acq
                btkCloseAcquisition(acq);

                if isfield(angles, kneeVariable) && isfield(events, strikeField) && isfield(events, footOffField)
                    strikes = events.(strikeField);
                    footOffs = events.(footOffField);

                    % Make sure we only have foot offs that are later than the first foot strike.
                    footOffs = footOffs(footOffs > strikes(1));

                    steps_N = length(strikes) - 1;
                    kneeStepRot = NaN(steps_N, 1);

                    % Get values for all gait cycles (only stance) of each c3d file.
                    for j = 1:steps_N
                        start = round(strikes(j) / (1/framerate));
                        stop  = round(footOffs(j) / (1/framerate));

                        tmp = angles.(kneeVariable)(:, varNameKneeAngle_c3d.posTrans);
                        start = max(1, min(start, length(tmp)));
                        stop = max(1, min(stop, length(tmp)));

                        if stop > start
                            kneeStepRot(j) = median(tmp(start:stop), 'omitnan');
                        end
                    end

                    % Store median of medians for this file
                    kneeRots(i) = round(median(kneeStepRot, 'omitnan'), 0);
                else
                    warning('%s or %s not found in file: %s', kneeVariable, strikeField, c_File);
                end

            catch ME
                warning('Error reading %s: %s', c_File, ME.message);
            end
        end

        % Final calculation.
        overallMedianKneeRot = round(median(kneeRots, 'omitnan'), 0);
        tt_angle_used = overallMedianKneeRot * -1; % to have it similar to TorsionTool where ext is positive but in Clevenland ext. is negative.
        % fprintf('Overall median %s knee rotation: %.2f degrees\n', upper(char(side)), tt_angle_used);

    case 'fromStatic'

        %% Get TT angle directly from static direct kinematic model output from the c3d file.

        if lower(side) == 'l'
            kneeVariable = varNameKneeAngle_c3d.L;

        elseif lower(side) == 'r'
            kneeVariable = varNameKneeAngle_c3d.R;

        else
            error("The variable 'side' is not defined correctly in <getTTangleFromDirKinemStatic>")
        end

        acq = btkReadAcquisition(static);
        angles = btkGetAngles(acq);
        kneeRot = median(angles.(kneeVariable)(:,varNameKneeAngle_c3d.posTrans));

        tt_angle_used = round(kneeRot,0) * -1; % to have it similar to TorsionTool where ext is positive but in Clevenland ext. is negative.

        % Close acq
        btkCloseAcquisition(acq);
end
%% Clear variables except output to prevet memory leak.
clearvars -except tt_angle_used
end