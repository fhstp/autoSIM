function [successfull, finalModel] = adaptWrappingObjectsToCorrectMuscleMomentArms(origModelFilename, motionFilename, path2JamPlugIn)
% This function checks the muscle moment arms of a scaled model for a
% movement task given by the motionFile. In case it detects discontinuities
% it iteratively decreses the radius of the wrapping object until the
% discontinuities are resolved.

% Written by: Willi Koller
% Downloaded and modified by B. Horsak (24.09.2024) from
% https://github.com/WilliKoller/OpenSimMatlabBasic
% -------------------------------------------------------------------------

%% Settings.
stepSize = 0.001;  % 1 mm
threshold = 0.004; % 4 mm

% Decide if plots should be created.
createPlot = false;

%% Prepare models

% Load opensim API and comak-plugin.
import org.opensim.modeling.*
if ~isempty(path2JamPlugIn)
    opensimCommon.LoadOpenSimLibrary(path2JamPlugIn);
end

% Load motion file.
motion = Storage(motionFilename);

% Load model and create working copy.
tmpFilename = strrep(origModelFilename, '.osim', '_modWO.osim');
copyfile(origModelFilename, tmpFilename);

% Set diary log file.
logFile = fullfile([tmpFilename(1:end-11) '_modifiedWrapObjects.log']);

% Delete old log files, in case there are any.
if exist(logFile, 'file')
    delete(logFile);
end

% Start logging.
diary(fullfile([tmpFilename(1:end-5) '_modifiedWrapObjects.log']));

%% Start while loop to test for discontinuities.
hasDiscontinuities = 1;
failed = 0; % solving discontinuities did not work
discConFlag = 0; % discontinuities detected and solved

% Get initial model.
model = Model(tmpFilename);
model.initSystem();
state = model.initSystem();

% Initialize log files to track changes.
wrapObjectsModified = {};
wrapObjectsOrigRadius = [];
wrapObjectsModifiedRadius = [];
wrapEllipsoidsModified = {};
wrapEllipsoidsOrigDimension = [];
wrapEllipsoidsModifiedDimension = [];

% Safety measure to prevent endless while loop
timeout = 60*30;  % default = 30 minutes; Max allowed duration (in seconds)
tic  % Start timing

while(hasDiscontinuities)

	% Check elapsed time
    if toc > timeout
        warning('While loop in check-moment-arms terminated after exceeding time limit of %.0f minutes.', timeout/60);
        failed = 1;
        break;
    end

    % Check which model is used
    if contains(origModelFilename, {'lernergopal', 'rajagopal'}, 'IgnoreCase', true)
    % Get Coordinate indices.
    hip_flexIndexL = model.getCoordinateSet().getIndex('hip_flexion_l');
    hip_flexCoordL = model.updCoordinateSet().get(hip_flexIndexL);

    hip_rotIndexL = model.getCoordinateSet().getIndex('hip_rotation_l');
    hip_rotCoordL = model.updCoordinateSet().get(hip_rotIndexL);

    hip_addIndexL = model.getCoordinateSet().getIndex('hip_adduction_l');
    hip_addCoordL = model.updCoordinateSet().get(hip_addIndexL);

    hip_flexIndexR = model.getCoordinateSet().getIndex('hip_flexion_r');
    hip_flexCoordR = model.updCoordinateSet().get(hip_flexIndexR);

    hip_rotIndexR = model.getCoordinateSet().getIndex('hip_rotation_r');
    hip_rotCoordR = model.updCoordinateSet().get(hip_rotIndexR);

    hip_addIndexR = model.getCoordinateSet().getIndex('hip_adduction_r');
    hip_addCoordR = model.updCoordinateSet().get(hip_addIndexR);

    knee_flexIndexL = model.getCoordinateSet().getIndex('knee_angle_l');
    knee_flexCoordL = model.updCoordinateSet().get(knee_flexIndexL);

    knee_flexIndexR = model.getCoordinateSet().getIndex('knee_angle_r');
    knee_flexCoordR = model.updCoordinateSet().get(knee_flexIndexR);
    
    elseif contains(origModelFilename, 'COMAK', 'IgnoreCase', true)
    % Get Coordinate indices.
    hip_flexIndexL = model.getCoordinateSet().getIndex('hip_flex_l');
    hip_flexCoordL = model.updCoordinateSet().get(hip_flexIndexL);

    hip_rotIndexL = model.getCoordinateSet().getIndex('hip_rot_l');
    hip_rotCoordL = model.updCoordinateSet().get(hip_rotIndexL);

    hip_addIndexL = model.getCoordinateSet().getIndex('hip_add_l');
    hip_addCoordL = model.updCoordinateSet().get(hip_addIndexL);

    hip_flexIndexR = model.getCoordinateSet().getIndex('hip_flex_r');
    hip_flexCoordR = model.updCoordinateSet().get(hip_flexIndexR);

    hip_rotIndexR = model.getCoordinateSet().getIndex('hip_rot_r');
    hip_rotCoordR = model.updCoordinateSet().get(hip_rotIndexR);

    hip_addIndexR = model.getCoordinateSet().getIndex('hip_add_r');
    hip_addCoordR = model.updCoordinateSet().get(hip_addIndexR);

    knee_flexIndexL = model.getCoordinateSet().getIndex('knee_flex_l');
    knee_flexCoordL = model.updCoordinateSet().get(knee_flexIndexL);

    knee_flexIndexR = model.getCoordinateSet().getIndex('knee_flex_r');
    knee_flexCoordR = model.updCoordinateSet().get(knee_flexIndexR);
    else
        Error('The used *.osim model is not implemented yet in <adaptWrappingObjectsTo CorrectMuscleMomentArms.m> function!')
    end


    % Get main muscles.
    numMuscles = model.getMuscles().getSize();
    muscleIndices = []; muscleNames = {};
    for i = 0 : numMuscles - 1
        tmp_muscleName = char(model.getMuscles().get(i).getName());
        if contains(tmp_muscleName, 'add') || contains(tmp_muscleName, 'gl') ...
                || contains(tmp_muscleName, 'semi') || contains(tmp_muscleName, 'bf') ...
                || contains(tmp_muscleName, 'grac') || contains(tmp_muscleName, 'piri') ...
                || contains(tmp_muscleName, 'sart') || contains(tmp_muscleName, 'tfl') ...
                || contains(tmp_muscleName, 'iliacus') || contains(tmp_muscleName, 'psoas') ...
                || contains(tmp_muscleName, 'rec') || contains(tmp_muscleName, 'gas') ...
                || contains(tmp_muscleName, 'vas') || contains(tmp_muscleName, 'sol') ...
                || contains(tmp_muscleName, 'pec')
            muscleIndices = [muscleIndices, i];
            muscleNames{end+1} = tmp_muscleName;
        end
    end

    hip_flexMomentArms = zeros(motion.getSize(), length(muscleIndices));
    hip_addMomentArms = zeros(motion.getSize(), length(muscleIndices));
    hip_rotMomentArms = zeros(motion.getSize(), length(muscleIndices));
    knee_flexMomentArms = zeros(motion.getSize(), length(muscleIndices));

    for i = 1 : motion.getSize()
        hip_flexAngleL = motion.getStateVector(i-1).getData().get(hip_flexIndexL) / 180 * pi;
        hip_rotAngleL = motion.getStateVector(i-1).getData().get(hip_rotIndexL) / 180 * pi;
        hip_addAngleL = motion.getStateVector(i-1).getData().get(hip_addIndexL) / 180 * pi;
        knee_flexAngleL = motion.getStateVector(i-1).getData().get(knee_flexIndexL) / 180 * pi;

        hip_flexAngleR = motion.getStateVector(i-1).getData().get(hip_flexIndexR) / 180 * pi;
        hip_rotAngleR = motion.getStateVector(i-1).getData().get(hip_rotIndexR) / 180 * pi;
        hip_addAngleR = motion.getStateVector(i-1).getData().get(hip_addIndexR) / 180 * pi;
        knee_flexAngleR = motion.getStateVector(i-1).getData().get(knee_flexIndexR) / 180 * pi;

        % Update the state with the hip joint angle.
        coordSet = model.updCoordinateSet();

        coordSet.get(hip_flexIndexL).setValue(state, hip_flexAngleL);
        coordSet.get(hip_rotIndexL).setValue(state, hip_rotAngleL);
        coordSet.get(hip_addIndexL).setValue(state, hip_addAngleL);
        coordSet.get(knee_flexIndexL).setValue(state, knee_flexAngleL);

        coordSet.get(hip_flexIndexR).setValue(state, hip_flexAngleR);
        coordSet.get(hip_rotIndexR).setValue(state, hip_rotAngleR);
        coordSet.get(hip_addIndexR).setValue(state, hip_addAngleR);
        coordSet.get(knee_flexIndexR).setValue(state, knee_flexAngleR);

        % Realize the state to compute dependent quantities.
        model.computeStateVariableDerivatives(state);
        model.realizeVelocity(state);

        % Compute the moment arm for each muscle.
        for j = 1:length(muscleIndices)
            muscleIndex = muscleIndices(j);

            % For left side.
            if strcmp(muscleNames{j}(end), 'l')
                hip_flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_flexCoordL);
                hip_flexMomentArms(i,j) = hip_flexMomentArm;

                hip_rotMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_rotCoordL);
                hip_rotMomentArms(i,j) = hip_rotMomentArm;

                hip_addMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_addCoordL);
                hip_addMomentArms(i,j) = hip_addMomentArm;

                knee_flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, knee_flexCoordL);
                knee_flexMomentArms(i,j) = knee_flexMomentArm;

            else % For right side.
                hip_flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_flexCoordR);
                hip_flexMomentArms(i,j) = hip_flexMomentArm;

                hip_rotMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_rotCoordR);
                hip_rotMomentArms(i,j) = hip_rotMomentArm;

                hip_addMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, hip_addCoordR);
                hip_addMomentArms(i,j) = hip_addMomentArm;

                knee_flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, knee_flexCoordR);
                knee_flexMomentArms(i,j) = knee_flexMomentArm;
            end
        end
    end

    % Now check for discontinuities.
    musclesWithDiscontinuities = {};

    for i = 1 : length(muscleIndices)

        % Hip flexion.
        dy = diff(hip_flexMomentArms(:, i));
        discontinuity_indices = find(abs(dy) > threshold);
        if size(discontinuity_indices, 1) > 0
            musclesWithDiscontinuities{end + 1} = muscleNames{i};
        end

        % Hip adduction.
        dy = diff(hip_addMomentArms(:, i));
        discontinuity_indices = find(abs(dy) > threshold);
        if size(discontinuity_indices, 1) > 0
            musclesWithDiscontinuities{end + 1} = muscleNames{i};
        end

        % Hip rotation.
        dy = diff(hip_rotMomentArms(:, i));
        discontinuity_indices = find(abs(dy) > threshold);
        if size(discontinuity_indices, 1) > 0
            musclesWithDiscontinuities{end + 1} = muscleNames{i};
        end

        % Knee flexion.
        dy = diff(knee_flexMomentArms(:, i));
        discontinuity_indices = find(abs(dy) > threshold);
        if size(discontinuity_indices, 1) > 0
            musclesWithDiscontinuities{end + 1} = muscleNames{i};
        end
    end

    % Now modify wrapping objects in case of discontinuities.
    if size(musclesWithDiscontinuities, 1) > 0

        % Set flag if discontinuities were detected at the beginning.
        discConFlag = 1;

        % Modify corresponding wrapping object
        musclesWithDiscontinuities = unique(musclesWithDiscontinuities);

        % Tell user that there are discontinuities.
        disp(' ');
        disp(['>>>> ' num2str(numel(musclesWithDiscontinuities)) ' muscle(s) with discontinuities detected. Adapting wrapping object(s) ...']);

        % Iterate through discontinuities and adapt.
        for m = 1 : numel(musclesWithDiscontinuities)
            muscleName = musclesWithDiscontinuities{m};
            muscle = model.getMuscles().get(muscleName);

            geoPath = muscle.get_GeometryPath;
            wrapSet = geoPath.getWrapSet;
            for w = 1 : wrapSet.getSize()
                wrapObj = wrapSet.get(w-1);
                wrapObjName = wrapObj.getWrapObjectName;

                for i = 1 : model.getBodySet.getSize
                    currBody = model.getBodySet.get(i-1);

                    for j = 1 : currBody.getWrapObjectSet.getSize

                        if strcmp(currBody.getWrapObjectSet.get(j-1).getName, wrapObjName)
                            wrapObject = currBody.getWrapObjectSet.get(j-1);

                            % The Lenhart model does not have wrapping zylinders, it has ellipsoids. This will make sure that in case of
                            % ellipsoides an alternative approach willbe used.
                            try
                                wrapCylinder = WrapCylinder.safeDownCast(wrapObject);
                                radius = wrapCylinder.get_radius;

                                if radius - stepSize > 0
                                    wrapCylinder.set_radius(radius - stepSize)
                                    model.print(tmpFilename);
                                    pause(0.5);                                

                                    % Display changes made.
                                    disp([muscleName ': radius decreased to ' num2str(radius - stepSize)]);
                                    
                                    % Track changes in log.
                                    wrapObjectsModified{end+1} = char(wrapObjName);
                                    wrapObjectsOrigRadius(end+1) = radius;
                                    wrapObjectsModifiedRadius(end+1) = radius - stepSize;

                                    % Update model to pick up changes in the model for the next loop.
                                    model = Model(tmpFilename);
                                    model.initSystem();
                                    state = model.initSystem();
                                    
                                    break;
                                else
                                    disp([muscleName ' moment arm wrong but wrap object is already too small! Check manually']);
                                    failed = 1;
                                    hasDiscontinuities = 0;
                                end
                            catch
                                % This makes sure that in case of OpenSim 4.0 WO are still changed since for
                                % ellipsoids there is no API in Open Sim 4.0
                                wrapEllipsoid = WrapEllipsoid.safeDownCast(wrapObject);
                                dimensionsString = wrapEllipsoid.getDimensionsString;
                                dimensions = strsplit(char(dimensionsString), ' ');
                                newDimension(1) = str2double(dimensions(2))-stepSize;
                                newDimension(2) = str2double(dimensions(3))-stepSize;
                                newDimension(3) = str2double(dimensions(4))-stepSize;
                                if sum(newDimension > 0) == 3
                                    newDimensionString = [num2str(newDimension(1)) ...
                                        ' ' num2str(newDimension(2)) ' ' num2str(newDimension(3))];

                                    % Replace dimensions string - no API exists to do this in OpenSim 4.0
                                    fid  = fopen(tmpFilename,'r');
                                    f=fread(fid,'*char')';
                                    fclose(fid);
                                    indexStartOfWO = strfind(f, ['name="' char(wrapObjName)]);
                                    indexDimensionString_0 = strfind(f(indexStartOfWO:end), '<dimensions>');
                                    indexDimensionString_1 = strfind(f(indexStartOfWO:end), '</dimensions>');
                                    indexDimensionString_0 = indexStartOfWO + indexDimensionString_0(1);
                                    indexDimensionString_1 = indexStartOfWO + indexDimensionString_1(1);
                                    newF = [f(1 : indexDimensionString_0 + 10) newDimensionString f(indexDimensionString_1 - 2 : end)];
                                    fid  = fopen(tmpFilename,'w');
                                    fprintf(fid,'%s',newF);
                                    fclose(fid);
                                    pause(0.5);

                                    % Display changes made.
                                    disp([muscleName, ': radius decreased to: ', num2str(newDimension)]);

                                    % Track changes in log.
                                    wrapEllipsoidsModified{end+1} = char(wrapObjName);
                                    wrapEllipsoidsOrigDimension{end+1} = strjoin(dimensions(2:end), ' ');
                                    wrapEllipsoidsModifiedDimension{end+1} = newDimensionString;
                                    
                                    % Update model to pick up changes in the model for the next loop.
                                    model = Model(tmpFilename);
                                    model.initSystem();
                                    state = model.initSystem();

                                    break;
                                else
                                    disp([muscleName ' moment arm wrong but wrap object is already too small! Check manually']);
                                    failed = 1;
                                    hasDiscontinuities = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        % No discontinuities detected.
        disp('No muscles with discontinuities were detected.');
        hasDiscontinuities = 0;
    end
end


%% Set final output
if failed % Discontinuities detected but not solved.
    fprintf(2, 'Something went wrong!');
    successfull = 'Discontinuities detected but not solved - original model used.';
    delete(tmpFilename) % delete temp model and keep original since corrections did not work.
    finalModel = origModelFilename; % the original model will be the output.
else
    if discConFlag % Discontinuities detected and solved
        successfull = 'Discontinuities detected and solved.';
        finalModel = tmpFilename; % the changed model will be the output but we also keep the original one.

        % Report change summary.
        uniqueWrapObjects = unique(wrapObjectsModified);
        uniqueWrapEllipsoids = unique(wrapEllipsoidsModified);
        disp(' ');
        disp('SUMMARY: ');
        if numel(uniqueWrapObjects) > 0
            for m = 1 : numel(uniqueWrapObjects)
                indizes = find(contains(wrapObjectsModified, uniqueWrapObjects{m}));
                disp([uniqueWrapObjects{m} ' radius reduced from ' num2str(wrapObjectsOrigRadius(indizes(1))) ' to ' num2str(wrapObjectsModifiedRadius(indizes(end)))])
            end
        end
        if numel(uniqueWrapEllipsoids) > 0
            for m = 1 : numel(uniqueWrapEllipsoids)
                indizes = find(contains(wrapEllipsoidsModified, uniqueWrapEllipsoids{m}));
                disp([uniqueWrapEllipsoids{m} ' dimensions reduced from ' char(wrapEllipsoidsOrigDimension{indizes(1)}) ' to ' char(wrapEllipsoidsModifiedDimension{indizes(end)}) ])
            end
        end

    else % No disconinuities detected
        successfull = 'No disconinuities detected.';
        finalModel = origModelFilename; % the original model will be the output.
        delete(tmpFilename) % delete temp model and keep original since no corrections were made.
    end
end

diary off;

%% Plotting
if createPlot
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    plot(hip_flexMomentArms);
    title(['All muscle moment arms in motion ' motionFilename]);
    legend(muscleNames, 'Interpreter', 'none');
    ylabel('Hip Flexion Moment Arm (m)');
    xlabel('Frame (after start time)');
    drawnow;
end

%% Clear variables except output to prevet memory leak.
clearvars -except successfull finalModel
end

