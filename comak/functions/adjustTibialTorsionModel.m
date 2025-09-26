function [path2adjustedModel, tibTorsion_out] = adjustTibialTorsionModel(path2scaledModel, path2bin, side, path2static, tf_torsion_fromStatic, tib_torsion_Markers)
% This file adapts a *.osim model and changes the tibial torsion based on
% the torsion seen in the static trial. It uses the tibia_tibia_prox joint
% from the lenhart model for this prupose.
%
% INPUT:
%   - path2scaledModel: full path to sclaed model to adjust
%   - path2bin: full path to the bin folder with the comak executables
%   - side: string, specifies which model should be created and scaled
%     'left' or 'right'
%   - path2static: path to static files
%   - tf_torsion_fromStatic: true or false
%   - tib_torsion_Markers: knee and ankle markers. note order is important!
%     KneeLateral, KneeMedial, AnkleLateral, AnkleMedial, e.g. {'LKNE', 'LMKNE', 'LANK', 'LMMA'};

% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
%
% Last changed:         09/2022
% -------------------------------------------------------------------------

if tf_torsion_fromStatic

    % Some housekeeping variables
    side = char(lower(side(1)));

    % Get the tibial torsion
    acq = btkReadAcquisition(path2static);
    markers = btkGetMarkers(acq);

    u = mean(markers.(tib_torsion_Markers{2}) - markers.(tib_torsion_Markers{1})); % Knee
    v = mean(markers.(tib_torsion_Markers{4}) - markers.(tib_torsion_Markers{3})); % Ankle

    cosOfAngle = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    tibialTorsionStatic = real(acosd(cosOfAngle));

    % Hardcoded settings from the comak Lenhart Model - the
    % original Model's tibial torsion
    LatKnee = [-0.0676933 0.515165 -0.141643];
    MedKnee = [-0.067021 0.515165 -0.0404022];
    LatAnkle = [-0.0911261 0.109623 -0.140404];
    MedAnkle = [-0.0601182 0.126559 -0.0450755];
    u = MedKnee - LatKnee;
    v = MedAnkle - LatAnkle;

    cosOfAngle = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    tibialTorsionOriginal = real(acosd(cosOfAngle));

    % Torsion to adjust
    tibTorsion = round(tibialTorsionStatic - tibialTorsionOriginal); 

    % Identify if int or ext correction
    if tibTorsion == 0
        corrName = '';
    elseif tibTorsion < 0
        corrName = strcat('int_', num2str(abs(tibTorsion))); % replace decimal with underline
    elseif tibTorsion > 0
        corrName = strcat('ext_', num2str(abs(tibTorsion))); % replace decimal with underline
    end

    % Prepare angle correction depending if it is a right or a left model
    % Right side, negative values = external rot
    % Left side, positive values = external rot.
    if tibTorsion < 0 && strcmp(side(1), 'r') % if neg. tibTors (I want to int. rot tibia) I need a pos. value for the right side
        tibTorsion = tibTorsion * -1;
    elseif tibTorsion > 0 && strcmp(side(1), 'r') % if pos. tibTors (I want to ext. rot tibia) I need a neg. value for the right side
        tibTorsion = tibTorsion * -1;
    end

    % Calculate radians from deg
    tibTorsion = tibTorsion * pi/180;

    % Load the API, and jam plugin
    import org.opensim.modeling.*
    opensimCommon.LoadOpenSimLibrary(strcat(path2bin,'jam_plugin'));

    % Create a copy of the scaled model, change name, and create output path
    [filepath,name,ext] = fileparts(path2scaledModel);
    path2Model2Adjust = strcat(filepath,'\',name,'-tt_',corrName, ext); % <-- This is the output of this function
    copyfile(path2scaledModel,path2Model2Adjust);

    % Change the tibio-femoral alignment in the specified model.
    % Load model
    model = Model(path2Model2Adjust);

    % Get tibia proximal body and change rot orientation. 
    tibJointSet = model.get_JointSet.get(strcat('tibia_tibia_proximal_', side));
    tibFrame = tibJointSet.get_frames(0); % there is only one frame
    tt_angleOsim = ArrayDouble.createVec3([0, tibTorsion, 0]);
    tibFrame.set_orientation(tt_angleOsim);

    % convert tibtorsion back to degress
    tibTorsion_out = tibTorsion *(180/pi);

    % Save model
    model.print(path2Model2Adjust);
    path2adjustedModel = path2Model2Adjust;
else 
    % Do nothing in case of false
    path2adjustedModel = path2scaledModel;
    tibTorsion_out = NaN;    
end

%% Clear variables except output to prevet memory leak.
clearvars -except path2adjustedModel tibTorsion_out
end