%-------------------------------------------------------------------------%
% Copyright (c) 2021 % Kirsten Veerkamp, Hans Kainz, Bryce A. Killen,     %
%    Hulda Jónasdóttir, Marjolein M. van der Krogt     		              %
%                                                                         %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         %
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
%                                                                         %
%    Authors: Hulda Jónasdóttir & Kirsten Veerkamp                        %
%                            February 2021                                %
%    email:    k.veerkamp@amsterdamumc.nl                                 %
%                                                                         %
%  modified by Elias Wallnöfer & Willi Koller   (Mai 2023)                %
%                                                                         %
%  email: willi.koller@univie.ac.at                                       %
%                                                                         %
% ----------------------------------------------------------------------- %
%%%  main script to create opensim model with personalised geometries   %%%
% Give the subject-specific femoral anteversion (AV) and neck-shaft (NS) angles,
% 	as well as the tibial torsion (TT) angles, as input for the right and left leg.
% 	Lines which require these inputs are indicated by a % at the end of the line.
% The final model with personalised torsions is saved in the DEFORMED_MODEL
% 	folder, and is called FINAL_PERSONALISEDTORSIONS.osim.
% 	The adjusted markerset can also be found in this folder.
%
% note1: The angle definitions for AV and TT are as follows:
% 	- AV: positive: femoral anteversion; negative: femoral retroversion.
% 	- TT: positive: external rotation; negative: internal rotation.
% note2: Adjust the MarkerSet.xml in the main folder to your marker set,
% 	when using markers for the greater trochanter (when adjusting
% 	femur) and/or when using markers on the feet (when adjusting tibia).
% note3: If you only wish to adjust the femoral geometry (and not the tibial
% 	torsion), set the input to the tibial torsion to 0 degrees (=default
% 	tibial torsion in generic femur).
% note4: Default angles of the generic OpenSim model geometry should be
%   measured with the same method (e.g. Hernandez, ...) which you use for your
%   partipants to ensure consistency.
% note5: in models with WrapObjects it is MANDATORY to check muscle moment
%   arms for each movement before running muscle specific analysis (e.g.
%   Static Optimization).
%   It can happen that muscles pull through WrapObjects, this has to be
%   avoided. One possible solution to tackle this problem is to resize the
%   WrapObjects.
%   You can use the "checkMuscleMomentArms" script to identfiy problems.


% applyTibiaTorsionToJointOffset = 0 is the original method where torsion
% is applied via translation and rotation axis and not via body coordinate
% system rotation. This method is not applicable with Rajagopal model
% because it does not have these elements...

%
%
% 12/26/2022
%  changes Elias Wallnoefer:
%  Femur torsion should now work with RajagopalModel
%  FinalModel = "leftNSA*"
%
%  ! Tibia torsion does not yet work - hard coded line-numbers of XML need to
%  be fixed in tibia.m and tibia_locationInParent-rotation retested + adapted for both models
%
%
% 14/03/2023
%  changes by Willi Koller
%     Femur and Tibia Torsion should now work with all models, testet with gait2392, Hamner, Rajagopal, Lernagopal
%     WrapObjects in the proximal part of the femur (part which is rotated) are not supported yet, you have to adjust the location and rotation manually.
%     A message is written to the console if this is the case!
%
% 15/03/2023
%  changes by Willi Koller for Lenhart model
%     Lenhart model - need to adjust Ligament locations!
%     also rotate additional geometries of foot
%
% 16/03/2023
%  changes by Willi Koller for Lenhart model
%     should work now with ligaments and additional geometries

% Adapted by bhorsak 08/2024 to make it work within my automated SIM (autoSIM) workflow.
% ----------------------------------------------------------------------- %
function [path2adjustedModel, markerSet_rotated, torsionAdjustOut] = TorsionToolAllModels(torsiontool, persInfo, path, path2Model, rootWorkingDirectory, workingDirectory, side)

%% Initialize output values - if used they will be overwritten
% Store adaption settings

% Logicals
if strcmp(side, 'bothSides')
    torsionAdjustOut.tibTorsionAdaption = mat2str(torsiontool.tibTorsionAdaption);
    torsionAdjustOut.femurAntetorsionAdaption = mat2str(torsiontool.femurAntetorsionAdaption);
    torsionAdjustOut.neckShaftAdaption = mat2str(torsiontool.neckShaftAdaption);
else
    torsionAdjustOut.(strcat('tibTorsionAdaption_',lower(side(1)))) = mat2str(torsiontool.tibTorsionAdaption);
    torsionAdjustOut.(strcat('femurAntetorsionAdaption_',lower(side(1)))) = mat2str(torsiontool.femurAntetorsionAdaption);
    torsionAdjustOut.(strcat('neckShaftAdaption_',lower(side(1)))) = mat2str(torsiontool.neckShaftAdaption);
end

% Values
torsionAdjustOut.AVR = NaN;
torsionAdjustOut.AVL = NaN;
torsionAdjustOut.NSAR = NaN;
torsionAdjustOut.NSAL = NaN;
torsionAdjustOut.TTR = NaN;
torsionAdjustOut.TTL = NaN;

% Set default flag to false
defaultFlag_TT_right = false;
defaultFlag_TT_left = false;
defaultFlag_NS_right = false;
defaultFlag_NS_left = false;
defaultFlag_AV_right = false;
defaultFlag_AV_left = false;

%% Get rotations from setup file
if torsiontool.tibTorsionAdaption

    % Get tibiatorsion angle right side
    if isnumeric(persInfo.TTR_degree) && ~isnan(persInfo.TTR_degree); angle_TT_right = persInfo.TTR_degree; else; defaultFlag_TT_right = true; end

    % Get tibiatorsion angle left side
    if isnumeric(persInfo.TTL_degree) && ~isnan(persInfo.TTL_degree); angle_TT_left = persInfo.TTL_degree; else; defaultFlag_TT_left = true; end

else
    %Set to default values
    defaultFlag_TT_right = true;
    defaultFlag_TT_left = true;
end

if torsiontool.neckShaftAdaption
    % Get neckshaft angle right side
    if isnumeric(persInfo.NSAR_degree) && ~isnan(persInfo.NSAR_degree); angle_NS_right = persInfo.NSAR_degree; else; defaultFlag_NS_right = true; end

    % Get neckshaft angle left side
    if isnumeric(persInfo.NSAL_degree) && ~isnan(persInfo.NSAL_degree); angle_NS_left = persInfo.NSAL_degree; else; defaultFlag_NS_left = true; end

else
    %Set to default values
    defaultFlag_NS_right = true;
    defaultFlag_NS_left = true;
end

if torsiontool.femurAntetorsionAdaption
    % Get femur anteversion angle right side
    if isnumeric(persInfo.AVR_degree) && ~isnan(persInfo.AVR_degree); angle_AV_right = persInfo.AVR_degree; else; defaultFlag_AV_right = true; end

    % Get femur anteversion angle left side
    if isnumeric(persInfo.AVL_degree) && ~isnan(persInfo.AVL_degree); angle_AV_left = persInfo.AVL_degree; else; defaultFlag_AV_left = true; end

else
    %Set to default values
    defaultFlag_AV_right = true;
    defaultFlag_AV_left = true;
end

% Set Model side
if strcmp(side,'right')
    sideModel = '_RM';
elseif strcmp(side,'left')
    sideModel = '_LM';
else
    sideModel = ''; % used for models were both sides will be chnaged in the same model. COMAK/Lenhart needs side-specific models.
end

% Set paths
mfile_name = mfilename('fullpath');
[pathstr,name,ext] = fileparts(mfile_name);
cd(pathstr);
addpath(genpath(pwd))

% Delete old temp foldertry
try
    delete('DEFORMED_MODEL/*');
end

model = path2Model;
GeometryFolder = path.geometry;
% applyTibiaTorsionToJointOffset = 0 is the original method where torsion
% is applied via translation and rotation axis and not via body coordinate
% system rotation. This method is not applicable with Rajagopal model
% because it does not have these elements...

[~,nameMarkerSet,ext] = fileparts(path.markerSet);
markerSet_rotated = [path.scaling, nameMarkerSet, '_rotated', sideModel, ext];

copyfile(path.markerSet, markerSet_rotated)
markerset = markerSet_rotated;

%% Find the block for your model, and uncomment it.

if contains(path.genModel, 'lernergopal', 'IgnoreCase', true)
    %% Values for Lernagopal (Lerner + Rajagopal) as base model, TT, AV external is positive
    default_Anteversion = 10;       % checked by Fabian Unglaube(10) / Willi (10.1), based on Hernandez
    default_NeckShaftAngle = 120;   % checked by Fabian Unglaube ~120
    default_TibiaTorsion = 26;      % checked by Brian (26) and Willi Koller (26.5)
    applyTibiaTorsionToJointOffset = 1;

    % Set to default if no data available for adaption.
    if defaultFlag_TT_right; angle_TT_right = default_TibiaTorsion; end
    if defaultFlag_TT_left; angle_TT_left = default_TibiaTorsion; end
    if defaultFlag_NS_right; angle_NS_right = default_NeckShaftAngle; end
    if defaultFlag_NS_left; angle_NS_left = default_NeckShaftAngle; end
    if defaultFlag_AV_right; angle_AV_right = default_Anteversion; end
    if defaultFlag_AV_left; angle_AV_left = default_Anteversion; end

elseif contains(path.genModel, 'rajagopal', 'IgnoreCase', true)
    %% values for Rajagopal as base model   
    default_Anteversion = 10;           % checked by Fabian Unglaube(10) / Willi (10.1), based on Hernandez
    default_NeckShaftAngle = 121;       % checked by Fabian Unglaube ~120
    default_TibiaTorsion = 27;          % checked by Willi Koller (27.1)
    applyTibiaTorsionToJointOffset = 1;    

    % Set to default if no data available for adaption.
    if defaultFlag_TT_right; angle_TT_right = default_TibiaTorsion; end
    if defaultFlag_TT_left; angle_TT_left = default_TibiaTorsion; end
    if defaultFlag_NS_right; angle_NS_right = default_NeckShaftAngle; end
    if defaultFlag_NS_left; angle_NS_left = default_NeckShaftAngle; end
    if defaultFlag_AV_right; angle_AV_right = default_Anteversion; end
    if defaultFlag_AV_left; angle_AV_left = default_Anteversion; end

elseif contains(path.genModel, 'Lenhart', 'IgnoreCase', true) || contains(path.genModel, 'COMAK', 'IgnoreCase', true)
    %% values for Lenhart as base model
    applyTibiaTorsionToJointOffset = 0;
    default_Anteversion = 10;           % checked by Fabian Unglaube(10) / Willi (10.1), based on Hernandez
    default_NeckShaftAngle = 122;       % checked by Fabian Unglaube ~122
    default_TibiaTorsion = 27;          % checked by Willi Koller (27.1)

    % Set to default if no data available for adaption.
    if defaultFlag_TT_right; angle_TT_right = default_TibiaTorsion; end
    if defaultFlag_TT_left; angle_TT_left = default_TibiaTorsion; end
    if defaultFlag_NS_right; angle_NS_right = default_NeckShaftAngle; end
    if defaultFlag_NS_left; angle_NS_left = default_NeckShaftAngle; end
    if defaultFlag_AV_right; angle_AV_right = default_Anteversion; end
    if defaultFlag_AV_left; angle_AV_left = default_Anteversion; end

else
    warning('It seems the model was not yet implemented in <TorsionToolAllModels.m>!')
    pause();
end

%% More models.

%% values for gait2392 as base model
% model = 'gait2392_genericsimplOS4.osim';
% GeometryFolder = 'C:\OpenSim 4.2\Geometry';
% % applyTibiaTorsionToJointOffset = 0;
% applyTibiaTorsionToJointOffset = 1;
% not measured yet!!!!  you have to measure this on the OpenSim Geometry with the same method as with your participants!
% default_Anteversion = 21;
% default_NeckShaftAngle = 121;
% default_TibiaTorsion = 24;

%% values for Hamner as base model
% model = 'Hamner/Hamner_baseModel.osim';
% GeometryFolder = 'Hamner/Geometry';
% applyTibiaTorsionToJointOffset = 1;
% not measured yet!!!!  you have to measure this on the OpenSim Geometry with the same method as with your participants!
% default_Anteversion = 21;
% default_NeckShaftAngle = 121;
% default_TibiaTorsion = 24;
TorsToolUsedCheck = 0;
if torsiontool.neckShaftAdaption || torsiontool.femurAntetorsionAdaption
    %% right femur
    deform_bone = 'F';
    which_leg = 'R';

    if angle_AV_right - default_Anteversion ~= 0 || angle_NS_right - default_NeckShaftAngle ~= 0
        deformed_model = ['tmpTorsionTool_rightNSA' num2str(angle_NS_right) '_rightAVA' num2str(angle_AV_right) ];
        make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_AV_right - default_Anteversion, angle_NS_right - default_NeckShaftAngle, GeometryFolder, [], workingDirectory, path.genModel4Scaling, sideModel);

        % Track if TorsTool was used.
        TorsToolUsedCheck = TorsToolUsedCheck +1;

        %Write variable for setting info output
        torsionAdjustOut.AVR = angle_AV_right;
        torsionAdjustOut.NSAR = angle_NS_right;
    else
        deformed_model = model(1:end-5);
    end

    %% left femur
    model = [deformed_model '.osim'];
    deform_bone = 'F';
    which_leg = 'L';

    % Run tool only if values are different from standard
    if  angle_AV_left - default_Anteversion ~= 0 || angle_NS_left - default_NeckShaftAngle ~= 0
        deformed_model = [ 'tmpTorsionTool_leftNSA' num2str(angle_NS_left) '_leftAVA' num2str(angle_AV_left)];
        make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_AV_left - default_Anteversion, angle_NS_left - default_NeckShaftAngle, GeometryFolder, [], workingDirectory, path.genModel4Scaling, sideModel);

        % Track if TorsTool was used.
        TorsToolUsedCheck = TorsToolUsedCheck +1;

        %Write variable for setting info output
        torsionAdjustOut.AVL = angle_AV_left;
        torsionAdjustOut.NSAL = angle_NS_left;

        model = [deformed_model '.osim'];

    else
        deformed_model = model(1:end-5);
    end

end

if torsiontool.tibTorsionAdaption
    %% right tibia
    deform_bone = 'T';
    which_leg = 'R';

    % Run tool only if values are different from standard
    if angle_TT_right - default_TibiaTorsion ~= 0
        deformed_model = [ 'tmpTorsionTool_rightTT' num2str(angle_TT_right) ];
        make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_TT_right -default_TibiaTorsion, [], GeometryFolder, applyTibiaTorsionToJointOffset, workingDirectory, path.genModel4Scaling, sideModel);

        % Track if TorsTool was used.
        TorsToolUsedCheck = TorsToolUsedCheck +1;

        %Write variable for setting info output
        torsionAdjustOut.TTR = angle_TT_right;
    else
        deformed_model = model(1:end-5);
    end

    %% left tibia
    model = [deformed_model '.osim'];
    deform_bone = 'T';
    which_leg = 'L';

    % Run tool only if values are different from standard
    if angle_TT_left - default_TibiaTorsion ~= 0
        deformed_model = [ 'tmpTorsionTool_leftTT' num2str(angle_TT_left) ];
        make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_TT_left -default_TibiaTorsion, [], GeometryFolder, applyTibiaTorsionToJointOffset, workingDirectory, path.genModel4Scaling, sideModel);

        % Track if TorsTool was used.
        TorsToolUsedCheck = TorsToolUsedCheck +1;

        %Write variable for setting info output
        torsionAdjustOut.TTL = angle_TT_left;
    else
        deformed_model = model(1:end-5);
    end
end

% Model rotation done. Close all figures
close all

%% Add path to adjusted model as output

[pathstr,name,ext] = fileparts(path2Model);
if ~contains(name, 'tors') && TorsToolUsedCheck > 0
    path2adjustedModel = strcat(pathstr,'\',name, '-tors', ext); % <-- This is the output of this function
else
    % In case the model ran through the torision tool already and to make
    % sure I do not get something like "-tors-tors"
    path2adjustedModel = strcat(pathstr,'\',name, ext); % <-- This is the output of this function

end
% In some cases deformed model will be either the full path OR just a file
% name. This makes sure it will work for both ways.
[~,name,~] = fileparts(deformed_model);

% Only copy if target and adjusted model are not the same.
target = [workingDirectory, name,'.osim'];
if ~isequal(path2adjustedModel, target)
    movefile(target, path2adjustedModel)
end

%% Clear vars to prevent memory leak
clearvars -except path2adjustedModel markerSet_rotated torsionAdjustOut
end