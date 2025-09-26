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

% ----------------------------------------------------------------------- %
clear all; close all
mfile_name = mfilename('fullpath');
[pathstr,name,ext] = fileparts(mfile_name);
cd(pathstr);
try
    delete('DEFORMED_MODEL/*');
end
addpath(genpath(pwd))

%% find the block for your model, and uncomment it.
% don't forget to comment all other blocks!
%% values for Rajagopal as base model
model = 'Rajagopal/Rajagopal2015.osim';
GeometryFolder = 'Rajagopal/Geometry';
applyTibiaTorsionToJointOffset = 1;
% measured by Hans, Basilio and Willi with Sangeux 2015 (Femur, doi:10.1097/RCT.0000000000000161) and Yan 2019 (Tibia,
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6803189/)
default_Anteversion = 21; 
default_NeckShaftAngle = 121;
default_TibiaTorsion = 24;
% %Method after Hernandez - measured by Willi
% default_Anteversion = 13; 

%% values for Lernagopal (Lerner + Rajagopal) as base model
% model = 'Lernagopal/Lernagopal_41_OUF.osim';
% GeometryFolder = 'Lernagopal/Geometry';
% applyTibiaTorsionToJointOffset = 1;
% % same as Rajagopal - measured by Hans, Basilio and Willi with Sangeux
% % 2015 (Femur, doi:10.1097/RCT.0000000000000161) and Yan 2019 (Tibia,
% % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6803189/)
% default_Anteversion = 21; 
% default_NeckShaftAngle = 121;
% default_TibiaTorsion = 24;
% %Method after Hernandez - measured by Willi
% default_Anteversion = 13; 


%% values for gait2392 as base model
% model = 'gait2392_genericsimplOS4.osim';
% GeometryFolder = 'C:\OpenSim 4.2\Geometry';
% % applyTibiaTorsionToJointOffset = 0;
% applyTibiaTorsionToJointOffset = 1;
% not measured yet!!!!  you have to measure this on the OpenSim Geometry with the same method as with your participants!
% default_Anteversion = 21; 
% default_NeckShaftAngle = 121;
% default_TibiaTorsion = 24;

%% values for Lenhart as base model
% model = 'Lenhart/OActiveLeft_modified.osim';
% model = 'Lenhart/OActiveRight_modified.osim';
% GeometryFolder = 'Lenhart/Geometry';
% applyTibiaTorsionToJointOffset = 0;
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

markerset = 'MarkerSet.xml'; 

deform_bone = 'F'; 
which_leg = 'R'; 
angle_AV_right = 32; % right anteversion angle (in degrees) %
angle_NS_right = 165; % right neck-shaft angle (in degrees) %
deformed_model = ['rightNSA' num2str(angle_NS_right) '_rightAVA' num2str(angle_AV_right) ];

make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_AV_right - default_Anteversion, angle_NS_right - default_NeckShaftAngle, GeometryFolder);

%% left femur
model = [deformed_model '.osim']; 
markerset = [deformed_model '_' markerset]; 

deform_bone = 'F'; 
which_leg = 'L'; 
angle_AV_left = 38; % left anteversion angle (in degrees) %
angle_NS_left = 121; % left neck-shaft angle (in degrees) %
deformed_model = [ 'leftNSA' num2str(angle_NS_left) '_leftAVA' num2str(angle_AV_left)]; 
make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_AV_left - default_Anteversion, angle_NS_left - default_NeckShaftAngle, GeometryFolder);

%% right tibia
model = [deformed_model '.osim']; 
markerset = [deformed_model '_' markerset]; 

deformed_model = 'RT15'; 
deform_bone = 'T'; 
which_leg = 'R';
angle_TT_right = 24; % right tibial torsion angle (in degrees) %
deformed_model = [ 'rightTT' num2str(angle_TT_right) ];

make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_TT_right -default_TibiaTorsion, [], GeometryFolder, applyTibiaTorsionToJointOffset);

%% left tibia
model = [deformed_model '.osim']; 
markerset = [deformed_model '_' markerset]; 

deformed_model = 'LT5';
deform_bone = 'T';
which_leg = 'L'; 
angle_TT_left = 24; % left tibial torsion angle (in degrees) %
deformed_model = [ 'leftTT' num2str(angle_TT_left) ];

make_PEmodel( model, deformed_model, markerset, deform_bone, which_leg, angle_TT_left -default_TibiaTorsion, [], GeometryFolder, applyTibiaTorsionToJointOffset);
