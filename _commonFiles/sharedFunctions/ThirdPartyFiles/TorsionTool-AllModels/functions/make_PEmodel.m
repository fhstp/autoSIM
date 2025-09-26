%-------------------------------------------------------------------------%
% Copyright (c) 2021 % Kirsten Veerkamp, Hans Kainz, Bryce A. Killen,     %
%    Hulda Jónasdóttir, Marjolein M. van der Krogt      		          %
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
% ----------------------------------------------------------------------- %
% make model with personalised torsions
%   inputs:
% answerModel = %i.e. gait2392_simbody.osim
% deformed_model = % string
% answerMarkerset =  % string
% which_leg = % string (L or R)
% angle = % in degrees
function [ ready ] = make_PEmodel( answerModel, deformed_model, answerMarkerSet, deform_bone, which_leg, angle, angle_NS, geometryFolder, applyTibiaTorsionToJointOffset, workingDirectory, genModel4Scaling, sideModel)


% % the path to save the deformed model
% place = [cd '\DEFORMED_MODEL\'];
% 
% % what model you want to deform
% if strcmp(which_leg, 'R') == 1 &&  strcmp(deform_bone, 'F') == 1;
%         answerModel_tmp = [ answerModel];
%         answerMarkerSet_tmp = [ answerMarkerSet];
% else
%     answerModel_tmp = [place answerModel];
%     answerMarkerSet_tmp = [place answerMarkerSet];
% end
% dataModel = xml2struct(answerModel_tmp);
% 
% % 
% if ~exist("DEFORMED_MODEL\Geometry", 'dir')
%     mkdir("DEFORMED_MODEL\Geometry")
% end
% The path to save the deformed model
place = workingDirectory;

% Check if paths are full paths or only file names
if exist(fileparts(answerModel)) == 0
   answerModel_tmp = fullfile(place, answerModel);
else
    answerModel_tmp = answerModel;
end

if exist(fileparts(answerMarkerSet)) == 0
   answerMarkerSet_tmp = fullfile(place, answerMarkerSet);
else
    answerMarkerSet_tmp = answerMarkerSet;
end

dataModel = xml2struct(answerModel_tmp);

% Write markerSet in model - Added by B.Guggenberger
dataMarker = xml2struct(answerMarkerSet_tmp);
dataModel.OpenSimDocument.Model.MarkerSet = dataMarker.OpenSimDocument.MarkerSet;


if strcmp(deform_bone, 'F') && strcmp(which_leg, 'R')
    for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry, 'Mesh')
                if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                        .Mesh, 2) > 1
                    for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                            .Mesh, 2)
                        try
                            vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                                .Mesh{j}.mesh_file.Text;
                            copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                        end
                    end
                else
                    try
                        vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                            .Mesh.mesh_file.Text;
                        copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                    end
                end
            end
        else
            try
                if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                        .GeometrySet.objects.DisplayGeometry, 2) > 1
                    for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                        .GeometrySet.objects.DisplayGeometry, 2)
                        try
                            vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                                .GeometrySet.objects.DisplayGeometry{j}.geometry_file.Text;
                            copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                        end
                    end

                else
                    try
                        vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                            .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
                        copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                    end
                end
            end
        end

        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'components')
            if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components, 'PhysicalOffsetFrame')
                if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame, 'attached_geometry')
                    if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry, 'Mesh')
                        if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                .Mesh, 2) > 1
                            for j = 1 : size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                    .Mesh, 2)
                                try
                                    vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                        .Mesh{j}.mesh_file.Text;
                                    copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                                end
                            end
                        else
                            try
                                vtp_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.components.PhysicalOffsetFrame.attached_geometry...
                                    .Mesh.mesh_file.Text;
                                copyfile(fullfile(geometryFolder, vtp_filename), fullfile([place 'Geometry'], vtp_filename));
                            end
                        end
                    end
                end
            end
        end
    end
end
% what you want to name the deformed model
answerNameModel = deformed_model;

% the marker set for this model.
markerset = xml2struct(answerMarkerSet_tmp);
answerLegFemur = which_leg;
answerDegFemur = angle;

for i = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body)
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['femur_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh.mesh_file.Text;
        else
            femur_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
        end
    end
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['tibia_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            tibia_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh{1,1}.mesh_file.Text;
        else
            if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                    .GeometrySet.objects.DisplayGeometry, 2) > 1
                for j = 1 : numel(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                        .GeometrySet.objects.DisplayGeometry, 2)
                    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                            .GeometrySet.objects.DisplayGeometry{j}.geometry_file.Text, 'tibia')
                        tibia_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                            .GeometrySet.objects.DisplayGeometry{j}.geometry_file.Text;
                        break;
                    end
                end
            else
                tibia_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                    .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
            end
        end
    end
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['calcn_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh, 2) > 1
                calcn_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh{1, 1}.mesh_file.Text;
            else
            calcn_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh.mesh_file.Text;
            end
        else
            calcn_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
        end
    end
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['toes_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            if size(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh, 2) > 1
                toes_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                    .Mesh{1,1}.mesh_file.Text;
            else
                toes_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                    .Mesh.mesh_file.Text;
            end
        else
            toes_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
        end
    end
    if contains(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.Attributes.name, ['talus_' lower(which_leg)])
        if isfield(dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}, 'attached_geometry')
            talus_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.attached_geometry...
                .Mesh.mesh_file.Text;
        else
            talus_filename = dataModel.OpenSimDocument.Model.BodySet.objects.Body{1,i}.VisibleObject...
                .GeometrySet.objects.DisplayGeometry.geometry_file.Text;
        end
    end
end

bone = 'T';
if strcmp(deform_bone, bone) == 1; % Rotation of the tibia
    % Ask the user if they want to rotate the left or right leg.
    
    answerLegTibia = which_leg;
    % Ask the user how large the torsion angle is.
    
    answerDegTibia = angle;
    rightboneTibia = 'R';
    if strcmp(answerLegTibia, rightboneTibia) == 1;
        %Right torsion angle is defined for the rotation
        TT_angle = -(answerDegTibia*(pi/180));
        % the geometry of the right calcn, talus and toes is imported

        dataTibia = xml2struct(fullfile(geometryFolder, tibia_filename));
        dataCalcn = xml2struct(fullfile(geometryFolder, calcn_filename));
        dataTalus = xml2struct(fullfile(geometryFolder, talus_filename));
        dataToes = xml2struct(fullfile(geometryFolder, toes_filename));
    else
        %The left torsion angle is defined for the rotation
        TT_angle = answerDegTibia*(pi/180);
        % the geometry of the left calcn, talus and toes are imported
        
        dataTibia = xml2struct(fullfile(geometryFolder, tibia_filename));
        dataCalcn = xml2struct(fullfile(geometryFolder, calcn_filename));
        dataTalus = xml2struct(fullfile(geometryFolder, talus_filename));
        dataToes = xml2struct(fullfile(geometryFolder, toes_filename));
    end
    % the script for the rotation of the tibia is called
    placeNameModel = tibia_Rajag(dataModel, markerset, answerLegTibia, rightboneTibia, TT_angle,...
        answerNameModel, answerMarkerSet, dataTibia, dataCalcn, dataTalus,...
        dataToes, place, applyTibiaTorsionToJointOffset, geometryFolder, sideModel);
    %rotation of the femur
else
    
    % femoral anteversion
    if strcmp(answerLegFemur, 'R') == 1; % Rotation of the right foot
%         FA_preAngle = 17.6;
%         NS_preAngle = 123.3;
        % The added anteversion angle is definded
        angleCorrection = answerDegFemur;% - FA_preAngle;
        FA_angle = -(angleCorrection*(pi/180));
        NS_angle = -((angle_NS)*(pi/180));
        % The geomerty of the right femur is imported
        dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
    else % Rotation of the left foot
%         FA_preAngle = 17.6;
%         NS_preAngle = 123.3;
        % The added anteversion angle is definded
        angleCorrection = answerDegFemur;% - FA_preAngle;
        FA_angle = angleCorrection*(pi/180);
        NS_angle = ((angle_NS)*(pi/180));
        % The geometry of the left femur is imported
        dataFemur = xml2struct(fullfile(geometryFolder, femur_filename));
    end
    % the script for the rotation of the femur is called.
    femur_ns(dataModel, markerset, answerLegFemur, 'R', FA_angle, NS_angle,...
        answerNameModel,answerMarkerSet, dataFemur, place, sideModel, geometryFolder);
end

% correct order of path points before writing osim file
%%
if strcmp(deform_bone, 'T')
    [~, answerNameModel] = fileparts(placeNameModel);
end
%filetext = fileread(['DEFORMED_MODEL/' answerNameModel '.osim']);
filetext = fileread([workingDirectory answerNameModel '.osim']);

GeometryPathStarts = strfind(filetext, '<GeometryPath');
GeometryPathEnds = strfind(filetext, '</GeometryPath');

newFileText = filetext;

for i = 1 : numel(GeometryPathStarts)
    section = filetext(GeometryPathStarts(i) : GeometryPathEnds(i));
    i1 = strfind(section, '<objects');
    i2 = strfind(section, '</objects');
    if ~isempty(i2)

        % Ignore "PathWrapSets" in case there is one.
        i1 = i1(1); % adjusted by bhorsak (08/2024): this makes sure that always the correct (first) PathPoint is used.
        i2 = i2(1); % adjusted by bhorsak (08/2024): this makes sure that always the correct (first) PathPoint is used.

        section = section(i1:i2);
        nameIdx = strfind(section, 'name=');
        absStart = 0;
        absEnd = 0;
        names = [];
        for j = 1 : numel(nameIdx)
            sep = section(nameIdx(j) + 5);
            sep2 = strfind(section(nameIdx(j) + 6 : end), sep);
            sep2 = sep2(1);
            names{j} = section(nameIdx(j) + 6: nameIdx(j) + 4 + sep2);
            tmp1 = strfind(section(nameIdx(j) - 25 : end), '<');
            tmp1 = nameIdx(j) - 25 + tmp1(1);
            objType = section(tmp1 : nameIdx(j)-2);
            endIdx = strfind(section(tmp1 : end), ['</' objType]);
            endIdx = tmp1 + endIdx(1); % adjusted by bhorsak (08/2024): this makes sure that always the correct (first) PathPoint is used. 
            section(tmp1-1 : endIdx + length(objType) +1);
            fullText{j} = section(tmp1-1 : endIdx + length(objType) +1);
            if j == 1
                absStart = tmp1-1;
            end
            if j == numel(nameIdx)
                absEnd = endIdx + length(objType) +1;
            end
        end
        [~, order] = sort(names);
        sectionNewOrder = '';
        for j = 1 : numel(order)
            sectionNewOrder = [sectionNewOrder fullText{order(j)}];
        end
        newFileText = strrep(newFileText, section(absStart : absEnd), sectionNewOrder);
    end
end

% fid = fopen(['DEFORMED_MODEL/' answerNameModel '.osim'],'w'); % adjusted by bhorsak (08/2024): to make work with the autoSIM workflow
fid = fopen([workingDirectory answerNameModel '.osim'],'w');
fprintf(fid, newFileText);
fclose(fid);