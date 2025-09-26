%-------------------------------------------------------------------------%
% Copyright (c) 2021 % Kirsten Veerkamp, Hans Kainz, Bryce A. Killen,     %
%    Hulda J�nasd�ttir, Marjolein M. van der Krogt      		          %
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
%    Authors: Hulda J�nasd�ttir & Kirsten Veerkamp                        %
%                            February 2021                                %
%    email:    k.veerkamp@amsterdamumc.nl                                 %
% ----------------------------------------------------------------------- %
% The muscle attachments for the femur are put in one matrix.
% 
% 12/26/2022
% changes Elias Wallnoefer:
% try catch blocks -> either work with architecture of thelen or millard muscles
%
%----------------------------------------------------------------------- %

function [femurMuscle, femurPlace1, femurNR, femurMuscleType] = femur_MA(dataModel, answerLeg, rightbone)

%% Find the muscles in the model
muscleType = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);
muscleType = muscleType{1};
muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);

%Identify the left and right leg
if strcmp(answerLeg, rightbone) == 1;
    femurMA = 'femur_r';
else
    femurMA = 'femur_l';
end

femurMuscle =[];
femurPlace1 = {};
femurNR = [];
femurMuscleType = {};

muscleTypes = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects);

for k = 1 : numel(muscleTypes)
    muscleType = muscleTypes{k};
    muscles = dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType);

    for i = 1:size(muscles,2)
        if isfield(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}, 'GeometryPath')
            AttachmentSize = fieldnames(dataModel.OpenSimDocument.Model.ForceSet.objects.(muscleType){1,i}.GeometryPath.PathPointSet.objects);
        else
            AttachmentSize = [];
        end
        
        for j = 1 : numel(AttachmentSize)
            MuscleAttachments = muscles{1,i}.GeometryPath.PathPointSet.objects.(AttachmentSize{j});
    
            if size(MuscleAttachments, 2) == 1 % only one PathPoint of this kind in this muscle (e.g. rect_fem PathPoint in gait2392)
                try
                    CompareStrings1_femur = strcmp(femurMA, MuscleAttachments.body.Text);
                catch ME
                    if contains(ME.message, 'body')
                        CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments.socket_parent_frame.Text,'/bodyset/',''));
                    else
                        error('Unknown Error in femur_MA check try-catch blocks')
                    end
                end
                if CompareStrings1_femur == 1;
                    try
                        femurMuscle = [femurMuscle; str2num(MuscleAttachments.location.Text)];
                    catch
                        femurMuscle = [femurMuscle; str2num(MuscleAttachments.socket_parent_frame.Text)];
                    end
                    femurMuscleType = [femurMuscleType; muscleType];
                    femurNR = [femurNR; i];
                    place1 = AttachmentSize{j};
                    femurPlace1 = [femurPlace1;place1];
                end
            else    % more of a kind - we need to add {1,%d} after the type
                for ii = 1:size(MuscleAttachments,2)
                    try
                        CompareStrings1_femur = strcmp(femurMA, MuscleAttachments{1,ii}.body.Text);
                    catch ME
                        if contains(ME.message, 'body')
                            CompareStrings1_femur = strcmp(femurMA, regexprep(MuscleAttachments{1,ii}.socket_parent_frame.Text,'/bodyset/',''));
                        else
                            error('Unknown Error in femur_MA check try-catch blocks')
                        end
                    end
                    if CompareStrings1_femur == 1;
                        try
                            femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.location.Text)];
                        catch
                            femurMuscle = [femurMuscle; str2num(MuscleAttachments{1,ii}.socket_parent_frame.Text)];
                        end
                        femurMuscleType = [femurMuscleType; muscleType];
                        femurNR = [femurNR; i];
                        place1 = [AttachmentSize{j} '{1,' num2str(ii) '}'];
                        femurPlace1 = [femurPlace1;place1];
                    end
                end
            end
        end
    % %     % The attachments are found depending on how many types of attachments there are for each muscle, the for loop is devided into 3 part.
    % %     try
    % %         AttachmentSize = size(struct2cell(dataModel.OpenSimDocument.Model.ForceSet.objects.Thelen2003Muscle{1,i}.GeometryPath.PathPointSet.objects),1);
    % %     catch ME
    % %         if contains(ME.message,'Thelen2003Muscle')
    % %             AttachmentSize = size(struct2cell(dataModel.OpenSimDocument.Model.ForceSet.objects.Millard2012EquilibriumMuscle{1,i}.GeometryPath.PathPointSet.objects),1);
    % %         else
    % %             error('Unknown Error in femur_MA regarding ForceSet-Muscles check try-catch blocks')
    % %         end
    % %     end
    % %     % if the muscle only has pathpoints
    % %     if AttachmentSize == 1
    % %         MuscleAttachments1_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.PathPoint;
    % %         for ii = 1:size(MuscleAttachments1_femur,2)
    % %             %CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %             try
    % %                 CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings1_femur = MuscleAttachments1_femur{1,ii}.socket_parent_frame.Text;
    % %                     CompareStrings1_femur = strcmp(femurMA,extractAfter(modCompareStrings1_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings1_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments1_femur{1,ii}.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('PathPoint{1,%d}',ii);
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %         end
    % %         % if the muscle has pathpoints and conditional path points
    % %     elseif AttachmentSize ==2
    % %         % The muscle attachemtns of the type pathpoint
    % %         MuscleAttachments1_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.PathPoint;
    % %         for ii = 1:size(MuscleAttachments1_femur,2)
    % %             %CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %             try
    % %                 CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings1_femur = MuscleAttachments1_femur{1,ii}.socket_parent_frame.Text;
    % %                     CompareStrings1_femur = strcmp(femurMA,extractAfter(modCompareStrings1_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings1_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments1_femur{1,ii}.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('PathPoint{1,%d}',ii);
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %         end
    % %         % The muscle attachment of the type conditional path point
    % %         MuscleAttachments2_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.ConditionalPathPoint;
    % %         if size(MuscleAttachments2_femur) == 1;
    % %             %CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur.body.Text);
    % %             try
    % %                 CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings2_femur = MuscleAttachments2_femur.socket_parent_frame.Text;
    % %                     CompareStrings2_femur = strcmp(femurMA,extractAfter(modCompareStrings2_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings2_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments2_femur.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('ConditionalPathPoint');
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %         else
    % %             for ii = 1:size(MuscleAttachments2_femur,2)
    % %                 %CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur{1,ii}.body.Text);
    % %                 try
    % %                     CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur{1,ii}.body.Text);
    % %                 catch ME
    % %                     if contains(ME.message, 'body')
    % %                         modCompareStrings2_femur = MuscleAttachments2_femur{1,ii}.socket_parent_frame.Text;
    % %                         CompareStrings2_femur = strcmp(femurMA,extractAfter(modCompareStrings2_femur,"/bodyset/"));
    % %                     else
    % %                         error('Unknown Error in femur_MA check try-catch blocks')
    % %                     end
    % %                 end
    % %                 if CompareStrings2_femur == 1;
    % %                     femurMuscle = [femurMuscle; str2num(MuscleAttachments2_femur{1,ii}.location.Text)];
    % %                     femurNR = [femurNR; i];
    % %                     place1 = sprintf('ConditionalPathPoint{1,%d}',ii);
    % %                     femurPlace1 = [femurPlace1;place1];
    % %                 end
    % %             end
    % %         end
    % %         % if the msucle has pathpoints, conditional path point and moving
    % %         % path point
    % %     elseif AttachmentSize == 3
    % %         %The muscle attachments of the type pathpoint
    % %         MuscleAttachments1_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.PathPoint;
    % %         % for the muscle that only have one pathpoint
    % %         if size(MuscleAttachments1_femur) == 1;
    % %             %CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur.body.Text);
    % %             try
    % %                 CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings1_femur = MuscleAttachments1_femur.socket_parent_frame.Text;
    % %                     CompareStrings1_femur = strcmp(femurMA,extractAfter(modCompareStrings1_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings1_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments1_femur.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('PathPoint');
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %             % for the muscle that have more than one pathpoint
    % %         else
    % %             for ii = 1:size(MuscleAttachments1_femur,2)
    % %                 %CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %                 try
    % %                     CompareStrings1_femur = strcmp(femurMA, MuscleAttachments1_femur{1,ii}.body.Text);
    % %                 catch ME
    % %                     if contains(ME.message, 'body')
    % %                         modCompareStrings1_femur = MuscleAttachments1_femur{1,ii}.socket_parent_frame.Text;
    % %                         CompareStrings1_femur = strcmp(femurMA,extractAfter(modCompareStrings1_femur,"/bodyset/"));
    % %                     else
    % %                         error('Unknown Error in femur_MA check try-catch blocks')
    % %                     end
    % %                 end
    % %                 if CompareStrings1_femur == 1;
    % %                     femurMuscle = [femurMuscle; str2num(MuscleAttachments1_femur{1,ii}.location.Text)];
    % %                     femurNR = [femurNR; i];
    % %                     place1 = sprintf('PathPoint{1,%d}',ii);
    % %                     femurPlace1 = [femurPlace1;place1];
    % %                 end
    % %             end
    % %         end
    % %         %The muscle attachments of the type conditional path point
    % %         MuscleAttachments2_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.ConditionalPathPoint;
    % %         %The muscles that only have one conditional path point
    % %         if size(MuscleAttachments2_femur) == 1
    % %             %CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur.body.Text);
    % %             try
    % %                 CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings2_femur = MuscleAttachments2_femur.socket_parent_frame.Text;
    % %                     CompareStrings2_femur = strcmp(femurMA,extractAfter(modCompareStrings2_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings2_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments2_femur.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('ConditionalPathPoint');
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %             % The muscle that have more than one conditional path point
    % %         else
    % %             for ii = 1:size(MuscleAttachments2_femur,2)
    % %                 %CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur{1,ii}.body.Text);
    % %                 try
    % %                     CompareStrings2_femur = strcmp(femurMA, MuscleAttachments2_femur{1,ii}.body.Text);
    % %                 catch ME
    % %                     if contains(ME.message, 'body')
    % %                         modCompareStrings2_femur = MuscleAttachments2_femur{1,ii}.socket_parent_frame.Text;
    % %                         CompareStrings2_femur = strcmp(femurMA,extractAfter(modCompareStrings2_femur,"/bodyset/"));
    % %                     else
    % %                         error('Unknown Error in femur_MA check try-catch blocks')
    % %                     end
    % %                 end
    % %                 if CompareStrings2_femur == 1;
    % %                     femurMuscle = [femurMuscle; str2num(MuscleAttachments2_femur{1,ii}.location.Text)];
    % %                     femurNR = [femurNR; i];
    % %                     place1 = sprintf('ConditionalPathPoint{1,%d}',ii);
    % %                     femurPlace1 = [femurPlace1;place1];
    % %                 end
    % %             end
    % %         end
    % %         % The muscle attachments of the type moveing path points
    % %         MuscleAttachments3_femur = muscles{1,i}.GeometryPath.PathPointSet.objects.MovingPathPoint;
    % %         % The muscles with one moveing pathpoint
    % %         if size(MuscleAttachments3_femur) == 1
    % %             %CompareStrings3_femur = strcmp(femurMA, MuscleAttachments3_femur.body.Text);
    % %             try
    % %                 CompareStrings3_femur = strcmp(femurMA, MuscleAttachments3_femur.body.Text);
    % %             catch ME
    % %                 if contains(ME.message, 'body')
    % %                     modCompareStrings3_femur = MuscleAttachments3_femur.socket_parent_frame.Text;
    % %                     CompareStrings3_femur = strcmp(femurMA,extractAfter(modCompareStrings3_femur,"/bodyset/"));
    % %                 else
    % %                     error('Unknown Error in femur_MA check try-catch blocks')
    % %                 end
    % %             end
    % %             if CompareStrings3_femur == 1;
    % %                 femurMuscle = [femurMuscle; str2num(MuscleAttachments3_femur.location.Text)];
    % %                 femurNR = [femurNR; i];
    % %                 place1 = sprintf('MovingPathPoint');
    % %                 femurPlace1 = [femurPlace1;place1];
    % %             end
    % %             % The muscles with more than one moving path point
    % %         else
    % %             for ii = 1:size(MuscleAttachments3_femur,2)
    % %                 %CompareStrings3_femur = strcmp(femurMA, MuscleAttachments3_femur{1,ii}.body.Text);
    % %                 try
    % %                     CompareStrings3_femur = strcmp(femurMA, MuscleAttachments3_femur.body.Text);
    % %                 catch ME
    % %                     if contains(ME.message, 'body')
    % %                         modCompareStrings3_femur = MuscleAttachments3_femur.socket_parent_frame.Text;
    % %                         CompareStrings3_femur = strcmp(femurMA,extractAfter(modCompareStrings3_femur,"/bodyset/"));
    % %                     else
    % %                         error('Unknown Error in femur_MA check try-catch blocks')
    % %                     end
    % %                 end
    % %                 if CompareStrings3_femur == 1;
    % %                     femurMuscle = [femurMuscle; str2num(MuscleAttachments3_femur{1,ii}.location.Text)];
    % %                     femurNR = [femurNR; i];
    % %                     place1 = sprintf('MovingPathPoint{1,%d}',ii);
    % %                     femurPlace1 = [femurPlace1;place1];
    % %                 end
    % %             end
    % %         end
    % %     end
    end
end
disp('The muscle attachments have been rotated')