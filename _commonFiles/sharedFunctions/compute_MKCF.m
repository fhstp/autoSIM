function [MCFy_r, MCFy_l, LCFy_r, LCFy_l] = compute_MKCF(ScaleFactorsFile, JointReactionLoadsFile)
% COMPUTEJOINTCONTACTFORCES Computes medial and lateral knee joint contact forces
%   This function calculates medial contact force (MCF) for the right and 
%   left knee, based on OpenSim joint reaction load results and applied 
%   femur scaling factors. The implementation is adapted from OpenCap
%   Processing (see: https://github.com/stanfordnmbl/opencap-processing,
%   OpenSimPipeline/JointReaction/computeJointLoading.py, starting at line 638).
%
%   INPUTS:
%       ScaleFactorsFile - full path to the OpenSim scale factors XML file
%                          (e.g., '...\RajagopalLaiUhlrich2023-ScaleFactors.xml')
%       JointReactionLoadsFile - full path to the JointReaction loads .sto file
%                          (e.g., '...\standard-tf-2_WalkA02_1_r_JointReaction_ReactionLoads.sto')
%
%   OUTPUTS:
%       MCF_r  - medial contact force of the right knee
%       MCF_l  - medial contact force of the left knee
%
%   METHOD:
%       1. Extract femur scaling factors from the scale factors XML file.
%       2. Compute subject-specific intercondylar distances by multiplying 
%          a base intercondylar distance (0.04 m, Lerner 2015) with the
%          medial-lateral scaling factor (z-direction) of each femur.
%       3. Load joint reaction results from the .sto file.
%       4. Compute medial contact force for each knee:
%              MCF = -Fy/2 ± KAM / intercondylarDistance
%
%   REFERENCES:
%       OpenCap Processing (Github), Stanford NMBL, 2023.:
%       https://github.com/stanfordnmbl/opencap-processing/blob/main/OpenSimPipeline/JointReaction/computeJointLoading.py
%       (Line: 638, 09/2025)
%
%
%   Author: Willi Koller, adapted by B. Horsak
%   Date:   09/2025%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read scale factors XML
fid = fopen(ScaleFactorsFile, 'r');
lines = {};
tline = fgetl(fid);
while ischar(tline)
    lines{end+1} = tline; %#ok<AGROW>
    tline = fgetl(fid);
end
fclose(fid);

scale_femur_r = [];
scale_femur_l = [];

for i = 2:length(lines)
    if contains(lines{i}, '<segment>femur_r</segment>')
        scaleLine = strtrim(lines{i-1});
        if contains(scaleLine, '<scales>')
            nums = regexp(scaleLine, '<scales>\s*([\d\.\-eE\s]+)\s*</scales>', 'tokens');
            if ~isempty(nums)
                scale_femur_r = str2double(strsplit(strtrim(nums{1}{1})));
            end
        end
    end
    if contains(lines{i}, '<segment>femur_l</segment>')
        scaleLine = strtrim(lines{i-1});
        if contains(scaleLine, '<scales>')
            nums = regexp(scaleLine, '<scales>\s*([\d\.\-eE\s]+)\s*</scales>', 'tokens');
            if ~isempty(nums)
                scale_femur_l = str2double(strsplit(strtrim(nums{1}{1})));
            end
        end
    end
end

if isempty(scale_femur_r) || isempty(scale_femur_l)
    error('Could not find femur scale factors in the provided file.');
end

%% Compute subject-specific intercondylar distances
base_intercondylarDistance = 0.04; % [m], from Lerner 2015
subjectspecific_intercondylarDistance_l = base_intercondylarDistance * scale_femur_l(3);
subjectspecific_intercondylarDistance_r = base_intercondylarDistance * scale_femur_r(3);

%% Load joint reaction results
tmp_data = load_sto_file(JointReactionLoadsFile); % <-- Assumes you have a function load_sto_file.m

KAM_r = tmp_data.walker_knee_r_on_tibia_r_in_tibia_r_mx;
KAM_l = tmp_data.walker_knee_l_on_tibia_l_in_tibia_l_mx;
Fy_r  = tmp_data.walker_knee_r_on_tibia_r_in_tibia_r_fy;
Fy_l  = tmp_data.walker_knee_l_on_tibia_l_in_tibia_l_fy;

%% Compute medial contact forces
MCFy_r = -Fy_r/2 - KAM_r/subjectspecific_intercondylarDistance_r;
LCFy_r = -Fy_r - MCFy_r;
MCFy_l = -Fy_l/2 + KAM_l/subjectspecific_intercondylarDistance_l;
LCFy_l = -Fy_l - MCFy_l;



end




















% %% code to compute medial and lateral joint contact force - similar to OpenCap Processing
% % https://github.com/stanfordnmbl/opencap-processing/blob/main/OpenSimPipeline/JointReaction/computeJointLoading.py
% % starting from line 638
% 
% % we account for the size of the femur and multiply the
% % intercondylar distance with the applied scaling factor for
% % the femur (in med/lat direction)
% 
% %% Input Variables
% ScaleFactorsFile = 'E:\OSS\RajagopalLaiUhlrich2023-standard-tf-2\Scaling\RajagopalLaiUhlrich2023-ScaleFactors.xml';
% jrlFile = 'E:\OSS\RajagopalLaiUhlrich2023-standard-tf-2\Simulation\standard-tf-2_WalkA02_1_r\analyze\standard-tf-2_WalkA02_1_r_JointReaction_ReactionLoads.sto';
% 
% %% open applied scale factors and find femur_r and femur_l factors
% fid = fopen(ScaleFactorsFile, 'r'); %  fopen(fullfile(sessionPath, 'OpenSim', 'applied_scaleset.xml'), 'r');
% lines = {};
% tline = fgetl(fid);
% while ischar(tline)
%     lines{end+1} = tline;
%     tline = fgetl(fid);
% end
% fclose(fid);
% 
% for i = 2:length(lines)
%     if contains(lines{i}, '<segment>femur_r</segment>')
%         scaleLine = strtrim(lines{i-1});
%         if contains(scaleLine, '<scales>')
%             nums = regexp(scaleLine, '<scales>\s*([\d\.\-eE\s]+)\s*</scales>', 'tokens');
%             if ~isempty(nums)
%                 scale_femur_r = str2double(strsplit(strtrim(nums{1}{1})));
%             else
%                 warning('No scale values found in the expected format.');
%             end
%         end
%     end
%     if contains(lines{i}, '<segment>femur_l</segment>')
%         scaleLine = strtrim(lines{i-1});
%         if contains(scaleLine, '<scales>')
%             nums = regexp(scaleLine, '<scales>\s*([\d\.\-eE\s]+)\s*</scales>', 'tokens');
%             if ~isempty(nums)
%                 scale_femur_l = str2double(strsplit(strtrim(nums{1}{1})));
%             else
%                 warning('No scale values found in the expected format.');
%             end
%         end
%     end
% end
% 
% % Code from OpenCap --> this is just a very rough approximation
% base_intercondylarDistance = 0.04; % (Lerner 2015)
% subjectspecific_intercondylarDistance_l = base_intercondylarDistance * scale_femur_l(3); % use value of z direction, i.e. med/lat
% subjectspecific_intercondylarDistance_r = base_intercondylarDistance * scale_femur_r(3); % use value of z direction, i.e. med/lat
% 
% tmp_data = load_sto_file(jrlFile);
% 
% % code from opencap processing - see comment above
% KAM_r = tmp_data.walker_knee_r_on_tibia_r_in_tibia_r_mx;
% KAM_l = tmp_data.walker_knee_l_on_tibia_l_in_tibia_l_mx;
% Fy_r = tmp_data.walker_knee_r_on_tibia_r_in_tibia_r_fy;
% Fy_l = tmp_data.walker_knee_l_on_tibia_l_in_tibia_l_fy;
% 
% MCF_r = -Fy_r/2 - KAM_r/subjectspecific_intercondylarDistance_r;
% MCF_l = -Fy_l/2 + KAM_l/subjectspecific_intercondylarDistance_l;
% 
