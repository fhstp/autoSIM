function subResults = analyzeSubjectFiles(workingDirectory, rootDirectory, path2setupFiles, condition, prefix, timeNormFlag, trialType, Model2Use)

%% Collect and analyze Simulation output from one subject
% This function searches for all relevant output files for one
% subject and single condition and creates summarized plots and a results
% summary file.
%
% INPUT:
%   - workingDirectory: full path of the working directory where all basic
%     c3d files are stored
%   - rootDirectory: e.g.  top level dir containing alls working dirs
%   - sujectname: string, names of the subject (e.g. 'ID01 or 'UncleSam')
%   - condition: string, name of the condition (e.g. walking, dynamic, ...)
%   - prefix: string, used prefix (e.g. 'test')
%   - indicates if all data will be time-normalized to 100%

% Written by:           Brian Horsak - brian.horsak@fhstp.ac.at
% Acknowledgements:     Many thx to Bryce Killian (KU Leuven) for his great
%                       support
%
% Last changed:         06/2025
% -------------------------------------------------------------------------

%% Some house keeping variables
line_width = 1.2;
path.NormKCF = strcat(path2setupFiles,'..\..\_commonFiles\normData\OrthoLoad\OrthoLoad_KCFnorm.mat');
path.NormEMG = strcat(path2setupFiles,'..\..\_commonFiles\normData\Lencioni_2019\Lencioni2019_EMGnorm.mat');

%% Define variables of interest

% Norm muscles
msls_norm = {'Soleus','Gastrocnemius_Medialis','Rectus_Femoris', ...
    'Gluteus_Maximus', 'Biceps_Femoris', 'Biceps_Femoris',  ...
    'Peroneus_Longus','Vastus_Medialis','Tibialis_Anterior'};

% GRF vars 2 plot
GRFVars2plot = {
    'calcn_#_FP1_Fx', 'calcn_#_FP1_Fy', 'calcn_#_FP1_Fz', ...
    'calcn_#_FP2_Fx', 'calcn_#_FP2_Fy', 'calcn_#_FP2_Fz', ...
    'calcn_#_FP3_Fx', 'calcn_#_FP3_Fy', 'calcn_#_FP3_Fz', ...
    'calcn_#_FP4_Fx', 'calcn_#_FP4_Fy', 'calcn_#_FP4_Fz', ...
    'calcn_#_FP5_Fx', 'calcn_#_FP5_Fy', 'calcn_#_FP5_Fz', ...
    'calcn_#_FP6_Fx', 'calcn_#_FP6_Fy', 'calcn_#_FP6_Fz', ...
    'calcn_#_FP7_Fx', 'calcn_#_FP7_Fy', 'calcn_#_FP7_Fz', ...
    'calcn_#_FP8_Fx', 'calcn_#_FP8_Fy', 'calcn_#_FP8_Fz', ...
    'calcn_#_FP9_Fx', 'calcn_#_FP9_Fy', 'calcn_#_FP9_Fz'};

switch Model2Use
    case {'lernergopal'}

        % Kinematics vars 2 plot
        kinemVars2plot = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', ...
            'hip_flexion_#', 'hip_adduction_#', 'hip_rotation_#', ...
            'knee_angle_#', 'ankle_angle_#', 'subtalar_angle_#', 'mtp_angle_#'};

        % Moment vars 2 plot
        InvDynVars2plot = {'pelvis_tilt_moment', 'pelvis_list_moment', 'pelvis_rotation_moment', ...
            'hip_flexion_#_moment', 'hip_adduction_#_moment', 'hip_rotation_#_moment', ...
            'knee_angle_#_moment', 'ankle_angle_#_moment', 'subtalar_angle_#_moment', 'mtp_angle_#_moment'};

        % Select muscle activations
        msls = {'addbrev_', 'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', ...
            'bflh_', 'bfsh_', 'edl_', 'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', ...
            'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', 'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', ...
            'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', 'recfem_', 'sart_', 'semimem_', ...
            'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_', 'addbrev_', ...
            'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', 'bflh_', 'bfsh_', 'edl_', ...
            'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', 'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', ...
            'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', 'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', ...
            'recfem_', 'sart_', 'semimem_', 'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_'};

        % Muscle subset for which I have normdata
        msls_subNorm = {'soleus_','gasmed_','recfem_', ...
            'glmax1_', 'bfsh_', 'bflh_', ...
            'perlong_', 'vasmed_', 'tibant_'};

        % Select contact forces (3)
        % NOTE: this might need adjustment if some seetings were changed in the set up files.
        tmp2Replace = '<##>'; % this is later needed to replace with the side tag
        %tf_contactF = {['femur_weld_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fx'],['femur_weld_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fy'],['femur_weld_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fz']};
        %tf_contactF = {['knee_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fx'],['knee_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fy'],['knee_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fz']};
        %tf_contactF = {['femur_weld_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fx'],['femur_weld_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fy'],['femur_weld_',tmp2Replace,'_on_femoral_cond_',tmp2Replace,'_in_femoral_cond_',tmp2Replace,'_fz']};
        tf_contactF = {['Lerner_knee_', tmp2Replace, '_on_sagittal_articulation_frame_', tmp2Replace, '_in_sagittal_articulation_frame_', tmp2Replace, '_fx'], ['Lerner_knee_', tmp2Replace, '_on_sagittal_articulation_frame_', tmp2Replace, '_in_sagittal_articulation_frame_', tmp2Replace, '_fy'], ['Lerner_knee_', tmp2Replace, '_on_sagittal_articulation_frame_', tmp2Replace, '_in_sagittal_articulation_frame_', tmp2Replace, '_fz']};
    
    case 'rajagopal'

        % Kinematics vars 2 plot
        kinemVars2plot = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', ...
            'hip_flexion_#', 'hip_adduction_#', 'hip_rotation_#', ...
            'knee_angle_#', 'ankle_angle_#', 'subtalar_angle_#', 'mtp_angle_#'};

        % Moment vars 2 plot
        InvDynVars2plot = {'pelvis_tilt_moment', 'pelvis_list_moment', 'pelvis_rotation_moment', ...
            'hip_flexion_#_moment', 'hip_adduction_#_moment', 'hip_rotation_#_moment', ...
            'knee_angle_#_moment', 'ankle_angle_#_moment', 'subtalar_angle_#_moment', 'mtp_angle_#_moment'};

        % Select muscle activations
        msls = {'addbrev_', 'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', ...
            'bflh_', 'bfsh_', 'edl_', 'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', ...
            'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', 'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', ...
            'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', 'recfem_', 'sart_', 'semimem_', ...
            'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_', 'addbrev_', ...
            'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', 'bflh_', 'bfsh_', 'edl_', ...
            'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', 'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', ...
            'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', 'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', ...
            'recfem_', 'sart_', 'semimem_', 'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_'};

        % Muscle subset for which I have normdata
        msls_subNorm = {'soleus_','gasmed_','recfem_', ...
            'glmax1_', 'bfsh_', 'bflh_', ...
            'perlong_', 'vasmed_', 'tibant_'};

        % Select contact forces (3)
        tmp2Replace = '<##>'; % this is later needed to replace with the side tag
        tf_contactF = {['walker_knee_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fx'], ['walker_knee_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fy'],['walker_knee_',tmp2Replace,'_on_femur_',tmp2Replace,'_in_femur_',tmp2Replace,'_fz']};

    case 'RajagopalLaiUhlrich2023' % in theory this is similar to rajagopal. BUT, there I have used express force in parent instead of child resulting in "walker_knee_r_on_tibia_r_in_tibia_r_" instead of walker_knee_r_on_femur_r_in_femur_r_.

        % Kinematics vars 2 plot
        kinemVars2plot = {'pelvis_tilt', 'pelvis_list', 'pelvis_rotation', ...
            'hip_flexion_#', 'hip_adduction_#', 'hip_rotation_#', ...
            'knee_angle_#', 'ankle_angle_#', 'subtalar_angle_#', 'mtp_angle_#'};

        % Moment vars 2 plot
        InvDynVars2plot = {'pelvis_tilt_moment', 'pelvis_list_moment', 'pelvis_rotation_moment', ...
            'hip_flexion_#_moment', 'hip_adduction_#_moment', 'hip_rotation_#_moment', ...
            'knee_angle_#_moment', 'ankle_angle_#_moment', 'subtalar_angle_#_moment', 'mtp_angle_#_moment'};

        % Select muscle activations
        msls = {'addbrev_', 'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', ...
            'bflh_', 'bfsh_', 'edl_', 'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', ...
            'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', 'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', ...
            'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', 'recfem_', 'sart_', 'semimem_', ...
            'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_', 'addbrev_', ...
            'addlong_', 'addmagDist_', 'addmagIsch_', 'addmagMid_', 'addmagProx_', 'bflh_', 'bfsh_', 'edl_', ...
            'ehl_', 'fdl_', 'fhl_', 'gaslat_', 'gasmed_', 'glmax1_', 'glmax2_', 'glmax3_', 'glmed1_', 'glmed2_', ...
            'glmed3_', 'glmin1_', 'glmin2_', 'glmin3_', 'grac_', 'iliacus_', 'perbrev_', 'perlong_', 'piri_', 'psoas_', ...
            'recfem_', 'sart_', 'semimem_', 'semiten_', 'soleus_', 'tfl_', 'tibant_', 'tibpost_', 'vasint_', 'vaslat_', 'vasmed_'};

        % Muscle subset for which I have normdata
        msls_subNorm = {'soleus_','gasmed_','recfem_', ...
            'glmax1_', 'bfsh_', 'bflh_', ...
            'perlong_', 'vasmed_', 'tibant_'};

        % Select contact forces (3)
        tmp2Replace = '<##>'; % this is later needed to replace with the side tag
        tf_contactF = {['walker_knee_',tmp2Replace,'_on_tibia_',tmp2Replace,'_in_tibia_',tmp2Replace,'_fx'], ['walker_knee_',tmp2Replace,'_on_tibia_',tmp2Replace,'_in_tibia_',tmp2Replace,'_fy'],['walker_knee_',tmp2Replace,'_on_tibia_',tmp2Replace,'_in_tibia_',tmp2Replace,'_fz']};
end


%% Load InputData File
InputFile = dir(strcat(workingDirectory,'Simulation\','*', condition,'*InputData*.*'));
load(horzcat(InputFile.folder,'\', InputFile.name));

% Load Model MetaData Info
SubjInfo = dir(strcat(workingDirectory,'Simulation\','*', condition,'*SubjInfo*.*'));
load(horzcat(SubjInfo.folder,'\', SubjInfo.name));

% Get file names
FileNames_InputData = fieldnames(InputData);
subjectName = char(strrep(InputData.(FileNames_InputData{1}).subjectName,' ','_'));

%% Find all Joint Reaction files
cd(strcat(workingDirectory,'Simulation\'));
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListForce = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*JointReaction_ReactionLoads.sto')));
end

filesForce = {}; % initialize cell
for i = 1 : length(FileListForce)
    filesForce{i,1} = horzcat(FileListForce(i).folder,'\', FileListForce(i).name);
end

%% Find all muscle activation files
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListMuscleAct = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*StaticOptimization_activation.sto')));
end

filesMuscleAct = {}; % initialize cell
for i = 1 : length(FileListMuscleAct)
    filesMuscleAct{i,1} = horzcat(FileListMuscleAct(i).folder,'\', FileListMuscleAct(i).name);
end

%% Find all muscle forces files
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListMuscleForce = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*StaticOptimization_force.sto')));
end

filesMuscleForce = {}; % initialize cell
for i = 1 : length(FileListMuscleForce)
    filesMuscleForce{i,1} = horzcat(FileListMuscleForce(i).folder,'\', FileListMuscleForce(i).name);
end

%% Find all inverse-dynamics files
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListInverseDynamics = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*inverse-dynamics.sto')));
end

filesInverseDynamics = {}; % initialize cell
for i = 1 : length(FileListInverseDynamics)
    filesInverseDynamics{i,1} = horzcat(FileListInverseDynamics(i).folder,'\', FileListInverseDynamics(i).name);
end


%% Find all IK kinematics files
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListInvKinematics = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*IK_motion_file.mot')));
end

filesInvKinematics = {}; % initialize cell
for i = 1 : length(FileListInvKinematics)
    filesInvKinematics{i,1} = horzcat(FileListInvKinematics(i).folder,'\', FileListInvKinematics(i).name);
end

%% Find all IK error location results
if strcmp(prefix, '') % if prefix is empyt
else
    FileListErrorIKLocs = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*ErrorStruct.mat')));
end

filesErrorIKLocs = {}; % initialize cell
for i = 1 : length(FileListErrorIKLocs)
    filesErrorIKLocs{i,1} = horzcat(FileListErrorIKLocs(i).folder,'\', FileListErrorIKLocs(i).name);
end

%% Find all Force Reporter files
cd(strcat(workingDirectory,'Simulation\'));
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListForceReporter = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*ForceReporter_forces.sto')));

end

filesForceReporter = {}; % initialize cell
for i = 1 : length(FileListForceReporter)
    filesForceReporter{i,1} = horzcat(FileListForceReporter(i).folder,'\', FileListForceReporter(i).name);
end

%% Find all Body Kinematics files
cd(strcat(workingDirectory,'Simulation\'));
Folder = cd;
if strcmp(prefix, '') % if prefix is empyt
else
    FileListBodyKinematics = dir(fullfile(Folder, '**', strcat('*',prefix,'*',condition,'*BodyKinematics_pos_global.sto')));
end

filesBodyKinematics = {}; % initialize cell
for i = 1 : length(FileListBodyKinematics)
    filesBodyKinematics{i,1} = horzcat(FileListBodyKinematics(i).folder,'\', FileListBodyKinematics(i).name);
end


%% Check if all file lists have the same size
% Deal unequal file lists

tmp = strrep(filesInvKinematics, strcat(workingDirectory,'Simulation\'), ''); % start with the results which has highest chances to have results. This will help that I do not miss data during dealing!

cnt = 1;
for k = 1 : length(tmp)
    idxSlash = strfind(tmp{k},'\');
    idxSlash = idxSlash(1) -1;
    currentTrialName = tmp{k}(1:idxSlash);

    if any(contains(filesForce, currentTrialName)) ...
            && any(contains(filesMuscleAct, currentTrialName)) && ...
            any(contains(filesMuscleForce, currentTrialName)) && ...
            any(contains(filesInvKinematics, currentTrialName)) && ...
            any(contains(filesInverseDynamics, currentTrialName)) && ...
            any(contains(filesErrorIKLocs, currentTrialName)) && ...
            any(contains(filesForceReporter, currentTrialName)) && ...
            any(contains(filesBodyKinematics, currentTrialName))

        % Deale lists
        filesForce_dealed{cnt,1} = filesForce{find(contains(filesForce, currentTrialName))};
        filesMuscleAct_dealed{cnt,1} = filesMuscleAct{find(contains(filesMuscleAct, currentTrialName))};
        filesMuscleForce_dealed{cnt,1} = filesMuscleForce{find(contains(filesMuscleForce, currentTrialName))};
        filesInvKinematics_dealed{cnt,1} = filesInvKinematics{find(contains(filesInvKinematics, currentTrialName))};
        filesInverseDynamics_dealed{cnt,1} = filesInverseDynamics{find(contains(filesInverseDynamics, currentTrialName))};
        filesErrorIKLocs_dealed{cnt,1} = filesErrorIKLocs{find(contains(filesErrorIKLocs, currentTrialName))};
        filesForceReporter_dealed{cnt,1} = filesForceReporter{find(contains(filesForceReporter, currentTrialName))};
        filesBodyKinematics_dealed{cnt,1} = filesBodyKinematics{find(contains(filesBodyKinematics, currentTrialName))};

        % Increase cnt
        cnt = cnt + 1;
    end
end

% check again to make sure
lengths2check = [length(filesForce_dealed), length(filesMuscleAct_dealed), length(filesMuscleForce_dealed), length(filesInvKinematics_dealed), ...
    length(filesInverseDynamics_dealed), length(filesForceReporter_dealed), length(filesBodyKinematics_dealed)];

if all(diff(lengths2check) == 0)
    disp('>>>>> File lists looks ok!')
else
    warning('Number of files is still uneuqBodyKinematicsal for results categories! Please check!')
    return
end


%% Initialize the results struct
resultsAll_perTrial = struct(); % initialize struct
results_perSub = struct(); % initialize struct
sideS = {'r','l'};
for ss = 1 : 2
    results_perSub.(sideS{ss}) = struct(); % initialize
    results_perSub.(sideS{ss}).BodyKinematics = struct(); % initialize
    results_perSub.(sideS{ss}).ForceReporter = struct(); % initialize
    results_perSub.(sideS{ss}).ContactForces = struct(); % initialize
    results_perSub.(sideS{ss}).MuscleActivations = struct(); % initialize
    results_perSub.(sideS{ss}).MuscleForces = struct(); % initialize
    results_perSub.(sideS{ss}).InverseKinematics = struct(); % initialize
    results_perSub.(sideS{ss}).InverseDynamics = struct(); % initialize
    %results_perSub.(sideS{ss}).MetaData = struct(); % initialize
    results_perSub.(sideS{ss}).IKerrors = struct(); % initialize
    results_perSub.(sideS{ss}).GRFs = struct(); % initialize
    results_perSub.(sideS{ss}).Events = struct(); % initialize
    results_perSub.(sideS{ss}).ReserveActivation = struct(); % initialize
    results_perSub.(sideS{ss}).ForceReserve = struct(); % initialize
    results_perSub.(sideS{ss}).ResidualForcesPelvis = struct(); % initialize
    % Note that the metaData (== subjInfo) are added at the end of the script
end


%% Collect all relevant data

% Get nFrames over all trials
standardSizeArray = 500; % set standard array size for padding when timeNorm is off.

nFramesfilesForce = standardSizeArray; %getNumFrames(filesForce_dealed);
nFramesfilesMuscleAct = standardSizeArray; %getNumFrames(filesMuscleAct_dealed);
nFramesfilesMuscleForce = standardSizeArray; %getNumFrames(filesMuscleForce_dealed);
nFramesfilesInvKinematics = standardSizeArray; %getNumFrames(filesInvKinematics_dealed);
nFramesfilesInverseDynamics = standardSizeArray; %getNumFrames(filesInverseDynamics_dealed);
nFramesfilesForceReporter = standardSizeArray; %getNumFrames(filesForceReporter_dealed);
nFramesfilesBodyKinematics = standardSizeArray; %getNumFrames(filesBodyKinematics_dealed);


% Now read the data
for i = 1 : length(filesInvKinematics_dealed)

    % ONLY for RajagopalLaiUhlrich2023 - add the computed the Med. and Lat. KJCF to the *.sto file.
    if strcmp(Model2Use, 'RajagopalLaiUhlrich2023')
        
        % Set file paths.
        ScaleFactorsFile = fullfile(workingDirectory, 'Scaling', strcat(Model2Use,'-ScaleFactors.xml'));
        JointReactionLoadsFile = filesForce_dealed{i};

        % Compute the medial and lateral KJCF
        [MCFy_r, MCFy_l, LCFy_r, LCFy_l] = compute_MKCF(ScaleFactorsFile, JointReactionLoadsFile);        

        % Now add data to sto file.
        newVarData = [MCFy_r, MCFy_l, LCFy_r, LCFy_l];
        addVariablesToSTO(JointReactionLoadsFile, {'MCFy_r', 'MCFy_l', 'LCFy_r', 'LCFy_l'}, newVarData)

    end
    
    % Read all files containing results
    [contactForce_data, contactForce_labels, contactForce_header] = read_opensim_mot(filesForce_dealed{i});
    [muscleAct_data, muscleAct_labels, muscleAct_header] = read_opensim_mot(filesMuscleAct_dealed{i});
    [muscleForce_data, muscleForce_labels, muscleForce_header] = read_opensim_mot(filesMuscleForce_dealed{i});
    [inverseKinematics_data, inverseKinematics_labels, inverseKinematics_header] = read_opensim_mot(filesInvKinematics_dealed{i});
    [inverseDyn_data, inverseDyn_labels, inverseDyn_header] = read_opensim_mot(filesInverseDynamics_dealed{i});
    [forceReporter_data, forceReporter_labels, forceReporter_header] = read_opensim_mot(filesForceReporter_dealed{i});
    [bodyKinematics_data, bodyKinematics_labels, bodyKinematics_header] = read_opensim_mot(filesBodyKinematics_dealed{i});

    tmpError = load(filesErrorIKLocs_dealed{i});
    IK_errorLocs = tmpError.errors;

    % Extract current file name of the loop and delete trailing and leading "_". Also get rid of prefix
    stop_strIdx = strfind(filesForce_dealed{i} ,'analyze')-2;
    tmp = strfind(filesForce_dealed{i}(1:stop_strIdx),'\');
    start_strIdx = tmp(end)+1;
    trial_name_tmp = strip(erase(filesForce_dealed{i}(start_strIdx:stop_strIdx),prefix),'_');
    name_tmp = trial_name_tmp(strfind(trial_name_tmp,condition):end); %get rid of prefix if prefix was not defined specifically

    % Find in InputData trial names
    for k = 1 : length(FileNames_InputData)
        if contains(FileNames_InputData{k},name_tmp);
            trialNameInputData = FileNames_InputData{k};
        end
    end

    % Get rid of forces from side which is not part of the analysis
    contactForce_labels_side = {}; % initialize
    contactForce_data_side = []; % initialize
    side_tmp = char(lower(InputData.(trialNameInputData).Side(1)));

    % Define side to ignore
    if strcmp(side_tmp, 'l')
        Side2ignore = 'r';
    elseif strcmp(side_tmp, 'r')
        Side2ignore = 'l';
    end

    cnt = 1;
    for s = 1 : length(contactForce_labels)
        if contains(contactForce_labels{s}, ['_',Side2ignore, '_'])
            % do nothing
        else
            % Only get vars of interest
            contactForce_labels_side(cnt,1) = contactForce_labels(s);
            contactForce_data_side(:,cnt) = contactForce_data(:,s);
            cnt = cnt + 1;
        end

    end

    % Update  vars
    contactForce_labels = contactForce_labels_side;
    contactForce_data = contactForce_data_side;

    % Shorten Var Names, Matlab only allows table names with max. of 65 characters
    contactForce_short = contactForce_labels;
    pat = {'articulation', 'sagittal'};
    new = {'art', 'sag'};
    for k = 1 : length(contactForce_short)
        for j = 1 : length(pat)
            contactForce_short{k,1} = strrep(contactForce_short{k}, pat{j}, new{j});
        end
    end

    % We also need to shorten the search strings.
    for k = 1 : length(tf_contactF)
        for j = 1 : length(pat)
            tf_contactF{k} = strrep(tf_contactF{k}, pat{j}, new{j});
        end
    end

    % Update  vars again after shorting the long one
    contactForce_labels = contactForce_short;

    % Test if <trial_name_tmp> starts with a number. If so add string to the start,
    % otherwise it is no a valid struct fieldname
    if isstrprop(trial_name_tmp(1),'digit')
        trial_name_tmp = strcat('trial_',trial_name_tmp);
    end

    % Update header
    contactForce_header.nColumns = length(contactForce_data_side);

    % Contact Forces
    resultsAll_perTrial.(trial_name_tmp).dataForceData = contactForce_data;
    resultsAll_perTrial.(trial_name_tmp).labelsForceData = contactForce_labels;
    resultsAll_perTrial.(trial_name_tmp).headerForceData = contactForce_header;
    resultsAll_perTrial.(trial_name_tmp).ForceDataTable = cell2table(num2cell(contactForce_data),'VariableNames',contactForce_labels');

    % Muscle activation
    resultsAll_perTrial.(trial_name_tmp).dataMuscleAct = muscleAct_data;
    resultsAll_perTrial.(trial_name_tmp).labelsMuscleAct = muscleAct_labels;
    resultsAll_perTrial.(trial_name_tmp).headerMuscleAct = muscleAct_header;
    resultsAll_perTrial.(trial_name_tmp).MuscleActTable = cell2table(num2cell(muscleAct_data),'VariableNames',muscleAct_labels');

    % Muscle force
    resultsAll_perTrial.(trial_name_tmp).dataMuscleForce = muscleForce_data;
    resultsAll_perTrial.(trial_name_tmp).labelsMuscleForce = muscleForce_labels;
    resultsAll_perTrial.(trial_name_tmp).headerMuscleForce = muscleForce_header;
    resultsAll_perTrial.(trial_name_tmp).MuscleForceTable = cell2table(num2cell(muscleForce_data),'VariableNames',muscleForce_labels');

    % Inverse Kinematics
    resultsAll_perTrial.(trial_name_tmp).dataInverseKinematics = inverseKinematics_data;
    resultsAll_perTrial.(trial_name_tmp).labelsInverseKinematics = inverseKinematics_labels;
    resultsAll_perTrial.(trial_name_tmp).headerInverseKinematics = inverseKinematics_header;
    resultsAll_perTrial.(trial_name_tmp).InverseKinematicsTable = cell2table(num2cell(inverseKinematics_data),'VariableNames',inverseKinematics_labels');

    % Inverse Dynamics
    resultsAll_perTrial.(trial_name_tmp).datainverseDyn = inverseDyn_data;
    resultsAll_perTrial.(trial_name_tmp).labelsinverseDyn = inverseDyn_labels;
    resultsAll_perTrial.(trial_name_tmp).headerinverseDyn = inverseDyn_header;
    resultsAll_perTrial.(trial_name_tmp).InverseDynamicsTable = cell2table(num2cell(inverseDyn_data),'VariableNames',inverseDyn_labels');

    % Force Reporter
    resultsAll_perTrial.(trial_name_tmp).dataforceReporter = forceReporter_data;
    resultsAll_perTrial.(trial_name_tmp).labelsforceReporter = forceReporter_labels;
    resultsAll_perTrial.(trial_name_tmp).headerforceReporter = forceReporter_header;
    resultsAll_perTrial.(trial_name_tmp).ForceReporterTable = cell2table(num2cell(forceReporter_data),'VariableNames',forceReporter_labels');

    % Body Kinematics
    resultsAll_perTrial.(trial_name_tmp).databodyKinematics = bodyKinematics_data;
    resultsAll_perTrial.(trial_name_tmp).labelsbodyKinematics = bodyKinematics_labels;
    resultsAll_perTrial.(trial_name_tmp).headerbodyKinematics = bodyKinematics_header;
    resultsAll_perTrial.(trial_name_tmp).BodyKinematicsTable = cell2table(num2cell(bodyKinematics_data),'VariableNames',bodyKinematics_labels');

    % IK error locs
    resultsAll_perTrial.(trial_name_tmp).IKerror = IK_errorLocs.IKerror;
    resultsAll_perTrial.(trial_name_tmp).ReserveActivation = IK_errorLocs.ReserveActivation;
    resultsAll_perTrial.(trial_name_tmp).ForceReserve = IK_errorLocs.ForceReserve;
    resultsAll_perTrial.(trial_name_tmp).IKerrorLocs = IK_errorLocs.IKerrorLocs;
    resultsAll_perTrial.(trial_name_tmp).ResidualForcesPelvis = IK_errorLocs.ResidualForcesPelvis;

    %% Write data to per subject struct

    % Write trial info
    if ~isfield(results_perSub.(side_tmp), 'trialName')
        results_perSub.(side_tmp).trialName{1,1} = name_tmp;
    else
        results_perSub.(side_tmp).trialName{1,end+1} = name_tmp;
    end

    % Collect ReserveActivation
    numFrames = size(IK_errorLocs.ReserveActivation,1); % get number of frames
    ReserveActivationVars = IK_errorLocs.ReserveActivation.Properties.VariableNames;
    for j = 1 : length(ReserveActivationVars)

        if ~isfield(results_perSub.(side_tmp).ReserveActivation, (ReserveActivationVars{j}))
            if timeNormFlag
                results_perSub.(side_tmp).ReserveActivation.(ReserveActivationVars{j})(:,1) = interpft(IK_errorLocs.ReserveActivation.(ReserveActivationVars{j}),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).ReserveActivation.(ReserveActivationVars{j})(:,1) = padarray(IK_errorLocs.ReserveActivation.(ReserveActivationVars{j}), nFramesfilesForce-numFrames,0,'post');
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).ReserveActivation.(ReserveActivationVars{j})(:,end+1) = interpft(IK_errorLocs.ReserveActivation.(ReserveActivationVars{j}),101); % add data always to the right
            else
                results_perSub.(side_tmp).ReserveActivation.(ReserveActivationVars{j})(:,end+1) = padarray(IK_errorLocs.ReserveActivation.(ReserveActivationVars{j}), nFramesfilesForce-numFrames,0,'post');
            end
        end
    end

    % Collect ForceReserve
    numFrames = size(IK_errorLocs.ForceReserve,1); % get number of frames
    ForceReserveVars = IK_errorLocs.ForceReserve.Properties.VariableNames;
    for j = 1 : length(ForceReserveVars)

        if ~isfield(results_perSub.(side_tmp).ForceReserve, (ForceReserveVars{j}))
            if timeNormFlag
                results_perSub.(side_tmp).ForceReserve.(ForceReserveVars{j})(:,1) = interpft(IK_errorLocs.ForceReserve.(ForceReserveVars{j}),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).ForceReserve.(ForceReserveVars{j})(:,1) = padarray(IK_errorLocs.ForceReserve.(ForceReserveVars{j}), nFramesfilesForce-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).ForceReserve.(ForceReserveVars{j})(:,end+1) = interpft(IK_errorLocs.ForceReserve.(ForceReserveVars{j}),101); % add data always to the right
            else
                results_perSub.(side_tmp).ForceReserve.(ForceReserveVars{j})(:,end+1) = padarray(IK_errorLocs.ForceReserve.(ForceReserveVars{j}), nFramesfilesForce-numFrames,0,'post'); %create variable sthe first time
            end
        end
    end

    % Collect ResidualForcesPelvis
    numFrames = size(IK_errorLocs.ResidualForcesPelvis,1); % get number of frames
    ResidualForcesPelvisVars = IK_errorLocs.ResidualForcesPelvis.Properties.VariableNames;
    for j = 1 : length(ResidualForcesPelvisVars)

        if ~isfield(results_perSub.(side_tmp).ResidualForcesPelvis, (ResidualForcesPelvisVars{j}))
            if timeNormFlag
                results_perSub.(side_tmp).ResidualForcesPelvis.(ResidualForcesPelvisVars{j})(:,1) = interpft(IK_errorLocs.ResidualForcesPelvis.(ResidualForcesPelvisVars{j}),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).ResidualForcesPelvis.(ResidualForcesPelvisVars{j})(:,1) = padarray(IK_errorLocs.ResidualForcesPelvis.(ResidualForcesPelvisVars{j}), nFramesfilesForce-numFrames,0,'post');
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).ResidualForcesPelvis.(ResidualForcesPelvisVars{j})(:,end+1) = interpft(IK_errorLocs.ResidualForcesPelvis.(ResidualForcesPelvisVars{j}),101); % add data always to the right
            else
                results_perSub.(side_tmp).ResidualForcesPelvis.(ResidualForcesPelvisVars{j})(:,end+1) = padarray(IK_errorLocs.ResidualForcesPelvis.(ResidualForcesPelvisVars{j}), nFramesfilesForce-numFrames,0,'post');
            end
        end
    end

    % Collect IK errors
    numFrames = size(IK_errorLocs.IKerrorLocs,1); % get number of frames
    IKerrorsLocsMrksVars = IK_errorLocs.IKerrorLocs.Properties.VariableNames;
    for j = 1 : length(IKerrorsLocsMrksVars)

        if ~isfield(results_perSub.(side_tmp).IKerrors, (IKerrorsLocsMrksVars{j}))
            if timeNormFlag
                results_perSub.(side_tmp).IKerrors.(IKerrorsLocsMrksVars{j})(:,1) = interpft(IK_errorLocs.IKerrorLocs.(IKerrorsLocsMrksVars{j}),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).IKerrors.(IKerrorsLocsMrksVars{j})(:,1) = padarray(IK_errorLocs.IKerrorLocs.(IKerrorsLocsMrksVars{j}), nFramesfilesForce-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).IKerrors.(IKerrorsLocsMrksVars{j})(:,end+1) = interpft(IK_errorLocs.IKerrorLocs.(IKerrorsLocsMrksVars{j}),101); % add data always to the right
            else
                results_perSub.(side_tmp).IKerrors.(IKerrorsLocsMrksVars{j})(:,end+1) = padarray(IK_errorLocs.IKerrorLocs.(IKerrorsLocsMrksVars{j}), nFramesfilesForce-numFrames,0,'post'); %create variable sthe first time
            end
        end
    end

    % Write Events: IC, ICi, ...
    eventVars = {'cTO','TO','IC','cIC','ICi'};

    for j = 1 : length(eventVars)
        if ~isfield(results_perSub.(side_tmp).Events, (eventVars{j}))
            results_perSub.(side_tmp).Events.(eventVars{j})(1) = InputData.(trialNameInputData).(eventVars{j});
        else
            results_perSub.(side_tmp).Events.(eventVars{j})(end+1) = InputData.(trialNameInputData).(eventVars{j});
        end
    end

    % Contact forces
    numFrames = size(contactForce_data,1); %height(contactForce_data); % get number of frames
    for j = 1:length(contactForce_labels)
        fn = char(strrep(contactForce_labels(j),'.','_'));
        if ~isfield(results_perSub.(side_tmp).ContactForces,fn)
            if timeNormFlag
                results_perSub.(side_tmp).ContactForces.(fn)(:,1) = interpft(contactForce_data(:,j),101); %create variable the first time
            else
                results_perSub.(side_tmp).ContactForces.(fn)(:,1) = padarray(contactForce_data(:,j),nFramesfilesForce-numFrames,0,'post'); % pad with nan if necessary
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).ContactForces.(fn)(:,end+1) = interpft(contactForce_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).ContactForces.(fn)(:,end+1) = padarray(contactForce_data(:,j),nFramesfilesForce-numFrames,0,'post'); % add data always to the right
            end
        end
    end

    % Muscle Activations
    numFrames = size(muscleAct_data,1); %height(muscleAct_data); % get number of frames
    for j = 1:length(muscleAct_labels)
        fn = char(strrep(muscleAct_labels(j),'/','_'));
        if strcmp(fn(1), '_'); fn = fn(2:end); end % get rid of the first '_'
        if ~isfield(results_perSub.(side_tmp).MuscleActivations,fn)
            if timeNormFlag
                results_perSub.(side_tmp).MuscleActivations.(fn)(:,1) = interpft(muscleAct_data(:,j),101);
            else
                results_perSub.(side_tmp).MuscleActivations.(fn)(:,1) = padarray(muscleAct_data(:,j),nFramesfilesMuscleAct-numFrames,0,'post');
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).MuscleActivations.(fn)(:,end+1) = interpft(muscleAct_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).MuscleActivations.(fn)(:,end+1) = padarray(muscleAct_data(:,j),nFramesfilesMuscleAct-numFrames,0,'post'); % add data always to the right
            end
        end
    end

    % Muscle Forces
    numFrames = size(muscleForce_data,1); % height(muscleForce_data); % get number of frames
    for j = 1:length(muscleForce_labels)
        fn = char(strrep(muscleForce_labels(j),'/','_'));
        if strcmp(fn(1), '_'); fn = fn(2:end); end % get rid of the first '_'
        if ~isfield(results_perSub.(side_tmp).MuscleForces,fn)
            if timeNormFlag
                results_perSub.(side_tmp).MuscleForces.(fn)(:,1) = interpft(muscleForce_data(:,j),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).MuscleForces.(fn)(:,1) = padarray(muscleForce_data(:,j),nFramesfilesMuscleForce-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).MuscleForces.(fn)(:,end+1) = interpft(muscleForce_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).MuscleForces.(fn)(:,end+1) = padarray(muscleForce_data(:,j),nFramesfilesMuscleForce-numFrames,0,'post'); % add data always to the right
            end
        end
    end


    % Inverse Kinematics
    numFrames = size(inverseKinematics_data,1); % height(inverseKinematics_data); % get number of frames
    for j = 1:length(inverseKinematics_labels)
        fn = char(inverseKinematics_labels(j));

        if ~isfield(results_perSub.(side_tmp).InverseKinematics,fn)
            if timeNormFlag
                results_perSub.(side_tmp).InverseKinematics.(fn)(:,1) = interpft(inverseKinematics_data(:,j),101); %create variable the first time
            else
                results_perSub.(side_tmp).InverseKinematics.(fn)(:,1) = padarray(inverseKinematics_data(:,j),nFramesfilesInvKinematics-numFrames,0,'post'); %create variable the first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).InverseKinematics.(fn)(:,end+1) = interpft(inverseKinematics_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).InverseKinematics.(fn)(:,end+1) = padarray(inverseKinematics_data(:,j),nFramesfilesInvKinematics-numFrames,0,'post'); % add data always to the right
            end
        end
    end

    % Inverse Dynamics
    numFrames = size(inverseDyn_data,1); % height(inverseDyn_data); % get number of frames
    for j = 1:length(inverseDyn_labels)
        fn = char(inverseDyn_labels(j));

        if ~isfield(results_perSub.(side_tmp).InverseDynamics,fn)
            if timeNormFlag
                results_perSub.(side_tmp).InverseDynamics.(fn)(:,1) = interpft(inverseDyn_data(:,j),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).InverseDynamics.(fn)(:,1) = padarray(inverseDyn_data(:,j),nFramesfilesInverseDynamics-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).InverseDynamics.(fn)(:,end+1) = interpft(inverseDyn_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).InverseDynamics.(fn)(:,end+1) = padarray(inverseDyn_data(:,j),nFramesfilesInverseDynamics-numFrames,0,'post'); % add data always to the right
            end
        end
    end

    % Force Reporter
    numFrames = size(forceReporter_data,1);
    for j = 1:length(forceReporter_labels)
        fn = char(forceReporter_labels(j));

        if ~isfield(results_perSub.(side_tmp).ForceReporter,fn)
            if timeNormFlag
                results_perSub.(side_tmp).ForceReporter.(fn)(:,1) = interpft(forceReporter_data(:,j),101); %create variable sthe first time,0,'post'
            else
                results_perSub.(side_tmp).ForceReporter.(fn)(:,1) = padarray(forceReporter_data(:,j),nFramesfilesForceReporter-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).ForceReporter.(fn)(:,end+1) = interpft(forceReporter_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).ForceReporter.(fn)(:,end+1) = padarray(forceReporter_data(:,j),nFramesfilesForceReporter-numFrames,0,'post'); % add data always to the right
            end
        end
    end

    % Force Reporter (GRFs)
    numFrames = size(forceReporter_data,1);
    idxSide = find(contains(forceReporter_labels, strcat('_', side_tmp, '_', 'FP'))); % get cells populating the force plate data of the side of interest

    % Init
    FPlabel2use = [];
    forceReporter_labels_GRF = [];
    forceReporter_data_GRF = [];

    if isempty(idxSide)
        % do nothing 4 now
    else
        % Check if one or more FPs are listed, one FP should have 9 variables 3xforce 3xcop 3xtorque
        if length(idxSide) == 9 % here I only assume one FP
            forceReporter_data_GRF = forceReporter_data(:,idxSide);
            forceReporter_labels_GRF = forceReporter_labels(idxSide);
            FPlabel2use = forceReporter_labels_GRF{1}(end-5:end-3);

        else % Here we have several FPs, only one is the correct one
            forceReporter_labels_tmp = forceReporter_labels(idxSide);
            forceReporter_data_tmp = forceReporter_data(:,idxSide);

            % Get all vars for the vertical force
            idxFy = contains(forceReporter_labels_tmp, 'Fy');
            forceReporter_labels_Fy = forceReporter_labels_tmp(idxFy);
            forceReporter_data_Fy = forceReporter_data_tmp(:,idxFy);

            % Here I assume that the FP with the greater avg Force is the correct one
            [~, idxFP] = max(mean(abs(forceReporter_data_Fy),1));
            FPlabel2use = forceReporter_labels_Fy{idxFP}(end-5:end-3);

            % Get labels for the correct FP
            idxFPsignals = find(contains(forceReporter_labels_tmp, FPlabel2use));
            forceReporter_labels_GRF = forceReporter_labels_tmp(idxFPsignals);

            % Reduce the data to match the labels array
            forceReporter_data_GRF = forceReporter_data_tmp(:,idxFPsignals);
        end

        for j = 1:length(forceReporter_labels_GRF)
            fn = forceReporter_labels_GRF{j}(end-1:end); % extract meaningful var name

            if ~isfield(results_perSub.(side_tmp).GRFs,fn)
                if timeNormFlag
                    results_perSub.(side_tmp).GRFs.(fn)(:,1) = interpft(forceReporter_data_GRF(:,j),101); %create variables the first time
                else
                    results_perSub.(side_tmp).GRFs.(fn)(:,1) = padarray(forceReporter_data_GRF(:,j),nFramesfilesForceReporter-numFrames,0,'post'); %create variable sthe first time
                end
                % Collect info which FP was used
                results_perSub.(side_tmp).GRFs.FPused{1} = FPlabel2use;
            else
                if timeNormFlag
                    results_perSub.(side_tmp).GRFs.(fn)(:,end+1) = interpft(forceReporter_data_GRF(:,j),101); % add data always to the right
                else
                    results_perSub.(side_tmp).GRFs.(fn)(:,end+1) = padarray(forceReporter_data_GRF(:,j),nFramesfilesForceReporter-numFrames,0,'post'); % add data always to the right
                end
                % Collect info which FP was used
                results_perSub.(side_tmp).GRFs.FPused{end+1} = FPlabel2use;
            end
        end
    end

    % Body Kinematics
    numFrames = size(bodyKinematics_data,1);
    for j = 1:length(bodyKinematics_labels)
        fn = char(bodyKinematics_labels(j));

        if ~isfield(results_perSub.(side_tmp).BodyKinematics,fn)
            if timeNormFlag
                results_perSub.(side_tmp).BodyKinematics.(fn)(:,1) = interpft(bodyKinematics_data(:,j),101); %create variable sthe first time
            else
                results_perSub.(side_tmp).BodyKinematics.(fn)(:,1) = padarray(bodyKinematics_data(:,j),nFramesfilesBodyKinematics-numFrames,0,'post'); %create variable sthe first time
            end
        else
            if timeNormFlag
                results_perSub.(side_tmp).BodyKinematics.(fn)(:,end+1) = interpft(bodyKinematics_data(:,j),101); % add data always to the right
            else
                results_perSub.(side_tmp).BodyKinematics.(fn)(:,end+1) = padarray(bodyKinematics_data(:,j),nFramesfilesBodyKinematics-numFrames,0,'post'); % add data always to the right
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot GRF
hfig_GRF = figure;
set(hfig_GRF,'units','centimeters','position',[0,0,29,21]);
orient(hfig_GRF,'landscape');
sgtitle(strcat('Results ForceReporter (GRFs only):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');

for i = 1:length(GRFVars2plot)

    % Plot data
    subaxis(6,3,i,'SpacingVert',0.02,'MR',0.1)
    label_r = (strrep(GRFVars2plot{i},'#', 'r'));
    label_l = (strrep(GRFVars2plot{i},'#', 'l'));

    dat_r = [];
    dat_l = [];
    % Both
    if isfield(results_perSub.('r').ForceReporter,label_r) && isfield(results_perSub.('l').ForceReporter,label_l)
        dat_r = results_perSub.('r').ForceReporter.(label_r);
        dat_l = results_perSub.('l').ForceReporter.(label_l);
        hold on
        plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        plot(dat_l,'LineWidth',line_width, 'Color', 'b', 'LineStyle', '--');
        s_tmp = '(L/R)';

        % Only right
    elseif isfield(results_perSub.('r').ForceReporter,label_r) && ~isfield(results_perSub.('l').ForceReporter,label_l)
        dat_r = results_perSub.('r').ForceReporter.(label_r);
        dat_l = nan(size(dat_r));
        hold on
        plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        s_tmp = '(R)';

        % Only left
    elseif isfield(results_perSub.('l').ForceReporter,label_l) && ~isfield(results_perSub.('r').ForceReporter,label_r)
        dat_l = results_perSub.('l').ForceReporter.(label_l);
        dat_r = nan(size(dat_l));
        hold on
        plot(dat_l,'LineWidth',line_width, 'Color', 'b');
        s_tmp = '(L)';

    end

    box on;
    title(strrep(strrep(GRFVars2plot{i}, '_',' '),'#',''));
    ylabel('N');

    % Do this only if current force plate was used
    if ~isempty(dat_r)
        numFrames = size(dat_r,1); % get number of frames, both l/r will be always created, see above
        xlim([0 numFrames]);
    end

    if i < 16
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frames or time (%)');
    end
end

%% Plot kinematics
hfig_Kinem = figure;
set(hfig_Kinem,'units','centimeters','position',[0,0,29,21]);
orient(hfig_Kinem,'landscape');
sgtitle(strcat('Results Report (Kinematics):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');

for i = 1:length(kinemVars2plot)

    % Plot data
    subaxis(5,3,i,'SpacingVert',0.02,'MR',0.1)
    label_r = (strrep(kinemVars2plot{i},'#', 'r'));
    label_l = (strrep(kinemVars2plot{i},'#', 'l'));

    % Both
    if isfield(results_perSub.('r').InverseKinematics,label_r) && isfield(results_perSub.('l').InverseKinematics,label_l)
        dat_r = results_perSub.('r').InverseKinematics.(label_r);
        dat_l = results_perSub.('l').InverseKinematics.(label_l);
        hold on
        plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        plot(dat_l,'LineWidth',line_width, 'Color', 'b', 'LineStyle', '--');
        s_tmp = '(L/R)';

        cIC_r = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO_r = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO_r =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cIC_l = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO_l = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO_l =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

        cIC = mean([cIC_r, cIC_l]);
        cTO = mean([cTO_r, cTO_l]);
        TO = mean([TO_r, TO_l]);

        % Only right
    elseif isfield(results_perSub.('r').InverseKinematics,label_r) && ~isfield(results_perSub.('l').InverseKinematics,label_l)
        dat_r = results_perSub.('r').InverseKinematics.(label_r);
        dat_l = nan(size(dat_r));
        hold on
        plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        s_tmp = '(R)';
        cIC = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time

        % Only left
    elseif isfield(results_perSub.('l').InverseKinematics,label_l) && ~isfield(results_perSub.('r').InverseKinematics,label_r)
        dat_l = results_perSub.('l').InverseKinematics.(label_l);
        dat_r = nan(size(dat_l));
        hold on
        plot(dat_l,'LineWidth',line_width, 'Color', 'b');
        s_tmp = '(L)';

        cIC = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
    end

    box on;
    title(strrep(strrep(kinemVars2plot{i}, '_',' '),'#',''));
    ylabel('degrees')
    % Add cTO, TO
    hold on;
    y_min = min(min([dat_r dat_l]));
    y_max = max(max([dat_r dat_l]));
    line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);

    numFrames = size(dat_r,1); % get number of frames
    xlim([0 numFrames]);
    if y_max > y_min
        ylim([y_min y_max]);
    end

    if i < 11
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frames or time (%)');
    end
end

%% Plot Moments
hfig_InvDyn = figure;
set(hfig_InvDyn,'units','centimeters','position',[0,0,29,21]);
orient(hfig_InvDyn,'landscape');
sgtitle(strcat('Results Report (Moments):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');

for i = 1:length(kinemVars2plot)

    % Plot data
    subaxis(5,3,i,'SpacingVert',0.02,'MR',0.1)
    label_r = (strrep(InvDynVars2plot{i},'#', 'r'));
    label_l = (strrep(InvDynVars2plot{i},'#', 'l'));

    % Both
    if isfield(results_perSub.('r').InverseDynamics,label_r) && isfield(results_perSub.('l').InverseDynamics,label_l)
        dat_r = results_perSub.('r').InverseDynamics.(label_r);
        dat_l = results_perSub.('l').InverseDynamics.(label_l);
        hold on
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b', 'LineStyle', '--');
        s_tmp = '(L/R)';

        cIC_r = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO_r = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO_r =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cIC_l = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO_l = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO_l =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

        cIC = mean([cIC_r, cIC_l]);
        cTO = mean([cTO_r, cTO_l]);
        TO = mean([TO_r, TO_l]);

        % Only right
    elseif isfield(results_perSub.('r').InverseDynamics,label_r) && ~isfield(results_perSub.('l').InverseDynamics,label_l)
        dat_r = results_perSub.('r').InverseDynamics.(label_r);
        dat_l = nan(size(dat_r));
        hold on
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        s_tmp = '(R)';
        cIC = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time

        % Only left
    elseif isfield(results_perSub.('l').InverseDynamics,label_l) && ~isfield(results_perSub.('r').InverseDynamics,label_r)
        dat_l = results_perSub.('l').InverseDynamics.(label_l);
        dat_r = nan(size(dat_l));
        hold on
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b');
        s_tmp = '(L)';

        cIC = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
    end

    box on;
    title(strrep(strrep(InvDynVars2plot{i}, '_',' '),'#',''));
    ylabel('Nm')
    % Add cTO, TO
    hold on;
    y_min = min(min([dat_r dat_l]));
    y_max = max(max([dat_r dat_l]));
    line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);

    numFrames = size(dat_r,1); % get number of frames
    xlim([0 numFrames]);
    if y_max > y_min
        ylim([y_min y_max]);
    end

    if i < 11
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frames or time (%)');
    end
end
% ---------------------------------------------------------------------

%% Plot Muscle Forces
hfig_MF = figure;
set(hfig_MF,'units','centimeters','position',[0,0,29,21]);
orient(hfig_MF,'landscape');
sgtitle(strcat('Results Report (Muscle Forces):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');
msls_names = strrep(msls,'_', ' ');

for i = 1:length(msls)

    % Plot data
    subaxis(9,9,i,'SpacingVert',0.02,'MR',0.1)
    label_r = (strcat(msls{i},'r'));
    label_l = (strcat(msls{i},'l'));

    % Both
    if isfield(results_perSub.('r').MuscleForces,label_r) && isfield(results_perSub.('l').MuscleForces,label_l)
        dat_r = results_perSub.('r').MuscleForces.(label_r);
        dat_l = results_perSub.('l').MuscleForces.(label_l);
        hold on
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b', 'LineStyle', '--');
        s_tmp = '(L/R)';

        cIC_r = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO_r = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO_r =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cIC_l = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO_l = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO_l =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

        cIC = mean([cIC_r, cIC_l]);
        cTO = mean([cTO_r, cTO_l]);
        TO = mean([TO_r, TO_l]);

        % Only right
    elseif isfield(results_perSub.('r').MuscleForces,label_r) && ~isfield(results_perSub.('l').MuscleForces,label_l)
        dat_r = results_perSub.('r').MuscleForces.(label_r);
        dat_l = nan(size(dat_r));
        hold on
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');
        s_tmp = '(R)';
        cIC = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time

        % Only left
    elseif isfield(results_perSub.('l').MuscleForces,label_l) && ~isfield(results_perSub.('r').MuscleForces,label_r)
        dat_l = results_perSub.('l').MuscleForces.(label_l);
        dat_r = nan(size(dat_l));
        hold on
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b');
        s_tmp = '(L)';

        cIC = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
    end

    box on;
    ylabel(msls_names{i});

    % Add cTO, TO
    hold on;
    y_min = min(min([dat_r dat_l]));
    y_max = max(max([dat_r dat_l]));
    hcTO = line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hTO = line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);

    numFrames = size(dat_r,1); % get number of frames
    xlim([0 numFrames]);
    ylim([y_min y_max]);

    if i < 72
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frames or time (%)');
    end
end

%% Plot Muscle Activations
hfig_MA = figure;
set(hfig_MA,'units','centimeters','position',[0,0,29,21]);
orient(hfig_MA,'landscape');
sgtitle(strcat('Results Report (sel. Muscle Activ. with Norm):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');
msls_names = strrep(msls_subNorm,'_', ' ');

% Load NormData
load(path.NormEMG);

for i = 1:length(msls_subNorm)

    % Plot data
    subplot(3,3,i);
    hold on
    label_r = (strcat(msls_subNorm{i},'r'));
    label_l = (strcat(msls_subNorm{i},'l'));

    % Both
    if isfield(results_perSub.('r').MuscleActivations,label_r) && isfield(results_perSub.('l').MuscleActivations,label_l)
        dat_r = results_perSub.('r').MuscleActivations.(label_r);
        dat_l = results_perSub.('l').MuscleActivations.(label_l);

        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r'); %[233,163,201]/255)
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b'); %, 'LineStyle', '--'); %[203,123,161]/255
        s_tmp = '(L/R)';

        cIC_r = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO_r = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO_r =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cIC_l = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO_l = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO_l =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

        cIC = mean([cIC_r, cIC_l]);
        cTO = mean([cTO_r, cTO_l]);
        TO = mean([TO_r, TO_l]);

        % Only right
    elseif isfield(results_perSub.('r').MuscleActivations,label_r) && ~isfield(results_perSub.('l').MuscleActivations,label_l)
        dat_r = results_perSub.('r').MuscleActivations.(label_r);
        dat_l = nan(size(dat_r));
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');% [233,163,201]/255);
        s_tmp = '(R)';

        cIC = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time

        % Only left
    elseif isfield(results_perSub.('l').MuscleActivations,label_l) && ~isfield(results_perSub.('r').MuscleActivations,label_r)
        dat_l = results_perSub.('l').MuscleActivations.(label_l);
        dat_r = nan(size(dat_l));
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b');% [203,123,161]/255);
        s_tmp = '(L)';

        cIC = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

    end

    box on;
    ylim([0 1]);
    ylabel(msls_names{i});

    % Add cTO, TO
    hold on;
    y_min = (min(min([dat_r dat_l])));
    y_max = (max(max([dat_r dat_l])));
    hcTO = line([cTO cTO],[0 1],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hTO  = line([TO TO],[0 1],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    line([cIC cIC],[0 1],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    ylim([0 1]);
    numFrames = size(dat_r,1); % get number of frames
    xlim([0 numFrames]); % get number of frames]);
    hold off

    % Plot Norm Values
    hold on;
    % Plot Norm Values only if walking trial
    if strcmp(trialType, 'walking')
        if strcmp(msls_norm{i},'nan')
        else
            % Amplitude normalization
            for m = 1:size(EMGnorm.(msls_norm{i}),2)
                tmp_signal = EMGnorm.(msls_norm{i})(:,m);
                emg_dat(:,m) = (tmp_signal - (min(tmp_signal))) / (((max(tmp_signal))) - (min(tmp_signal))); % Min-Max normalization
            end
            norm_mean = mean(emg_dat,2);
            sf = y_max/max(norm_mean); % calculate scaling factor to have similar peak value between norm and muscle activations
            norm_std = std(emg_dat,0,2);
            norm_std = interpft(norm_std,length(dat_r)); % apply scaling factor and normalize to data length (IC-ICi)

            norm_curve = interpft(norm_mean*sf,length(dat_r)); % apply scaling factor
            plot(norm_curve,'k','LineWidth',1.5);
            plot(norm_curve+1*norm_std,'k--','LineWidth',0.5);
            plot(norm_curve-1*norm_std,'k--','LineWidth',0.5);
        end
    end

    % Last plotting tweaks
    if i > 6
        xlabel('Frames or time (%)');
    end

end

%% Plot Contact Force from Joint Reaction
hfig_CF = figure;
set(hfig_CF,'units','centimeters','position',[5,15,30,8]);
orient(hfig_CF,'landscape');
sgtitle(strcat('Results Report (Joint Reaction Forces):',{' '}, prefix,{' - '}, condition,{' - '}, strrep(subjectName,'_',{' '}), '; R = red | L = blue'),'FontWeight','bold');

for i = 1:length(tf_contactF)

    % Plot data
    subplot(1,3,i);
    label_r = strrep(tf_contactF{i},tmp2Replace,'r');
    label_l = strrep(tf_contactF{i},tmp2Replace,'l');

    % Both
    if isfield(results_perSub.('r').ContactForces,label_r) && isfield(results_perSub.('l').ContactForces,label_l)
        dat_r = results_perSub.('r').ContactForces.(label_r)./(SubjInfo.MetaData.bodymassFromC3D*9.81); % Normalize to Bodymass * 9.81; Use only the Body mass from first node, should be the same anyway, and in case of different amount of l/r there will be otherwise a problem
        dat_l = (results_perSub.('l').ContactForces.(label_l)./(SubjInfo.MetaData.bodymassFromC3D*9.81));

        % Normalize z on the left side only to display all with same sign conventions
        if contains(label_l, 'l_fz'); dat_l = dat_l * -1; end

        hold on
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');% [221,28,119]/255);
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b');% [201,08,109]/255, 'LineStyle', '--');
        s_tmp = '(L/R)';

        cIC_r = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO_r = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO_r =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cIC_l = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO_l = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO_l =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

        cIC = mean([cIC_r, cIC_l]);
        cTO = mean([cTO_r, cTO_l]);
        TO = mean([TO_r, TO_l]);



        % Only right
    elseif isfield(results_perSub.('r').ContactForces,label_r) && ~isfield(results_perSub.('l').ContactForces,label_l)
        dat_r = results_perSub.('r').ContactForces.(label_r)./(SubjInfo.MetaData.bodymassFromC3D*9.81); % Normalize to Bodymass * 9.81;
        dat_l = nan(size(dat_r));
        hre = plot(dat_r,'LineWidth',line_width, 'Color', 'r');%  [221,28,119]/255);
        s_tmp = '(R)';

        cIC = mean((results_perSub.('r').Events.cIC - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('r').Events.cTO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('r').Events.TO - results_perSub.('r').Events.IC)./(results_perSub.('r').Events.ICi - results_perSub.('r').Events.IC))*100; % calculate in % of gait cycle time


        % Only left
    elseif isfield(results_perSub.('l').ContactForces,label_l) && ~isfield(results_perSub.('r').ContactForces,label_r)
        dat_l = results_perSub.('l').ContactForces.(label_l)./(SubjInfo.MetaData.bodymassFromC3D*9.81);

        % Normalize z on the left side only tor display all with same sign conventions
        if strcmp(label_l(end-2:end), 'z_l'); dat_l = dat_l * -1; end

        dat_r = nan(size(dat_l));
        hle = plot(dat_l,'LineWidth',line_width, 'Color', 'b');% [201,08,109]/255/255);
        s_tmp = '(L)';

        cIC = mean((results_perSub.('l').Events.cIC - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        cTO = mean((results_perSub.('l').Events.cTO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time
        TO =  mean((results_perSub.('l').Events.TO - results_perSub.('l').Events.IC)./(results_perSub.('l').Events.ICi - results_perSub.('l').Events.IC))*100; % calculate in % of gait cycle time

    end

    box on;
    ylim([0 1]);
    if i == 1
        tmpLeg = 'anterior-posterior';
    elseif i == 2
        tmpLeg = 'vertical';
    elseif i == 3
        tmpLeg = 'medio-lateral';
    end

    ylabel({char(strcat(strrep(strcat(label_r(1:2),{' '},label_r(end-16:end-2),{' ' }),'_',' '),'[Body weight]'));tmpLeg});

    % Add cTO, TO, cIC
    hold on;

    % Plot the OrthoLoad data

    % Load the norm data
    load(path.NormKCF);
    KFCfn = {'Fy','F','Fx'};
    curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),round(TO));
    %     % Adjust sign for norm data
    %     if contains(tf_contactF{i},'f_z')
    %         curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),round(TO))*-1;
    %     else
    %         curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),round(TO));
    %     end

    % Adjust sign for norm data
    %     if strcmp(KFCfn{i}, 'F')
    %         curveKFC_mean_stance = (interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),round(TO)));
    %     else
    %         curveKFC_mean_stance = interpft(mean(KFCnorm.NormBW.(KFCfn{i}),2),round(TO));
    %     end

    if strcmp(trialType, 'walking')
        % Prepare data
        curveKFC_mean_stance = (padarray(curveKFC_mean_stance,length(dat_r) - round(TO),'post')/100)*-1; % padd zeros to end, change to multipel of BW
        KFC_sd_stance = interpft(std(KFCnorm.NormBW.(KFCfn{i}),0,2),round(TO))/100; % change to multipel of BW
        KFC_sd_stance = padarray(KFC_sd_stance,length(dat_r)-round(TO),'post'); % padd zeros to end

        % Plot norm
        hOrth = plot(curveKFC_mean_stance,'k','LineWidth',line_width);
        plot(curveKFC_mean_stance+3*KFC_sd_stance,'k--','LineWidth',0.5);
        plot(curveKFC_mean_stance-3*KFC_sd_stance,'k--','LineWidth',0.5);
    end

    % Last plotting tweaks
    %if i == 1 || i == 3; y_min = -1; y_max = 1; end
    %if i == 2; y_min = 0; y_max = 5; end
    y_min = (min(min([curveKFC_mean_stance-2*KFC_sd_stance, dat_r, dat_l])));
    y_max = (max(max([curveKFC_mean_stance+2*KFC_sd_stance, dat_r, dat_l])));
    hcTO = line([cTO cTO],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    hTO = line([TO TO],[y_min y_max],'Color','k', 'LineStyle', '-', 'LineWidth', 0.5);
    line([cIC cIC],[y_min y_max],'Color','r', 'LineStyle', '--', 'LineWidth', 0.5);
    ylim([y_min y_max]);
    numFrames = size(dat_r,1); % get number of frames
    xlim([0 numFrames]);
    xlabel('Frames or time (%)');

    if strcmp(trialType, 'walking')
        if i == 3
            switch s_tmp
                case '(L/R)'
                    hleg = legend([hle(1), hre(1), hOrth, hcTO, hTO], 'Left Knee Contact Force (%BW)','Right Knee Contact Force (%BW)','OrthoLoad (3SD) [%BW])', 'cont. lat. Toe-off/IC','Toe-off', 'Location','northeast');
                case '(L)'
                    hleg = legend([hle(1), hOrth, hcTO, hTO], 'Left Knee Contact Force (%BW)','OrthoLoad (3SD) [%BW])', 'cont. lat. Toe-off/IC','Toe-off', 'Location','northeast');
                case '(R)'
                    hleg = legend([hre(1), hOrth, hcTO, hTO], 'Right Knee Contact Force (%BW)','OrthoLoad (3SD) [%BW])', 'cont. lat. Toe-off/IC','Toe-off', 'Location','northeast');
            end
            legend boxoff
            set(hleg,'FontSize',8);
        end
    else
        if i == 3
            switch s_tmp
                case '(L/R)'
                    hleg = legend([hle(1), hre(1)], 'Left Knee Contact Force (%BW)','Right Knee Contact Force (%BW)', 'Location','northeast');
                case '(L)'
                    hleg = legend([hle(1)], 'Left Knee Contact Force (%BW)', 'Location','northeast');
                case '(R)'
                    hleg = legend(hre(1), 'Right Knee Contact Force (%BW)', 'Location','northeast');
            end
            legend boxoff
            set(hleg,'FontSize',8);
        end

    end
end

%% Save file
subResults.perTrial = resultsAll_perTrial;
subResults.perSubject = results_perSub;
subResults.InputData = InputData;
subResults.SubjInfo = SubjInfo; % Add Model Personalization & SubjInfo

% Add version information
try
    path2Version = path2setupFiles;

    % Get the parent directory three levels above.
    for i = 1:3
        path2Version = fileparts(path2Version);
    end

    curVersion = readLastLine(fullfile(path2Version, 'version.md'));
    subResults.Version = curVersion;
end

% Save files
if isempty(prefix); prefix2save = strcat(prefix,'-');
else
    prefix2save = prefix;
end

% Create output folder
outputFolder = char(strcat(workingDirectory,'Simulation\allResultsFigures\'));
if ~logical(exist(outputFolder, 'dir'))
    mkdir(outputFolder)
end

% Create output folder at top level
outputFolderGroupData = char(fullfile(strcat(rootDirectory, '_', Model2Use, '-', prefix2save, '-groupData\')));
if ~logical(exist(outputFolderGroupData, 'dir'))
    mkdir(outputFolderGroupData)
end

% Make sure that no other data from same individual are already
% stored in the folder, if so add custome label to "subjectName"
ls = dir(fullfile(outputFolderGroupData, strcat(subjectName,'*', prefix2save, condition, Model2Use)));
cfiles = {ls.name};
nFiles = length(cfiles);

if nFiles == 0
    subjectName2save = subjectName;
else
    suffix = {'a', 'b', 'c', 'd', 'e', 'f', 'g'};
    subjectName2save = strcat(subjectName, suffix{nFiles});
end

saveas(hfig_GRF,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-GRFs-', Model2Use,'-Results-AllTrials.png'))));
saveas(hfig_Kinem,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-Kinem-', Model2Use,'-Results-AllTrials.png'))));
saveas(hfig_InvDyn,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-InvDyn-', Model2Use,'-Results-AllTrials.png'))));
saveas(hfig_MF,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-MF-', Model2Use,'-Results-AllTrials.png'))));
saveas(hfig_MA,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-MA-', Model2Use,'-Results-AllTrials.png'))));
saveas(hfig_CF,char(fullfile(strcat(outputFolder, subjectName,'-', prefix2save, condition,'-KCF-', Model2Use,'-Results-AllTrials.png'))));
save(char(fullfile(strcat(workingDirectory,'Simulation\', subjectName,'-', prefix2save, condition,'-', Model2Use,'-Results-All-Trials.mat'))),'subResults');

save(char(fullfile(strcat(outputFolderGroupData, subjectName2save,'-', prefix2save, condition,'-', Model2Use,'-Results-All-Trials.mat'))),'subResults', '-v7.3');

% Now close all figures
close all;
end

