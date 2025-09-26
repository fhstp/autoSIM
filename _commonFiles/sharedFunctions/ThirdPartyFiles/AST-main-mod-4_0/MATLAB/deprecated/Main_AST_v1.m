function Main_AST_v1(path2BaseModel, path2trc, path2SetUpFile, workingDir4AST, modelData, pose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AST: Tool for Automatic Scaling of generic MSK OpenSim models                
% Authors: Andrea Di Pietro, University of Pisa, (Italy)                       
%          Alex Bersani, Alma Mater Studiorum - University of Bologna, (Italy) 
%          Uploaded on 23/12/2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Input:
%       - Unscaled OpenSim model
%       - TRC Static File
%       - Starting Scale Setup File
%
%   Output:
%       - Scaled MSK model with Marker Registering and all unlocked coordinates
%       - Scaled MSK model with Marker Registering
%       - Solving Setup Files
% 
%
% Created for OpenSim 4.3
%
% "AST: an OpenSim-based tool for the automatic scaling of generic musculoskeletal models" © 2024
% by Andrea Di Pietro and Alex Bersani is licensed under CC BY-NC 4.0 
% Please cite : "AST: an OpenSim-based tool for the automatic scaling of generic musculoskeletal models; 
%                Andrea Di Pietro, Alex Bersani, Cristina Curreli,Francesca Di Puccio; Computers in Biology and Medicine; 2024"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edited by Brian Horsak {brian.horsak@fhstp.ac.at}:

% The script can now be used as a function that accepts the paths to the *.trc static trial, *.osim base model, 
% and *.xml base scale setup files, and model infos as input. It then copies these files to
% a user specified working directory. If the folder does not exists, the srcipt will create it.
% Only these copies will be changed. All results will be placed within this working directory. 

% Input: 
% ---
% path2BaseModel = 'C:\dumpFolder\dat\AST_MH_K _R_04_adapted.osim';
% path2trc = 'C:\dumpFolder\dat\AST_04_neutral_marker.trc';
% path2SetUpFile = 'C:\dumpFolder\dat\AST_Scale_Setup_P04.xml';
% workingDir4AST = 'C:\dumpFolder\dat\ASTresults';
% modelData = struct() as

%       modelData.SubjectHeight = 182;        % Height of subject (cm)
%       modelData.SubjectWeight = 86;         % Weight of subject (Kg)
%       modelData.GenericModelHeight= 170;    % Height of generic model (cm)
%       modelData.GenericModelWeight= 75;     % Weight of generic model (Kg)

% pose = 1; default = 1; Match the experimental pose? Yes = 1, No =0
% ----

% Change log:
% #1 replaced the "inputAST.m with the following section between "Start Changes" and "No Further Changes"
% #2 load_trc.m: added one line to ignore empty lines in trc files.
% #3 AST_core_v1.m: added option to use trc files in m or mm.

%% ============================ Start Changes =============================

%% User Settings.

% Are the *.trc values in m oder mm?
trcScale = 1000; % 1000 or 1; Set to 1000 if values are in mm and set to 1 if values are in m.

% Set name of the finally scaled model.
name_ModelScaledAdj = 'ModelScaledMarkerAdj.osim'; 

%% Copy all three input files to AST-working dir (defined by user) and set new paths.
sourceFiles = {path2BaseModel, path2trc, path2SetUpFile};

if ~exist(workingDir4AST, 'dir')
    mkdir(workingDir4AST);
end

for i = 1:length(sourceFiles)
    copyfile(sourceFiles{i}, workingDir4AST);
end

%% Import OpenSim API.
import org.opensim.modeling.*;
% pluginPath='C:\OpenSim 4.3\plugins';
% opensimCommon.LoadOpenSimLibraryExact(fullfile(pluginPath,'osimJAMPlugin_osim_v4.2.dll')); % https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53089879/Using+Plugins

%% Importing data files to run operations, parameters to set, defining the manual scale factors

% Set timer
tic
tStart = cputime; % start timer for the analysis

% Get model and subject data.
SubjectHeight = modelData.SubjectHeight; % Height of subject (cm)
SubjectWeight = modelData.SubjectWeight; % Weight of subject (Kg)
GenericModelHeight = modelData.GenericModelHeight; % Height of generic model (cm): Rajagopal=170 ; TL=175 ; gait2392=180;
GenericModelWeight = modelData.GenericModelWeight; % Weight of generic model (Kg)

% Create Setup File path and Filename
[~, file, ext] = fileparts(path2SetUpFile);
startingSetupFolder = workingDir4AST;
SetupFile = [file, ext];

% Read base model *.osim file.
[~, file, ext] = fileparts(path2BaseModel);
modelFolder = workingDir4AST;
BaseModelFile = [file, ext];

model = Model(fullfile(modelFolder, BaseModelFile)); % Starting unscaled model
modelFile = regexprep(BaseModelFile, '.osim', '_itr.osim'); % Change name of the copy of the starting model.
model.print(fullfile(modelFolder,modelFile)); % Duplicate the generic model for running the iterations. 

for u=0:model.getJointSet.getSize-1 % For every joint.
    for v=0:model.getJointSet.get(u).numCoordinates-1 % Get every coordinate from every joint.
       model.getJointSet.get(u).get_coordinates(v).set_locked(model.getJointSet.get(u).get_coordinates(v).get_locked()) % Check each coordinate.
    end
end

% Read static *.trc file.
[~, file, ext] = fileparts(path2trc);
TRCFolder = workingDir4AST;
TRCFile = [file, ext];
[StatTRC,HeadTRC,HeadTRC_XYZ] = load_trc(fullfile(TRCFolder, TRCFile)); % Load the TRC file.

% Calculate the manual scale factor.
MeanScaleFactY = SubjectHeight/GenericModelHeight; % Calculate mean manual scale factor (height ratio).
MeanScaleFact = Vec3(MeanScaleFactY); % From vector to Vec3.

% Display status info.
disp('------------------------------')
disp(['Model unscaled : ', num2str(regexprep(BaseModelFile, '.osim', ''))])
disp('------------------------------')
disp(['Generic model height : ', num2str(GenericModelHeight), ' cm'])
disp(['Generic model weight : ', num2str(GenericModelWeight), ' kg'])
disp(['Subject height : ', num2str(SubjectHeight), ' cm'])
disp(['Subject weight : ', num2str(SubjectWeight), ' kg'])
disp('------------------------------')
disp(['Average scale factor (ASF): ', num2str(MeanScaleFactY)])
disp('------------------------------')

%% ========================== No Further Changes ==========================

%% Setup Tool user paramenters
Km=4; % iterations Threshold : number of iterations 
EndErr=0.004; % End condition for RMS error of loop while
ManualScaleErr=0.025; % error over which the scaling becomes manual for all bodies 
rep=4; % number of times to perform manual scaling just for detected segment 
rep2=8;%  number of times to perform manual scaling for all segments
%% start of algorithm's parameter : Don't modify these parameters
s=-1; % initially the sign for adding the position increment is negative
k=1; % first cycle
flag=0; % flag = 1 if the tool is scaling with Manual scaing factor just for some segments
ind=0; % counter for flag used as a controller
flag2=0;% flag2 =  if the tool is scaling with Manual scaing factor just for all segments
ind2=0;% counter for flag2 used as a controller

%% Creation of algorithm required tool
% Scale tool: loading pre-existent scaling setup file then uploading it
% Creation Ik tool for Static trial from Scaling setup
% Creation of Scale tool with Manual scale factor if RMS erorr exceeds ManualScaleErr Threshold
run createTools

%% determininig coordinate values
if pose==1 % pose =1 means you have chosen to match the experimental pose
    ScaleTool(path_manualScale).run; %Run the manual scaling tool
    ScaledModelFirst= Model(fullfile(modelFolder,ScaledFileName)); % calling the Scaled model
    ikCoord=InverseKinematicsTool(path_ik_static); % call back IK tool for static trial
    ikCoord.setModel(ScaledModelFirst); % set the Scaled model in the IK tool
    ikCoord.run; % run IK
    [CoordData, Coordhead]=load_mot(fullfile(modelFolder,CoordFileName));% load the just computed coordinates
    CoordData=CoordData(:,2:end); % exclude time column from the IK result file
    AvgCoordData = deg2rad(mean(CoordData));% averaging coordinates over time and convert to radians
    AvgCoordData(4:6)=zeros(1,3); % The translation of the pelvis are set to 0 !!!!! To modify in case of different coordinates sequence of the model !!!!!
    %putting the coordinate values inside scaling marker placer
    d=1;
    for u=0:model.getJointSet.getSize-1 % for every joint
        for v=0:model.getJointSet.get(u).numCoordinates-1 % get every coordinate from every joint
            model.getJointSet.get(u).get_coordinates(v).set_clamped(0); % not clamped
            model.getJointSet.get(u).get_coordinates(v).set_default_value(AvgCoordData(d));% insert in every coordinate the relative computed value from IK 
            d=d+1;
        end
    end
    model.print(fullfile(modelFolder,modelFile));%save the model with new coordinates
end
markerset=model.getMarkerSet; % getting the generic markerset from unscaled model
markerset.print(fullfile(modelFolder,'MarkerSet.xml')); % printing the markerset
Nmarkers=markerset.getSize;% retrive number of markers

%% Execution of Autoscaling
run AST_core_v1

%% create the final scaled model with marker placement
if ManualBodies~=0
    AdjScaler=ScaleTool(path_ScalerMix);
else 
    AdjScaler=ScaleTool(fullfile(modelFolder,SetupFile));%New scale tool with marker adjustments
end
AdjScaler.getGenericModelMaker.setModelFileName(modelFile);
AdjScaler.getMarkerPlacer.setApply(1);%repositioning markers after scaling
AdjScaler.getMarkerPlacer.setOutputModelFileName(name_ModelScaledAdj);%set the scaled model name
path_SetupScaleAdj=fullfile(modelFolder, 'ScalingSetupMarkerAdj.xml');
AdjScaler.print(path_SetupScaleAdj); %save setup file
ScaleTool(path_SetupScaleAdj).run; %create model
modelScaledUnlocked=UnlockModel(modelFolder,name_ModelScaledAdj);% unblocking coordinates to scaled model if at least one locked coordinate has been detected
tEnd=cputime;
ElapsedTime=tEnd-tStart;

tempo_exc = toc;
tempo_minuti_exc = tempo_exc/60;

end