import org.opensim.modeling.*;

%% Importing data files to run operations
tStart=cputime; % start of timer for the analysis
[BaseModelFile,modelFolder]=uigetfile('*.osim','Pick the base model file to use');% choose the baseline unscaled generic model
model = Model(fullfile(modelFolder, BaseModelFile)); %starting unscaled model
modelFile = regexprep(BaseModelFile, '.osim', '_itr.osim'); % changing name of the copy of the starting model
model.print(fullfile(modelFolder,modelFile));% duplicate the generic model for running the iterations 
for u=0:model.getJointSet.getSize-1 % for every joint
    for v=0:model.getJointSet.get(u).numCoordinates-1 % get every coordinate from every joint
       model.getJointSet.get(u).get_coordinates(v).set_locked(model.getJointSet.get(u).get_coordinates(v).get_locked())% check every coordinate
    end
end
[TRCFile,TRCFolder]=uigetfile('*.trc','Pick Static Trial file to use');% choose the TRC file to use for the scaling
[StatTRC,HeadTRC,HeadTRC_XYZ] = load_trc(fullfile(TRCFolder, TRCFile));%loading the TRC file
[SetupFile,startingSetupFolder]=uigetfile('*.xml','Pick the Scaling setup file to use');%loading the setup file to use
% Simres_folder = inputdlg({'Folder name'}, 'Name of the results folder');
% if ~exist(fullfile(modelFolder, Simres_folder{:}), 'dir')
%     mkdir(fullfile(modelFolder, Simres_folder{:}))
% end
% SetupFolder = inputdlg({'Folder name'}, 'Name of the setup files folder');
% if ~exist(fullfile(modelFolder, SetupFolder{:}), 'dir')
%     mkdir(fullfile(modelFolder, SetupFolder{:}))
% end

%% parameters to set
% Soggetto GC5: 180cm e 75Kg;   Soggetto GC6: 170cm e 70Kg
Info=inputdlg({'Subject height (cm):', 'Subject weight (kg):','Model height (cm):', 'Model weight (kg):'},...
    'Input data',[1 50; 1 50;1 50;1 50], {'180', '75','180', '75'});
SubjectHeight = str2double(Info{1}); %height of subject (cm)
SubjectWeight = str2double(Info{2}); %weight of subject (Kg) 
GenericModelHeight=str2double(Info{3});% Height of generic model (cm): Rajagopal=170 ; TL=175 ; gait2392=180;
GenericModelWeight=str2double(Info{4});%weight of generic model (Kg)
%input for pose estimation
poseAnsw = questdlg('Do you want to estimate the subject pose?');
if strcmp(poseAnsw,'Yes')
 pose=1; % do you want the generic model try to match the experiemntal pose? (default:yes)
else pose=0;
end
 %% defining the manual scale factors
MeanScaleFactY=SubjectHeight/GenericModelHeight; % mean manual scale factor = heights ratio
MeanScaleFact = Vec3(MeanScaleFactY); % from vector to Vec3

%% Display useful info
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
