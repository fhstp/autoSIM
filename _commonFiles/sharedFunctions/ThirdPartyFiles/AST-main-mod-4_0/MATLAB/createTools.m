import org.opensim.modeling.*;

%% Scale tool: loading pre-existent scaling setup file then uploading it
Scaler=ScaleTool(fullfile(startingSetupFolder,SetupFile)); % opening the Scale tool
Array2=ArrayStr(); % defining array object
Array2.append("measurements manualScale ");% initial setup has both manual and measurement scale
Scaler.getModelScaler.setScalingOrder(Array2);% in the scaling order we want the array just created ("manualScale")
ScaledFileName='ModelScaled_API.osim';%name of the Scaled model updated at every iteration
Scaler.getModelScaler.setOutputModelFileName(ScaledFileName);%setting the name of the new scaled model
Scaler.getGenericModelMaker.setModelFileName(modelFile); % set the generic model name in the scale tool
Scaler.getMarkerPlacer.setMarkerFileName(TRCFile);%loading the TRC file with exp markers to use to perform the scaling factors
time_range=Scaler.getModelScaler.getTimeRange; % get the time range frome initial scale setup file
Scaler.getMarkerPlacer.setTimeRange(time_range); % set the time range
Scaler.getMarkerPlacer.setApply(0); % make sure Markers won't be repositioned after scaling
Scaler.print(fullfile(startingSetupFolder,SetupFile)); % saving setup file

%% ========================= Start changes ================================

%% Creation Ik tool for Static trial from Scaling setup
IKSet=Scaler.getMarkerPlacer().getIKTaskSet();% getting IK Sets from Scaling setup
ikTool = InverseKinematicsTool();% define the IK setup tool
%ikTool.set_IKTaskSet(IKSet);%set Ik set previously retrieved ==> bhorsak commented
ikTool.setMarkerDataFileName(TRCFile);% set the trial data .TRC 
%ikTool.set_report_marker_locations(1);% true on "Report_marker_location"
CoordFileName=('Coord_Static.mot'); % name of Coordinates file
ikTool.setOutputMotionFileName(fullfile(modelFolder,CoordFileName));%setting path of output motion file
%TRCData=MarkerData(fullfile(TRCFolder,TRCFile)); % MarkerData object from TRC file
ikTool.setStartTime(time_range.get(0));% getting initial time of scaling trial
ikTool.setEndTime(time_range.get(1));% getting end time of scaling trial
path_ik_static=fullfile(modelFolder,'IkSetup(static_trial).xml');% path of ik setup file for scaling trial
ikTool.print(path_ik_static);%saving ik setup file for static trial

% Dirty hack #1 since not supported for OpenSim API 4.0
changeXML(path_ik_static,'report_marker_locations','true',1);
changeXML(path_ik_static,'report_errors','false',1);

% EXTRA dirty hack #2 since setting IKTaskSet is not supported in OpenSim API 4.0

% Write temp task set to file.
path2IKtaskSet = fullfile(modelFolder,'IKtaskTmp.xml');
IKSet.print(path2IKtaskSet);

% Open target ik static setup file
fid  = fopen(path_ik_static,'r');
f = fread(fid,'*char')';
fclose(fid);

% Find section for IKTaskSet
idxStart = strfind(f,'<IKTaskSet>');
idxEnd = strfind(f,'</IKTaskSet>');

% Find section to replace.
old = f(idxStart:idxEnd);

% Find string to replace old section with.
fid_new  = fopen(path2IKtaskSet,'r');
fnew = fread(fid_new,'*char')';
fclose(fid_new);
idxStart = strfind(fnew,'<IKTaskSet>');
idxEnd = strfind(fnew,'</IKTaskSet>');
new = fnew(idxStart:idxEnd);

% replace old section with new section.
fchanged = strrep(f,old,new);

% Save the model.
fid  = fopen(path_ik_static,'w');
fprintf(fid,'%s',fchanged);
fclose(fid);

%% ========================== END changes =================================

%% Creation of Scale tool with Manual scale factor if RMS erorr exceeds ManualScaleErr Threshold
ScalerManual=Scaler; % new Scale tool 
ScalerManual.setSubjectMass(SubjectWeight);%set Subject mass
ScalerManual.setSubjectHeight(SubjectHeight);% set subject height
NumBodies=model.getBodySet.getSize; % retrieving Number of bodies of the model
% inserting manual scale factors for each body of the subject
for m=0:NumBodies-1 % OpenSim starts from 0 not from 1
    scale=Scale(); % defining scale object
    scale.setScaleFactors(MeanScaleFact)% set the same scale factor for each body
    scale.setSegmentName(model.getBodySet.get(m).getName); % set body name
    scale.setApply(1);% apply: true
ScalerManual.getModelScaler.getScaleSet.cloneAndAppend(scale);% appending the scale factor on the scale set
end
%
path_scaledFile=fullfile(modelFolder,ScaledFileName);%path of scaled model
ScalerManual.getModelScaler.setOutputModelFileName(ScaledFileName);% setting the name of the new scaled model
ScalerManual.getModelScaler.setPreserveMassDist(1); % preserve mass: true
ScalerManual.getModelScaler.setMarkerFileName(TRCFile);%setting the TRC file name
ScalerManual.getMarkerPlacer.setApply(0);%  Markers won't be repositioned after scaling
ScalerManual.getGenericModelMaker.setModelFileName(modelFile); % set the generic model name in the scale tool
Array=ArrayStr(); % defining array object
Array.append("manualScale");% write inside array object
ScalerManual.getModelScaler.setScalingOrder(Array);% in the scaling order we want the array just created ("manualScale")
ScalerManual.getModelScaler.setTimeRange(time_range);%set the time range
path_manualScale=fullfile(modelFolder, 'ManualScaleSetup.xml');%define the manual scale setup file path 
ScalerManual.print(path_manualScale);  %save manual scale setup file
