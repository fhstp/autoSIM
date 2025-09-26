%% Example start script to call the ASTool.

%% Input information.
path2BaseModel = 'C:\dumpFolder\implementationTest\AstData\Rajagopal2015_Unscaled.osim'; 
path2trc = 'C:\dumpFolder\implementationTest\AstData\Static_for_scaling.trc'; 
path2SetUpFile = 'C:\dumpFolder\implementationTest\AstData\Setup_Scale.xml'; 
workingDir4AST = 'C:\dumpFolder\implementationTest\AstData\v2'; 

% Input data for model and subject mass and height for initial scaling.
modelData.SubjectHeight = 180;        % Height of subject (cm)
modelData.SubjectWeight = 75;         % Weight of subject (Kg)
modelData.GenericModelHeight= 170;         % Height of generic model (cm): Lenhart2015 ohne arme: 1560 Rajagopal=170 ; TL=175 ; gait2392=180;
modelData.GenericModelWeight= 75;     % Weight of generic model (Kg): Lenhart2015 ohne arme: 55.9486

% ---

%path2BaseModel = 'C:\dumpFolder\dat\AST_MH_K _R_04_adapted.osim';
%path2trc = 'C:\dumpFolder\dat\AST_04_neutral_marker.trc';
%path2SetUpFile = 'C:\dumpFolder\dat\AST_Scale_Setup_P04.xml';
%workingDir4AST = 'C:\dumpFolder\dat\ASTresults2';

% ---

%path2BaseModel = 'C:\dumpFolder\Example_FullBodyModel_AST\Example_FullBodyModel_AST\Rajagopal2015_Unscaled.osim';
%path2trc = 'C:\dumpFolder\Example_FullBodyModel_AST\Example_FullBodyModel_AST\Static_for_scaling.trc';
%path2SetUpFile = 'C:\dumpFolder\Example_FullBodyModel_AST\Example_FullBodyModel_AST\Setup_Scale.xml';
%workingDir4AST = 'C:\dumpFolder\Example_FullBodyModel_AST\Example_FullBodyModel_AST\ASTresults2';

% Input data for model and subject mass and height for initial scaling.
%modelData.SubjectHeight = 180;        % Height of subject (cm)
%modelData.SubjectWeight = 75;         % Weight of subject (Kg)
%modelData.GenericModelHeight= 170;    % Height of generic model (cm): Lenhart2015 ohne armer: 1560 Rajagopal=170 ; TL=175 ; gait2392=180;
%modelData.GenericModelWeight= 75;     % Weight of generic model (Kg): Lenhart2015 ohne armer: 55.9486

% API path 2 osimJAM
apiPath  = ''; % osimJAMPlugin_osim_v4.2.dll;

% User selection.
pose = 1;

trcScale = 1000;

% Call the ASTool.
path2model = ASTool(path2BaseModel, path2trc, path2SetUpFile, workingDir4AST, modelData, pose, apiPath, trcScale);

% Notes:
% I get differetn rsults depending on if I defined the static trial in the setup gile or not <Static.trc> ...CHECK
% Save Model with ATTACHED markers before use!
% The initial scalesetptup shoudl have no models defined ... set all to Unassigned

% (solved) OpenSim 4.3 works only with my models.
% (solved - changed code to work in both cases) it seems as if in the scalesetup file only setup should be included with true - if tehre some with false they are still used! thsi raises an error.
% Need to use newer  *.dll to make models work with OpenSim 4.3, old modle seems to work - I do not have the left copy of the new one! 
% Del arms from lenahert + but add mass to torso! use new dll and also update bin!

