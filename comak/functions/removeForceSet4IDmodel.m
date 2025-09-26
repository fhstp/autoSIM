function [path2only4IDmodel_noForceSet] = removeForceSet4IDmodel(path2scaledModel, path2bin, side)
% This file removes the entire force set from the model an save a new model copy without the forceset.
% The reason is, that Bryce suggested to rerun only the ID without the
% force set, as otherwise the ID results get biased.


%% Prepare variables
side = lower(side(1)); % take only first char of side and make lower case

%% Check if only4ID model exists or create a copy of the scaled model, change name, and create output path
[filepath,name,ext] = fileparts(path2scaledModel);
path2only4IDmodel_noForceSet = strcat(filepath,'\',name,'-only4ID_noForceSet', ext); % <-- This is the output of this function

if exist(path2only4IDmodel_noForceSet, 'file')
    disp('>>>>> Model without forceset for ID found. ****');
else
    copyfile(path2scaledModel,path2only4IDmodel_noForceSet);

    %% Adapt the model
    import org.opensim.modeling.*
    opensimCommon.LoadOpenSimLibrary(strcat(path2bin,'jam_plugin'));
    model = Model(path2scaledModel);
    
    % Create empty forceset
    efs = ForceSet();

    % set empty force set
    model.set_ForceSet(efs);

    % Update model - write forceset back
    model.upd_ForceSet();

    % Write new model
    model.print(path2only4IDmodel_noForceSet);

    % Display status
    disp(char(strcat('>>>>>',{' '}, strcat(upper(side(1)),lower(side(2:end))), {' '}, 'Model without forceset created for Inverse Dynamics!')));
end

%% Clear variables except output to prevet memory leak.
clearvars -except path2only4IDmodel_noForceSet
end
