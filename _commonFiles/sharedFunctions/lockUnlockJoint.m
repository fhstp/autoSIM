function lockUnlockJoint(modelFile, apiPath, joint, lockOrUnlock, dofs)
% This function uses the API to lock/unlock a specified joint.

% modelFile = full file path to model
% joint = name of the joint to lock, e.g. 'subtalar_r'
% lockOrUnlock = 'lock' or 'unlock'
% dofs to lock = 0:2 rotation (xyz), 3:5 trans, e.g. [0:2] or [3:5]


% Load the API
import org.opensim.modeling.* 
if ~isempty(apiPath)
	opensimCommon.LoadOpenSimLibrary(apiPath); % this is a comak specific plugin
end 

% Load model
model = Model(modelFile);

% change the dofs for the joint
Joint = model.get_JointSet().get(joint);
for i = 1 : length(dofs)
    dof2Change = Joint.upd_coordinates(dofs(i));
    switch lockOrUnlock
        case 'lock'
            dof2Change.set_locked(true);
        case 'unlock'
            dof2Change.set_locked(false);
    end      
end

% Save model    
model.print(modelFile);

%% Clear variables except output to prevet memory leak.
clearvars
end
