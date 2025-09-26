function [tf_angle_used] = getTFangleFromDirKinemStatic(static, varNameKneeAngle_c3d_posFront, varNameKneeAngle_c3d)


% Get TF angle directly from direct kinematic model output from the c3d file.
acq = btkReadAcquisition(static);
angles = btkGetAngles(acq);
% kneeVariable = strcat(upper(side), varNameKneeAngle_c3d);
kneeAdd = median(angles.(varNameKneeAngle_c3d)(:,varNameKneeAngle_c3d_posFront));

tf_angle_used = round(kneeAdd,1);

% Close acq
btkCloseAcquisition(acq);

%% Clear variables except output to prevent memory leak.
clearvars -except tf_angle_used
end