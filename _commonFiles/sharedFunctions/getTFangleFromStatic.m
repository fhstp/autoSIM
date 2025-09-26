function [tf_angle_used] = getTFangleFromStatic(static, side, bodymass, BodyHeight)

% Get tf angle from the static trial based on the work from Stief et al.
% 2020 https://doi.org/10.1016/j.gaitpost.2020.04.011.

acq = btkReadAcquisition(static);
angles = btkGetAngles(acq);
kneeVariable = strcat(upper(side),'KneeAngles');
kneeAdd = median(angles.(kneeVariable)(:,2));

% Calculate BMI
BMI = bodymass / (BodyHeight/1000)^2;

% Estimate the tf angle from the static.
if BMI <= 25
    tf_angle_used = 0.618 * kneeAdd + 0.009;
    tf_angle_used = round(tf_angle_used,1);

else
    tf_angle_used = 0.651 * kneeAdd + 0.280 * BMI -6.258;
    tf_angle_used = round(tf_angle_used,1);
end

% Close acq
btkCloseAcquisition(acq);

%% Clear variables except output to prevet memory leak.
clearvars -except tf_angle_used
end