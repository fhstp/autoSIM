function [InputData, node] = getEventsFF(mot_data, mot_labels, condition, FPs, k, i, trialCnt, stepCnt, paths, events, delta, InputData)

% Decide if walking trial or SEBT trial
% Walking
% Get the data
dat = mot_data(:,contains(mot_labels,strcat(num2str(i),'_vy')));
time = mot_data(:,contains(mot_labels,'time'));

fp_contact = time((dat>5));
IC = fp_contact(1);
TO = fp_contact(end);
side = char(FPs(i));

% Get events via c3d events and mot file
% IC ipsilateral
[~,loc_idx] = min(abs(events.(char(strcat(side,'_Foot_Strike')))-(IC + delta)));
ICi = events.(char(strcat(side,'_Foot_Strike')))(loc_idx + 1);

% contralateral TO % IC contralateral
if strcmp(side,'Right')
    [~,idx_cTO] = min(abs(events.Left_Foot_Off-(IC + delta)));
    cTO = events.Left_Foot_Off(idx_cTO);
    idx_cIC = find((events.Left_Foot_Strike > IC + delta),1, 'first');
    cIC = events.Left_Foot_Strike(idx_cIC);
end

if strcmp(side,'Left')
    [~,idx_cTO] = min(abs(events.Right_Foot_Off-(IC + delta)));
    cTO = events.Right_Foot_Off(idx_cTO);
    idx_cIC = find((events.Right_Foot_Strike > IC + delta),1, 'first');
    cIC = events.Right_Foot_Strike(idx_cIC);
end

% Write to final results
[~,c3dname,~] = fileparts(paths.c3d{k});
side_out = char(lower(side(1)));
node = (strcat('trial',num2str(trialCnt),'_', c3dname,'_',num2str(stepCnt),'_', side_out));

InputData.(node).name = char(strcat(c3dname,'_',num2str(stepCnt))); %char(c3dname);
InputData.(node).c3dPath = char(paths.c3d(k));
InputData.(node).enfPath = char(paths.enf(k));
InputData.(node).trcPath = char(paths.trc(k));
InputData.(node).motPath = char(paths.mot(k));
InputData.(node).static_trcPath = char(paths.trc_static);
InputData.(node).static_motPath = char(paths.mot_static);
InputData.(node).Side =  lower(char(side));
InputData.(node).IC = IC;
InputData.(node).TO = TO; 
InputData.(node).ICi = ICi-delta; % Vicon start frame delta to mot
InputData.(node).cTO = cTO-delta; % Vicon start frame delta to mot
InputData.(node).cIC = cIC-delta; % Vicon start frame delta to mot
InputData.(node).condition = condition;
InputData.(node).fromc3d = events;

%% Clear variables except output to prevet memory leak.
clearvars -except InputData node
end
