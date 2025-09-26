function [InputData, node] = getEventsOSS(mot_data, mot_labels, condition, FPs, k, i, trialCnt, stepCnt, paths, events, delta, InputData, markersTrial, cam_rate)

% Decide if walking trial or mSEBT trial
% It is assumed that when condition does not include the substring
% 'mSEBT' it is a walking trial!
if contains(condition, 'mSEBT')
    try
        % mSEBT
        side = char(FPs(i));

        % Write to final results (mSEBT)
        [~,c3dname,~] = fileparts(paths.c3d{k});

        % Get events
        if strcmp(side,'Right'); sideTmpSEBT = 'Left'; end % it is the event from the right foot but relevant for the left one
        if strcmp(side,'Left'); sideTmpSEBT = 'Right'; end

        startLabel = events.([sideTmpSEBT,'_Start']);
        touchDownLabel = events.([sideTmpSEBT,'_Touchdown']);
        stopLabel = events.([sideTmpSEBT,'_Stop']);
        side_out = char(lower(side(1)));
        node = (strcat('trial',num2str(trialCnt),'_', c3dname,'_',num2str(stepCnt),'_', side_out));

        % Get Start Event by marker distance from Start
        if strcmp(sideTmpSEBT, 'Left'); marker2Test = 'LTOE'; end
        if strcmp(sideTmpSEBT, 'Right'); marker2Test = 'RTOE'; end

        idx = 1;
        dist = 0;
        startFrameByLabel = round(startLabel*cam_rate);
        while dist < 100 %10cm
            % Calc the mdn position based on the frames of the second half of
            % frames from first to start Label.
            A = median(markersTrial.(marker2Test)(round(startFrameByLabel/2):startFrameByLabel,:));
            B = markersTrial.(marker2Test)(idx,:);
            dist = norm(A-B);
            idx = idx + 1;
        end

        % Convert frames to time
        startLabel = 1/cam_rate * idx;

    catch
        % Set nan in case an error occures during getting the events
        delta = nan;
        startLabel = nan;
        touchDownLabel = nan;
    end

    InputData.(node).name = char(strcat(c3dname,'_',num2str(stepCnt))); %char(c3dname);
    InputData.(node).c3dPath = char(paths.c3d(k));
    InputData.(node).enfPath = char(paths.enf(k));
    InputData.(node).trcPath = char(paths.trc(k));
    InputData.(node).motPath = char(paths.mot(k));
    InputData.(node).static_trcPath = char(paths.trc_static);
    InputData.(node).static_motPath = char(paths.mot_static);
    InputData.(node).Side =  lower(char(side));
    InputData.(node).IC = startLabel-delta; % Vicon start frame delta to mot
    InputData.(node).TO = 0; % touchDownLabel-delta; % Vicon start frame delta to mot
    InputData.(node).ICi = touchDownLabel-delta - (5*1/cam_rate); %stopLabel-delta; % Vicon start frame delta to mot; Redcue by five frames to make reduced chance we have ground contact
    InputData.(node).cTO = 0; % Vicon start frame delta to mot
    InputData.(node).cIC = 0; % Vicon start frame delta to mot
    InputData.(node).condition = condition;
    InputData.(node).fromc3d = events;

else  % Walking
    try
        % Get the data
        dat = mot_data(:,contains(mot_labels,strcat(num2str(i),'_vy')));
        time = mot_data(:,contains(mot_labels,'time'));

        fp_contact = time((dat>20)); % FP threshold of 20 N.
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

    catch
        % Set nan in case an error occures during getting the events
        delta = nan;
        IC = nan;
        TO = nan;
        ICi = nan;
        cTO = nan;
        cIC = nan;
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
end

%% Clear variables except output to prevet memory leak.
clearvars -except InputData node
end