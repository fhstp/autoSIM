function [fig, T] = reportIKmarkerErrors(stoFilePathIK, ExperimentalMarkersPath, workingDirectory, statesFileName);
% This file plots the distacne and RMSE between the *.trc marker locations
% and the ones from the IK error report file:

%% Read locations from IK and make table
stoFileIK = load_sto_file(stoFilePathIK);
IKlocsTable = struct2table(stoFileIK);

%% Read locations from experimental trial and make table with same headers
expMarkers = read_trcFile(ExperimentalMarkersPath);

% Create header
headerExpMarker = {};
direction =  {'x','y','z'};
cnt = 1;
for j = 1 : length(expMarkers.MarkerList)
    for i = 1 : 3
        str = direction{i};
        headerExpMarker{cnt,1} = strcat(expMarkers.MarkerList{j},'_t',str);
        cnt = cnt + 1;
    end
end

% Sync data
startTime = IKlocsTable.time(1);
startIdx = find(expMarkers.Data(:,2) == startTime);
endTime = IKlocsTable.time(end);
endIdx = find(expMarkers.Data(:,2) == endTime);

% Create array, exclude Frame Number
locationsIK = expMarkers.Data(startIdx:endIdx,3:end);
expMarkerLocsTable = array2table(locationsIK, 'VariableNames', headerExpMarker);
time = expMarkers.Data(startIdx:endIdx,2);

%% Now calculate the distance and plot it
fig = figure;
set(gcf, 'Units', 'normalized', 'Position', [0, 0, .8, .5]);

% Initialize T
T = table();
T.time = time;

% Prepare list of markers
vars2Check = IKlocsTable.Properties.VariableNames; 
vars2Check = vars2Check(~strcmp(vars2Check,'time')); % remove time
vars2Check = vars2Check(contains(vars2Check,'_tx')); % get only the ones with _tx
vars2Check = cellfun(@(x) x(1:end-3), vars2Check, 'un', 0); % remove substring _tx
vars2Check = sort(vars2Check);

% Prepare subplot
ncols = 10;
nrows = ceil(length(vars2Check) / ncols);

% Now do the calc
for i = 1 : length(vars2Check)
    var2Check = vars2Check{i};
    A = [expMarkerLocsTable.(strcat(var2Check,'_tx')), expMarkerLocsTable.(strcat(var2Check,'_ty')), expMarkerLocsTable.(strcat(var2Check,'_tz'))];
    B = [IKlocsTable.(strcat(var2Check,'_tx')), IKlocsTable.(strcat(var2Check,'_ty')), IKlocsTable.(strcat(var2Check,'_tz'))];
    d = sqrt(sum((A' - B') .^ 2)) * 100; % Distance, in cm
    RMSE = sqrt(mean((A' - B').^2)) * 100;  % Root Mean Squared Error, in cm

    subplot(nrows,ncols,i)
    hold on
    plot(d, 'b-');
    plot(RMSE, 'r-');
    if i == 1; yline([4 2],'r--',{'4cm','2cm'}); else; yline([4 2],'r--'); end
    ylim([0 5])
    ylabel(var2Check, 'FontWeight', 'bold')
    if i == 1; legend('distance','RMSE','location', 'northwest'); legend boxoff; end

    % Create a dynamic table here
    T.(var2Check) = d';
    
end
sgtitle({char(strrep(strrep(workingDirectory,{'_'},''),'\','/')),char(strcat(strrep(strrep(statesFileName,'_',' '),'/','\')))});

end