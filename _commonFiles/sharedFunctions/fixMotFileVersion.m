function fixMotFileVersion(filePath)
% This file loads a *.mot file an fixes wrong header name conventions:
% E.g., ground_torque1 (=wrong) vs. ground_torque_1

% Written by: Brian Horsak, 13.032023
%--------------------------------------------------------------------------

fid  = fopen(filePath,'r');
f=fread(fid,'*char')';
fclose(fid);

% Check if the headers are correct or not
if ~contains(f, 'ground_force_')

    % If not, replace strings
    f = strrep(f,'ground_torque','ground_torque_');
    f = strrep(f,'ground_force','ground_force_');

    % And overwrite file
    fid  = fopen(filePath,'w');
    fprintf(fid,'%s',f);
    fclose(fid);

    % Tell user that we changed the *.mot files
    warning off backtrace;
    warning('>>>>> An existing *.mot file with a wrong header was identified and corrected to have the header convention, e.g, <ground_force_1_y> instead of <ground_force1_y>.');
    warning on backtrace;
end

%% Clear variables except output to prevet memory leak.
clearvars
end