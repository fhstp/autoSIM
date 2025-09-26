function dataOk = checkIfExternalDataIsValid(persInfo, MRICTdate, EXAMINATIONdate, threshold)
% Checks if the personalization data in persInfo can be used
% based on the examination date. The date difference must be below 30 days for TT.
%
% Syntax: TorsionOk = checkTorsionOk(persInfo)
%
% Input:
%   persInfo - A struct with fields 'MRICTdate' and 'EXAMINATIONdate'
%   MRICTdate - string, Name of MRICTdate
%   EXAMINATIONdate - string, Name of EXAMINATIONdate
%   threshold - numeric, in days
%
% Output:
%   TorsionOk - Boolean indicating if the date difference is within the allowed range

    % Normalize date
    if ~isdatetime(persInfo.(MRICTdate))
        try
            tmpMRICTdate = datetime(persInfo.(MRICTdate), 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
        catch
            tmpMRICTdate = NaT; % Set to NaT if conversion fails
        end
    else
        tmpMRICTdate = persInfo.(MRICTdate); % Keep the original datetime if it's already valid
    end

    % Normalize date
    if ~isdatetime(persInfo.(EXAMINATIONdate))
        try
            tmpEXAMINATIONdate = datetime(persInfo.(EXAMINATIONdate), 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
        catch
            tmpEXAMINATIONdate = NaT; % Set to NaT if conversion fails
        end
    else
        tmpEXAMINATIONdate = persInfo.(EXAMINATIONdate); % Keep the original datetime if it's already valid
    end

    % Check if both dates are valid datetime and not NaT
    if ~isnat(tmpMRICTdate) && ~isnat(tmpEXAMINATIONdate)
        % Calculate the difference in days using normalized dates
        diffDays = abs(days(tmpMRICTdate - tmpEXAMINATIONdate));

        % Determine if Torsion is OK
        dataOk = diffDays < threshold;
    else
        % Handle the case where one or both dates are NaT
        dataOk = false;
    end
end