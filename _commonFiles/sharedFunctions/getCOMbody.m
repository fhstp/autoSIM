function COM = getCOMbody(path2model, segment2get)
% This function gets the COM of a body from an *.osim model.
% I use this in my workflo to i.e. get the COM positio of the pelivs for
% the reserve actuator files.

% segment2get = e.g. 'pelvis'
%--------------------------------------------------------------------------

% Load the API
import org.opensim.modeling.*

% Read model
model = Model(path2model);

%
obj = model.get_BodySet().get(segment2get).get_mass_center();

% Extract values from osim object
for i = 1:3
    COM(i) = obj.get(i-1); % The API indexing starts at 0!
end

%% Clear variables except output to prevet memory leak.
clearvars -except COM
end

