% Simple tool to measure the midpoint between two 3d coordinates, eg to get
% the midpoint between the RASI and LASI.

format long;

A = [-0.0593732, -0.06931, 0.0715861]
B = [-0.0593732, -0.06931, -0.0829339]

% Parameter t, where 0 <= t <= 1
t = 0.5;  % Example: midpoint

% Calculate the coordinates of the point P between A and B
P = (1 - t) * A + t * B;

% Display the result without scientific notation using sprintf
disp(['Coordinates of point P: (', sprintf('%.15f', P(1)), ' ', sprintf('%.15f', P(2)), ' ', sprintf('%.15f', P(3)), ')']);