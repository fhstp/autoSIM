function TT = estimateTibiaTorsion(LT, MT, MMAL, LMAL, TOE)
% This function calculates the angle between the malleoli markers and the
% femoral condyle markers to estimate the tibia torsion. It works for both
% left and right body sides. Internal tibia torsion (TT) results in negative
% and external in positive values. Zero TT is assumed when the vectors MT-LT and
% MMAL-LMAL are parralell. The script assumes xyz coordinates where the xy
% plane represents the transverse plane. Angles are calculated as projected
% angles to the xy plane. Note: all makers need to be in the same coordinate
% system.

% Input, coordinates in xyz:
% LT = lateral tibia
% MT = medial tibia
% MMAL = medial malleolus
% LMAL = lateral malleolus
% TOE = toe marker, used to determine internal/external TT

% Written by: Brian Horsak, last edited 09/2024

% Input example Left Tibia xyz:
% LT = [-0.0853849 -0.104111 0.409975];
% MT = [-0.0853849 -0.0536888 0.409975];
% MMAL = [-0.142224 -0.0736061 0.0605097];
% LMAL = [-0.0966324 -0.100558 0.0535023];
% TOE = [0.11 -0.0833836 0.0302173]
% 
% LT = markers.(tib_torsion_Markers_Left{1})(1,:);
% MT = markers.(tib_torsion_Markers_Left{2})(1,:);
% MMAL = markers.(tib_torsion_Markers_Left{4})(1,:);
% LMAL = markers.(tib_torsion_Markers_Left{3})(1,:);
% TOE = markers.RTOE(1,:)
% -------------------------------------------------------------------------

% Calculate the vectors
vector1 = MT - LT;
vector2 = MMAL - LMAL;

% Project the vectors onto the xz-plane
vector1_xz = [vector1(1), vector1(2), 0];
vector2_xz = [vector2(1), vector1(2), 0];

% Calculate the dot product and magnitudes of the projected vectors
dot_product = dot(vector1_xz, vector2_xz);
magnitude1 = norm(vector1_xz);
magnitude2 = norm(vector2_xz);

% Calculate the angle in radians
angle_rad = acos(dot_product / (magnitude1 * magnitude2));

% Convert the angle to degrees
TT = rad2deg(angle_rad);

% Ensure the angle is the inner angle
if TT > 180
    TT = 360 - TT;
end

% Determine if vector2 is pointing towards TOE
vector_toe = TOE - LMAL;
vector_toe_xz = [vector_toe(1), 0, vector_toe(3)];
dot_product_toe = dot(vector2_xz, vector_toe_xz);

if dot_product_toe < 0
    % Vector2 is pointing in the opposite direction of TOE
    TT = -TT;
end

end


