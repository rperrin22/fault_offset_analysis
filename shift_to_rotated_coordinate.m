function shifted_coordinate = shift_to_rotated_coordinate(x, y, angle) %#ok<DEFNU>
% Convert the angle to radians
theta = deg2rad(angle);

% Define the rotation matrix
rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Shift the coordinate to the new rotated coordinate system
rotated_coordinate = rotation_matrix .* [x; y];

% Return the shifted coordinate
shifted_coordinate = rotated_coordinate';
end