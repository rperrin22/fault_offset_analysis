function rotated_vector = rotate_to_vertical(vector) %#ok<DEFNU>
% Extract the coordinates
x1 = vector(1);
y1 = vector(2);
x2 = vector(3);
y2 = vector(4);

% Calculate the angle away from vertical
angle = angle_from_vertical(x1, y1, x2, y2);

% Rotate the line to make it vertical
theta = deg2rad(angle);
rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotated_coordinates = rotation_matrix * [x2 - x1; y2 - y1];

% Update the coordinates
rotated_vector = [x1, y1, x1 + rotated_coordinates(1), y1 + rotated_coordinates(2)];
end