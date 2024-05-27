function [rotated_x, rotated_y] = rotate_meshgrid_with_line(x, y, line_vector) %#ok<DEFNU>
% Extract line coordinates
x1 = line_vector(1);
y1 = line_vector(2);
x2 = line_vector(3);
y2 = line_vector(4);

% Calculate the angle of the line
angle = angle_from_vertical(x1, y1, x2, y2);

% Convert the angle to radians
theta = deg2rad(angle);

% Define the rotation matrix
rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Flatten x and y matrices
x_flat = x(:);
y_flat = y(:);

% Rotate the coordinates
rotated_coordinates = rotation_matrix * [x_flat'; y_flat'];

% Reshape the rotated coordinates
rotated_x = reshape(rotated_coordinates(1, :), size(x));
rotated_y = reshape(rotated_coordinates(2, :), size(y));
end