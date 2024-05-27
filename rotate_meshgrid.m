function [rotated_x, rotated_y] = rotate_meshgrid(x, y, angle, origin) %#ok<DEFNU>
% Convert the angle to radians
theta = deg2rad(angle);

% Define the rotation matrix
rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Shift the origin if provided
if nargin == 4
    x = x - origin(1);
    y = y - origin(2);
end

% Flatten x and y matrices
x_flat = x(:);
y_flat = y(:);

% Rotate the coordinates
rotated_coordinates = rotation_matrix * [x_flat'; y_flat'];

% Reshape the rotated coordinates
rotated_x = reshape(rotated_coordinates(1, :), size(x));
rotated_y = reshape(rotated_coordinates(2, :), size(y));

% Shift back to original origin
if nargin == 4
    rotated_x = rotated_x + origin(1);
    rotated_y = rotated_y + origin(2);
end
end