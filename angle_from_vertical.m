function angle = angle_from_vertical(x1, y1, x2, y2)
% Calculate the angle away from the vertical axis
dx = x2 - x1;
dy = y2 - y1;

% Calculate the angle using atan2
angle = atan2(dy, dx);

% Convert radians to degrees
angle = rad2deg(angle);

% Make sure the angle is positive
if angle < 0
    angle = angle + 360;
end

% Subtract 90 to get the angle away from vertical
angle = 90 - angle;
end