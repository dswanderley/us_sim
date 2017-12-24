function [ focus ] = getfocalpoints( sensor, center, u, v, orientation, max_angle )
%GETFOCALPOINTS Summary of this function goes here
%   Detailed explanation goes here

if (max_angle == 0)
    % If collinear each point has a focus
    ori = repmat(orientation, size(sensor,1), 1);
    focus = sensor - ori;
else
    % Angle is half of the max_ange (opening)
    angle = max_angle/2;
    % Calculate distance from the center
    aux = repmat(center, size(sensor,1), 1);
    d = sum((sensor - aux).^2, 2);
    % Get index for the maximun distance
    idx = find(d == max(d), 1);
    p_max = sensor(idx,:);
    % Get focal point
    pf = focalpoint(center, p_max, angle, u, v);
    % Replicate this point for n points in sensor
    focus = repmat(pf, size(sensor,1), 1);
end

end

