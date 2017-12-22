function [ focus ] = getfocalpoints( sensor, center, u, v, orientation, max_angle )
%GETFOCALPOINTS Summary of this function goes here
%   Detailed explanation goes here

if (max_angle == 0)
   
    ori = repmat(orientation, size(sensor,1), 1);
    focus = sensor - ori;
    
else

    aux = repmat(center, size(sensor,1), 1);
    d = sum((sensor - aux).^2, 2);

    idx = find(d == max(d), 1);
    p_max = sensor(idx,:);
    
    pf = focalpoint(center, p_max, max_angle, u, v);
    
    focus = repmat(pf, size(sensor,1), 1);
end


end

