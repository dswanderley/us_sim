function [ p_out ] = interception(line, plane)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p_a = line(1,:);
p_b = line(2,:);

p_dif = p_a - p_b;

p_0 = plane(1,:);
p_1 = plane(2,:);
p_2 = plane(3,:);

p_out = [];

M = [(p_dif)', (p_1 - p_0)', (p_2 - p_0)'];
if det(M) ~= 0
    t_u_v = M^-1 * (p_a - p_0)';
    
    p_cross = p_a - p_dif * t_u_v(1);
    
    if p_cross(1) <= max(plane(:,1)) && p_cross(1) >= min(plane(:,1))
        if p_cross(2) <= max(plane(:,2)) && p_cross(2) >= min(plane(:,2))
            if p_cross(3) <= max(plane(:,3)) && p_cross(3) >= min(plane(:,3))
                p_out = p_cross;
            end
        end
    end
    
end

end