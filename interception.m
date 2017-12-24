function [ p_out ] = interception(line, plane)
%INTERCEPTION Summary of this function goes here
%   Detailed explanation goes here

p_a = line(1,:);
p_b = line(2,:);

p_dif = p_a - p_b;

p_0 = plane(1,:);
p_1 = plane(2,:);
p_2 = plane(3,:);

p_out = [];

M = [(p_dif)', (p_1 - p_0)', (p_2 - p_0)'];
% Check if cross the plane
if abs(det(M)) > 10^-10 % det(M) ~= 0
    t_u_v = M^-1 * (p_a - p_0)';
    % Calculate point of interception
    p_cross = p_a - p_dif * t_u_v(1);
    % Verify if point inside the voxel triangle 
    if (pointintriangle(p_cross, p_0, p_1, p_2))
        p_out = p_cross;
    end 
    
else
    disp('')
end

end

function [valid] = sameside(p1,p2,a,b)

    cp1 = cross(b-a, p1-a);
    cp2 = cross(b-a, p2-a);
    
    if (dot(cp1, cp2) >= 0)
        valid = true;
    else
        valid = false;
    end
end

function [inside] = pointintriangle(p, a, b, c)

    if (sameside(p,a,b,c) && sameside(p,b,a,c) && sameside(p,c,a,b))
        inside = true;
    else        
        inside = false;
    end
end






















