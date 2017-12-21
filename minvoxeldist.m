function [ d_min, p_int,  face_int_idx ] = minvoxeldist( obj, line )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

p_ref = line(1,:);
d_ref = 10^9;
d_min = d_ref;
p_int = [];
face_int_idx = [];

for k = 1:length(obj.f.v)
    
    plane = obj.v(obj.f.v(k,:),:);
    %plot3(plane(:,1),plane(:,2),plane(:,3))
    
    p_out = interception(line, plane);
    
    if ~isempty(p_out)
        
        p_dist = pdist([p_out; p_ref]);
        
        if p_dist < d_min
            d_min = p_dist;
            p_int = p_out;
            face_int_idx = k;
        end
    end   
end

if d_min == d_ref
    d_min = [];
end

end

