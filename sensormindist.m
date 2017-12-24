function [min_dist, p_sensor, p_voxel] = sensormindist(sensor, focus, obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

min_dist = [];
p_sensor = [];
p_voxel = [];

%%% RAYS DISTANCE %%%
[idx_list, rays_distance, rays_intercept, ~] = sensorsbeams(sensor, focus, obj);

% Check if it is empty to avoid error
if (~isempty(rays_intercept))
    
    %%% MIM DISTANCE %%%
    min_ray = min(rays_distance);
    min_idx_list = find(rays_distance == min_ray);
    min_idx = min_idx_list(round(length(min_idx_list)/2));
    
    % Ray tracer
    p1 = sensor(idx_list,:);    % P sensor
    p2 = rays_intercept;        % P voxel
    plot3([p1(:,1),p2(:,1)]',[p1(:,2),p2(:,2)]',[p1(:,3),p2(:,3)]');
    
    % Outputs
    min_dist = rays_distance(min_idx);
    p_sensor = p1(min_idx);
    p_voxel = p2(min_idx);
    
end

end

