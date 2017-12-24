function [idx_list,rays_distance, rays_intercept, rays_voxel_id] = sensorsbeams(sensor, focus, obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%% RAYS DISTANCE %%%
n_rays = size(sensor,1);
rays_distance{n_rays} = [];
rays_intercept{n_rays} = [];
rays_voxel_id{n_rays} = [];
for k = 1:n_rays
    p_sensor = sensor(k,:);             % Transducer point
    p_focal = focus(k,:);    % Focal point
    line = [p_sensor; p_focal];
    
    [ d_min, p_int,  face_int_idx ] = minvoxeldist(obj, line);
    rays_distance{k} = d_min;
    rays_intercept{k} = p_int;
    rays_voxel_id{k} = face_int_idx;
end

idx_list = find(~cellfun(@isempty,rays_distance));
rays_distance = cell2mat(rays_distance(:));
rays_intercept = cell2mat(rays_intercept(:));
rays_voxel_id = cell2mat(rays_voxel_id(:));

end

