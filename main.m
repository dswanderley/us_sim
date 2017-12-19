close all; clear; clc;

obj_str = 'torus-188';

if  exist([obj_str ,'.mat'],'file')
    x_min = min(obj.v(:,1));
    x_offset = 5;
    obj.v(:,1) = obj.v(:,1) - x_min + x_offset;
    obj = load('torus.mat');
else
    obj = readObj([obj_str, '.obj']);
    save('obj.mat','-struct', 'obj');
end

scale = 0.05;
offsets = [1, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);

figure,
trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','red','FaceAlpha',0.1)
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
% Hold Trisurf

%     z
%      |__ y
%   x /

r_max = 2*10^-3; N_r = 10; N_ang = 36;
center = [0, 0.5, 0.75];
% sensor = circlesensor(r_max, N_r, N_ang, center);
sensor_shape = [10, 4];
sensor = rectsensor(sensor_shape, sensor_shape/10,  center);
% sensor = sortrows(sensor);

figure,
plot3(sensor(:,1),sensor(:,2),sensor(:,3),'.')
xlabel('x'); ylabel('y'); zlabel('z');

%%% RAYS DISTANCE %%%
n_rays = size(sensor,1);
rays_distance{n_rays} = [];
rays_intercept{n_rays} = [];
rays_voxel_id{n_rays} = [];
for k = 1:n_rays
    
    p_sensor = sensor(k,:);             % Transducer point
    p_focal = p_sensor - [ 1, 0, 0];    % Focal point
    line = [p_sensor; p_focal];
    
    [ d_min, p_int,  face_int_idx ] = minvoxeldist(obj, line);
    rays_distance{k} = d_min;
    rays_intercept{k} = p_int;
    rays_voxel_id{k} = face_int_idx;
end

rays_distance{1} = [];
rays_distance{10} = [];
idx_list = find(~cellfun(@isempty,rays_distance));

rays_distance = cell2mat(rays_distance(:));
rays_intercept = cell2mat(rays_intercept(:));
rays_voxel_id = cell2mat(rays_voxel_id(:));

figure,
stem3(sensor(idx_list,2), sensor(idx_list,3), rays_distance)
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

%[idx_l, idx_c] = ind2sub(sensor_shape, idx_list);
im_sensor = zeros(sensor_shape);
im_sensor(idx_list) = abs(rays_distance/max(rays_distance)-1)+0.5;
im_sensor = rot90(rot90(im_sensor));
figure, imshow(im_sensor)

%%% MIM DISTANCE %%%
min_ray = min(rays_distance);
min_idx = find(rays_distance == min_ray);
min_idx = min_idx(1);

disp(['dist min: ', num2str(min_ray)])
disp(['interception: ', num2str(rays_intercept(min_idx))])


