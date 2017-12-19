close all; clear; clc;

obj_str = 'torus-188';
% obj_str = 'two-blocks';

if  exist([obj_str ,'.mat'],'file')
    obj = load([obj_str, '.mat']);
else
    obj = readObj([obj_str, '.obj']);
    save([obj_str, '.mat'],'-struct', 'obj');
end

scale = 0.05;
offsets = [1, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);
%     z
%      |__ y
%   x /
r_max = 2*10^-3; N_r = 10; N_ang = 36;
center = [0, 0.5, 0.75];
% sensor = circlesensor(r_max, N_r, N_ang, center);
sensor_shape = [10, 4];
sensor = rectsensor(sensor_shape, sensor_shape/10,  center); % sensor = sortrows(sensor);

figure, % PLOT REAL WORLD
trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','red','FaceAlpha',0.1)
axis equal
hold on
plot3(sensor(:,1),sensor(:,2),sensor(:,3),'.')
xlabel('x'); ylabel('y'); zlabel('z');
%hold off

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

idx_list = find(~cellfun(@isempty,rays_distance));
rays_distance = cell2mat(rays_distance(:));
rays_intercept = cell2mat(rays_intercept(:));
rays_voxel_id = cell2mat(rays_voxel_id(:));

% PLOT RAYS
if (~isempty(rays_intercept))
    p1 = sensor(idx_list,:);
    p2 = rays_intercept;
    plot3([p1(:,1),p2(:,1)]',[p1(:,2),p2(:,2)]',[p1(:,3),p2(:,3)]');
end
hold off

% PLOT DISTANCES
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


%%% PLANE %%% 

MAX_RANGE = 3; % m
R_SAMPLES = 100;
range_arr = linspace(0, MAX_RANGE, R_SAMPLES);

n_sensors = sensor_shape(2);
% 3D Image, where each layer is a plan
V = zeros(R_SAMPLES, sensor_shape(2), sensor_shape(1));
for c = 1:length(idx_list)
    
    idx = idx_list(c);
    [nn, mm] = ind2sub(sensor_shape,idx);
    [val, index] = min(abs(range_arr - rays_distance(c)));
    
    depth = zeros(1, R_SAMPLES); 
    depth(index) = 1;
    
    V(:, mm, nn) = depth;
end

% Slice for 3D vizualization 
[xx, yy, zz] = meshgrid(1:sensor_shape(2) ,1:R_SAMPLES, 1:sensor_shape(1));
xslice = 1:sensor_shape(2); 
yslice = 1:R_SAMPLES;
zslice = 1:sensor_shape(1);
figure,
slice(xx,yy,zz, V,xslice,yslice, zslice)
colormap winter
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

% 3D Volume vizualization
xdata = [0,sensor_shape(2)]; 
ydata = [0,R_SAMPLES];
zdata = [0,sensor_shape(1)];

figure, 
vol3d('cdata', V, 'xdata', xdata, 'ydata', ydata, 'zdata', zdata);
colormap winter
%axis equal off
%set(gcf, 'color', 'w');
%xlabel('x'); ylabel('y'); zlabel('z');




