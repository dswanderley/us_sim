close all; clear; clc;
%% OBJ Settings %%
obj_str = 'torus-188';

if  exist([obj_str ,'.mat'],'file')
    obj = load([obj_str, '.mat']);
else
    obj = readObj([obj_str, '.obj']);
    save([obj_str, '.mat'],'-struct', 'obj');
end

scale = 0.5;
offsets = [8, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);
%%
%     z
%      |__ y
%   x /
%
%% Sensors and Array Settings %%
MAX_DEPTH = 30;         % cm
% Rectangular sensor
sensor_size = [3, 3];   % unit (number of beams)
% Angle in rad - 0 for collinear
max_angle = 0 * pi /180;
if (max_angle == 0)
    sensor_size = [1, 1];   % for collinear case
end
% Number of Beams in each sensor
sensor_beams = sensor_size(1) * sensor_size(2);
% Dimension in physical scale of the array
array_dim = [120, 120]/10;      % cm
% Sensors arrangement in array
array_size = [6, 6];            % unit (number of sensors in each direction)
% Number of sensors on array
num_sensors = array_size(1) * array_size(2);

% List of spatial informations settings
safe_guard = 10;
[list_pi, list_u, list_v ] = strides(obj.v, array_dim, safe_guard);

% Max volume
x_min = min(list_pi(:,1))-3*safe_guard;  x_max = max(list_pi(:,1))+3*safe_guard;
y_min = min(list_pi(:,2))-3*safe_guard;  y_max = max(list_pi(:,2))+3*safe_guard;
z_min = min(list_pi(:,3))-3*safe_guard;  z_max = max(list_pi(:,3))+3*safe_guard;

%% Loop walking with sensor
volume{size(list_pi,1)} = [];
s_points{size(list_pi,1)} = [];
tic;
for l = 1:size(list_pi,1)
    
    % Positions of the sensor array
    u = list_u(l,:);
    v = list_v(l,:);
    p0 = list_pi(l,:);    
    orientation = cross(u, v); % normal
    % orientation = list_ori(l,:);
    
    % Get sensors
    [sensors, centers, srect] = sensorarray(sensor_size, array_dim, array_size, u, v, p0);
    
    % GET FOCAL POINTS
    array_center = centers(:,:,1);
    sensor_base = sensors(:,:,1);
    
    % Focus based to be applied in each sensor
    focus_base = getfocalpoints(sensor_base, array_center, u, v, orientation, max_angle);
    % Remove offset from sensor_base
    sensor_base(:,1) = sensor_base(:,1) - array_center(1);
    sensor_base(:,2) = sensor_base(:,2) - array_center(2);
    sensor_base(:,3) = sensor_base(:,3) - array_center(3);
    % Remove offset from focus_base
    focus_base(:,1) = focus_base(:,1) - array_center(1);
    focus_base(:,2) = focus_base(:,2) - array_center(2);
    focus_base(:,3) = focus_base(:,3) - array_center(3);
    % Adjust focus position for each beam
    focus = repmat(focus_base, 1, 1, size(sensors,3));
    focus(:,1,:) = focus(:,1,:) + repmat(centers(:,1,:), size(focus(:,1,:),1), 1, 1);
    focus(:,2,:) = focus(:,2,:) + repmat(centers(:,2,:), size(focus(:,2,:),1), 1, 1);
    focus(:,3,:) = focus(:,3,:) + repmat(centers(:,3,:), size(focus(:,3,:),1), 1, 1);
        
    %%% RAYS DISTANCE %%%
    array_dists{size(sensors,3)} = [];
    spatial_points{size(sensors,3)} = [];
    % array_lines{size(sensors,3)} = [];
    for k = 1:size(sensors,3)
        % Get all point in each sensor, their focus and the sensor center
        sensor = sensors(:,:,k);
        fc = focus(:,:,k);
        c = centers(:,:,k);
        % Calculate minimal distance from center
        [min_dist, p_sensor, p_voxel] = sensormindist(sensor, fc, c, obj);
        array_dists{k} = min_dist;
        line = [c; p_voxel];
        spatial_points{k} = p_voxel;
    end
    
    % Get Positions with at least one interception
    idx_list = find(~cellfun(@isempty,array_dists))';
    rays_distance = cell2mat(array_dists(:));
    [idx_2d_y, idx_2d_x] = ind2sub(array_size, idx_list);
    s_points{l} = cell2mat(spatial_points(:));
    
    % Depth Image
    M_dist = MAX_DEPTH*ones(array_size);
    M_dist(idx_list) = rays_distance;
    % Show image
    im_sensor = M_dist; % rot90(rot90(M_dist));
    % im_sensor = padarray(im_sensor,[1,1],'symmetric','both');   % post
    clims = [0 MAX_DEPTH];
    
    %%% PLANE %%%
    R_SAMPLES = 100;
    range_arr = linspace(0, MAX_DEPTH, R_SAMPLES);
     
end
time = toc;

%% Plot points in 3D space

points = cell2mat(s_points(:));

figure('units','normalized','outerposition',[0 0 1 1])
plot3(points(:,1), points(:,2), points(:,3), '.')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
grid
hold on
% pause(5)
% k = boundary(points);
% % trisurf(k,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)
