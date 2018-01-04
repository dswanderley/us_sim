close all; clear; clc;

obj_str = 'torus-188';
% obj_str = 'two-blocks';

if  exist([obj_str ,'.mat'],'file')
    obj = load([obj_str, '.mat']);
else
    obj = readObj([obj_str, '.obj']);
    save([obj_str, '.mat'],'-struct', 'obj');
end

scale = 0.5;
offsets = [2, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);
%     z
%      |__ y
%   x /
% Positions of the sensor array
MAX_DEPTH = 30;         % cm
orientation = [1,0,0];  % Array beams orientation
u = [0, 1, 0];          % Array plan vector 1
v = [0, 0, 1];          % Array plan vector 2
% Rectangular sensor
sensor_size = [3, 3];   % unit (number of beams)
% Angle - 0 for collinear
max_angle = 0;  % 30 * pi /180;
if (max_angle == 0)
    sensor_size = [1, 1];   % for collinear case
end
% Number of Beams in each sensor
sensor_beams = sensor_size(1) * sensor_size(2);
% Dimension in physical scale of the array
array_dim = [90, 90]/10;      % cm
% Sensors arrangement in array
array_size = [9, 9];            % unit (number of sensors in each direction)
num_sensors = array_size(1) * array_size(2);    % Number of sensors on array
% Get sensors
[sensors, centers, srect] = sensorarray(sensor_size, array_dim, array_size, u, v, [0,0,0]);

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

%%% PLOT REAL WORLD %%%
figure, 
trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','red','FaceAlpha',0.1)
hold on
% PLOT SENSOR CENTERS
c_x = reshape(centers(:,1,:), 1, num_sensors)'; 
c_y = reshape(centers(:,2,:), 1, num_sensors)';
c_z = reshape(centers(:,3,:), 1, num_sensors)';
plot3(c_x, c_y, c_z, '*')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
% PLOT ARRAY BORDERS
srect = [srect; srect(1,:)];
plot3(srect(:,1),srect(:,2),srect(:,3))

%%% RAYS DISTANCE %%%
array_dists{size(sensors,3)} = [];
array_lines{size(sensors,3)} = [];
for k = 1:size(sensors,3)
    % Get all point in each sensor, their focus and the sensor center
    sensor = sensors(:,:,k);
    fc = focus(:,:,k);
    c = centers(:,:,k);
    % Calculate minimal distance from center
    [min_dist, p_sensor, p_voxel] = sensormindist(sensor, fc, c, obj);
	array_dists{k} = min_dist;
    line = [c; p_voxel]; 
    array_lines{k} = line;
    % PLOT MIN BEAM  
    plot3(line(:,1,:), line(:,2,:), line(:,3,:));
    pause(0.01)
end
hold off

% Get Positions with at least one interception
idx_list = find(~cellfun(@isempty,array_dists))';
rays_distance = cell2mat(array_dists(:));
[idx_2d_y, idx_2d_x] = ind2sub(array_size, idx_list);

% PLOT DISTANCES
figure,
stem3(idx_2d_x, idx_2d_y, rays_distance)
xlabel('W'); ylabel('H');
axis equal

% Depth Image
M_dist = MAX_DEPTH*ones(array_size);
M_dist(idx_list) = rays_distance;
% Show image
im_sensor = M_dist; % rot90(rot90(M_dist));
% im_sensor = padarray(im_sensor,[1,1],'symmetric','both');   % post
clims = [0 max(rays_distance)];
% 2D image coloured %
figure,
imagesc(im_sensor)%, clims)
set(gca,'YDir','reverse');
axis([1 6 1 4]);
colorbar
xlabel('W'); ylabel('H');
axis equal
% 3D Surf Image %
[X,Y] = meshgrid(1:array_size(2), 1:array_size(1));
figure, 
surf(X,Y,im_sensor);
colorbar
xlabel('W'); ylabel('H');
axis equal

%%% PLANE %%% 
R_SAMPLES = 100;
range_arr = linspace(0, MAX_DEPTH, R_SAMPLES);

% 3D Image, where each layer is a plan
V = zeros(R_SAMPLES, array_size(2), array_size(1));
for c = 1:length(idx_list)
    
    nn = idx_2d_y(c);
    mm = idx_2d_x(c);
    
    [val, index] = min(abs(range_arr - rays_distance(c)));
    
    depth = zeros(1, R_SAMPLES); 
    depth(index) = 1;
    
    V(:, mm, nn) = depth;
end

% Slice for 3D vizualization 
[xx, yy, zz] = meshgrid(1:array_size(2), 1:R_SAMPLES, 1:array_size(1));
xslice = 1:array_size(2); 
yslice = 1:R_SAMPLES;
zslice = 1:array_size(1);
figure,
slice(xx,yy,zz, V,xslice,yslice, zslice)
colormap winter
axis equal%([1, array_dim(1), 0, MAX_DEPTH, 1 array_dim(2)])
xlabel('x'); ylabel('y'); zlabel('z');



