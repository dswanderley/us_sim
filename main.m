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
orientation = [1,0,0];  % Array beams orientation
u = [0, 1, 0];          % Array plan vector 1
v = [0, 0, 1];          % Array plan vector 2
% Rectangular sensor
sensor_size = [5, 3];   % unit (number of beams)
sensor_beams = sensor_size(1) * sensor_size(2);
array_dim = [90, 160]/10;  % cm
array_size = [3, 5];    % unit (number of sensors in each direction)
num_sensors = array_size(1) * array_size(2);        % Number of sensors on array
% Get sensors
[sensors, centers] = sensorarray(sensor_size, array_dim, array_size, u, v, [0,0,0]);

% GET FOCAL POINTS
array_center = centers(:,:,1);
sensor_base = sensors(:,:,1);
% Angle - 0 for collinear
max_angle = 0;%30 * pi /180;
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

focus = repmat(focus_base, 1, 1, size(sensors,3));
focus(:,1,:) = focus(:,1,:) + repmat(centers(:,1,:), size(focus(:,1,:),1), 1, 1);
focus(:,2,:) = focus(:,2,:) + repmat(centers(:,2,:), size(focus(:,2,:),1), 1, 1);
focus(:,3,:) = focus(:,3,:) + repmat(centers(:,3,:), size(focus(:,3,:),1), 1, 1);

figure, % PLOT REAL WORLD
trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','red','FaceAlpha',0.1)
hold on
% PLOT SENSOR CENTERS
p_x = reshape(centers(:,1,:), 1, num_sensors); 
p_y = reshape(centers(:,2,:), 1, num_sensors);
p_z = reshape(centers(:,3,:), 1, num_sensors);
plot3(p_x, p_y, p_z, '*')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
% PLOT SENSOR BORDERS



%%% RAYS DISTANCE %%%
array_dists = zeros(1,size(sensors,3));
array_lines = zeros(2,3,size(sensors,3));
for k = 1:size(sensors,3)
    % Get all point in each sensor, their focus and the sensor center
    sensor = sensors(:,:,k);
    fc = focus(:,:,k);
    c = centers(:,:,k);
    % Calculate minimal distance from center
    [min_dist, p_sensor, p_voxel] = sensormindist(sensor, fc, c, obj);
	array_dists(:,k) = min_dist;
    line = [c; p_voxel]; 
    array_lines(:,:,k) = line;
    % PLOT MIN BEAM  
    plot3(line(:,1,:), line(:,2,:), line(:,3,:));
end



hold off

% PLOT DISTANCES
% % figure,
% % stem3(sensor(idx_list,2), sensor(idx_list,3), rays_distance)
% % axis equal
% % xlabel('x'); ylabel('y'); zlabel('z');
% 
% %[idx_l, idx_c] = ind2sub(sensor_shape, idx_list);
% im_sensor = zeros(sensor_shape);
% im_sensor(idx_list) = abs(rays_distance/max(rays_distance)-1)+0.5;
% im_sensor = rot90(rot90(im_sensor));
% figure, imshow(im_sensor)
% 
% %%% MIM DISTANCE %%%
% min_ray = min(rays_distance);
% min_idx = find(rays_distance == min_ray);
% min_idx = min_idx(1);
% 
% disp(['dist min: ', num2str(min_ray)])
% disp(['interception: ', num2str(rays_intercept(min_idx))])
% 
% 
% %%% PLANE %%% 
% MAX_RANGE = 3; % m
% R_SAMPLES = 100;
% range_arr = linspace(0, MAX_RANGE, R_SAMPLES);
% 
% n_sensors = sensor_shape(2);
% % 3D Image, where each layer is a plan
% V = zeros(R_SAMPLES, sensor_shape(2), sensor_shape(1));
% for c = 1:length(idx_list)
%     
%     idx = idx_list(c);
%     [nn, mm] = ind2sub(sensor_shape,idx);
%     [val, index] = min(abs(range_arr - rays_distance(c)));
%     
%     depth = zeros(1, R_SAMPLES); 
%     depth(index) = 1;
%     
%     V(:, mm, nn) = depth;
% end
% 
% % Slice for 3D vizualization 
% [xx, yy, zz] = meshgrid(1:sensor_shape(2) ,1:R_SAMPLES, 1:sensor_shape(1));
% xslice = 1:sensor_shape(2); 
% yslice = 1:R_SAMPLES;
% zslice = 1:sensor_shape(1);
% figure,
% slice(xx,yy,zz, V,xslice,yslice, zslice)
% colormap winter
% axis equal
% xlabel('x'); ylabel('y'); zlabel('z');
% 
% % 3D Volume vizualization
% xdata = [0,sensor_shape(2)]; 
% ydata = [0,R_SAMPLES];
% zdata = [0,sensor_shape(1)];
% 
% figure, 
% vol3d('cdata', V, 'xdata', xdata, 'ydata', ydata, 'zdata', zdata);
% colormap winter
% %axis equal off
% %set(gcf, 'color', 'w');
% %xlabel('x'); ylabel('y'); zlabel('z');




