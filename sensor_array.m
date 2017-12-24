clear; close all; clc

ori = [1,0,0];
u = [0, 1, 0];
v = [0, 0, 1];

sensor_edge = 0.1;% in percent

array_H = 9;       % unit
array_W = 16;      % unit

array_W_max = 160;  % cm
array_H_max = 90;   % cm

sensor_W_max = array_W_max / array_W;
sensor_H_max = array_H_max / array_H;

% Define first center
first_center = [sensor_H_max/2, sensor_W_max/2];
% Calculate centers references in each direction
y_centers = linspace(first_center(1), array_H_max - first_center(1), array_H)';
x_centers = linspace(first_center(2), array_W_max - first_center(2), array_W);
% Repeat to each col/row
y_centers = repmat(y_centers, 1, array_W);
x_centers = repmat(x_centers, array_H, 1);
% Convert to world coords
X = ori(1) + (u(1) * x_centers) + (v(1) * y_centers);
Y = ori(2) + (u(2) * x_centers) + (v(2) * y_centers);
Z = ori(3) + (u(3) * x_centers) + (v(3) * y_centers);


sensor_H_max_util = (1 - sensor_edge)*sensor_H_max;
sensor_W_max_util = (1 - sensor_edge)*sensor_W_max;

sensor_shape = [9, 9];
center = [0, 0, 0];
max_dim = [sensor_H_max_util, sensor_W_max_util]/2;
sensor0 = rectsensor(sensor_shape, max_dim, center, u, v);

sensors = repmat(sensor0, 1, 1, numel(X));
x_z = reshape(X, 1, 1, numel(X));
y_z = reshape(Y, 1, 1, numel(Y));
z_z = reshape(Z, 1, 1, numel(Z));

sensors(:,1,:) = sensors(:,1,:) + x_z;
sensors(:,2,:) = sensors(:,2,:) + y_z;
sensors(:,3,:) = sensors(:,3,:) + z_z;

% Reshape to plot
xs = reshape(X,1,numel(X));
ys = reshape(Y,1,numel(Y));
zs = reshape(Z,1,numel(Z));
% Plot
figure, hold on
plot3(xs, ys, zs,'*')
plot3([1,1,1,1,1],[0,0, 160, 160, 0],[0,90,90,0, 0],'-g')
xlabel('x'); ylabel('y'); zlabel('z');
view(40,35); % view(-30,10);
axis equal;

x_sensor = reshape(sensors(:,1,:),1,numel(sensors(:,1,:)));
y_sensor = reshape(sensors(:,2,:),1,numel(sensors(:,2,:)));
z_sensor = reshape(sensors(:,3,:),1,numel(sensors(:,3,:)));

plot3(x_sensor,y_sensor,z_sensor,'.m')




