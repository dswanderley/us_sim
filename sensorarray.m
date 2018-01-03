function [sensors, centers, srect] = sensorarray(sensor_size, array_dim, array_size, u, v, origin)
%SENSORARRAY Summary of this function goes here
%   Detailed explanation goes here

sensor_edge = 0.1;          % in percent
sensor_shape = sensor_size; % Pixels in each dimension of each sensor
% Full array discret dimensions
array_H = array_size(1);	% unit
array_W = array_size(2);	% unit
% Full array dimensions
array_W_max = array_dim(2); % cm
array_H_max = array_dim(1); % cm
% Dimensions of each sensor
sensor_W_max = array_W_max / array_W;
sensor_H_max = array_H_max / array_H;
% Active dimensions of each sensor
sensor_H_max_util = (1 - sensor_edge)*sensor_H_max;
sensor_W_max_util = (1 - sensor_edge)*sensor_W_max;

% Define first center
first_center = [sensor_H_max/2, sensor_W_max/2];
% Calculate centers references in each direction
y_centers = linspace(first_center(1), array_H_max - first_center(1), array_H)';
x_centers = linspace(first_center(2), array_W_max - first_center(2), array_W);
% Repeat to each col/row
y_centers = repmat(y_centers, 1, array_W);
x_centers = repmat(x_centers, array_H, 1);
% Convert to world coords
X = origin(1) + (u(1) * x_centers) + (v(1) * y_centers);
Y = origin(2) + (u(2) * x_centers) + (v(2) * y_centers);
Z = origin(3) + (u(3) * x_centers) + (v(3) * y_centers);

% Set variables for one sensor on origin
center = [0, 0, 0];
max_dim = [sensor_H_max_util, sensor_W_max_util]/2;
% Calculate sensor on origin
sensor0 = rectsensor(sensor_shape, max_dim, center, u, v);

% Sensors for all centers
sensors = repmat(sensor0, 1, 1, numel(X));
x_z = reshape(X, 1, 1, numel(X));
y_z = reshape(Y, 1, 1, numel(Y));
z_z = reshape(Z, 1, 1, numel(Z));
% Final sensors (in third dimension)
sensors(:,1,:) = sensors(:,1,:) + repmat(x_z, size(sensors(:,1,:),1),1,1);
sensors(:,2,:) = sensors(:,2,:) + repmat(y_z, size(sensors(:,2,:),1),1,1);
sensors(:,3,:) = sensors(:,3,:) + repmat(z_z, size(sensors(:,3,:),1),1,1);
% Centers of the sensors
centers =  [reshape(X, 1, 1, numel(X)), ...
            reshape(Y, 1, 1, numel(Y)), ...
            reshape(Z, 1, 1, numel(Z))];

recW = [- array_W_max/2; array_W_max/2; array_W_max/2; - array_W_max/2];
recH = [- array_H_max/2; - array_H_max/2; array_H_max/2; array_H_max/2;];

medium = [mean(centers(:,1,:)), mean(centers(:,2,:)), mean(centers(:,3,:))];
Xrec = medium(1) + (u(1) * recW) + (v(1) * recH);
Yrec = medium(2) + (u(2) * recW) + (v(2) * recH);
Zrec = medium(3) + (u(3) * recW) + (v(3) * recH);
        
srect = [Xrec, Yrec, Zrec];

end

