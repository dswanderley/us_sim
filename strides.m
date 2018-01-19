function [ positions, list_v, list_u ] = strides(obj_v, array_dim, safe_guard)
%STRIDES Summary of this function goes here
%   Detailed explanation goes here

% Main shape to be scaned
o_shape = [ min(obj_v(:,1)), max(obj_v(:,1));...
            min(obj_v(:,2)), max(obj_v(:,2));...
            min(obj_v(:,3)), max(obj_v(:,3))];
% Max length in each dimension
max_dim = [diff(o_shape(1,:)); diff(o_shape(2,:)); diff(o_shape(3,:))];
% Define limites
s_lim = o_shape;
s_lim(:,1) = o_shape(:,1) - safe_guard;
s_lim(:,2) = o_shape(:,2) + safe_guard;

% X axis postions
residue_x_l = floor(mod(max_dim(1), array_dim(1))/2);
residue_x_r = mod(max_dim(1), array_dim(1)) - residue_x_l;
xx1 = o_shape(1,1)+residue_x_l : +array_dim(1) : o_shape(1,2)-residue_x_r;
xx2 = o_shape(1,2)-residue_x_r : -array_dim(1) : o_shape(1,1)+residue_x_l;
% Y axis postions
residue_y_l = floor(mod(max_dim(2), array_dim(1))/2);
residue_y_r = mod(max_dim(2), array_dim(1)) - residue_y_l;
yy1 = o_shape(2,1)+residue_y_l : +array_dim(1) : o_shape(2,2)-residue_y_r-1;
yy2 = o_shape(2,2)-residue_y_r-1 : -array_dim(1) : o_shape(2,1)+residue_y_l;
% Y axis positions
zz = o_shape(3,1) : array_dim(2) : o_shape(3,2)-array_dim(2);
% Diagonal positions
diag_factor = array_dim(1)*sin(pi/4)/2;
u1 = [0,1,0];   u2 = [-1,0,0];   u3 = [0,-1,0];   u4 = [1,0,0];
d1 = [sqrt(2)/2, sqrt(2)/2, 0]; d2 = [-sqrt(2)/2, sqrt(2)/2, 0];
d3 = [-sqrt(2)/2, -sqrt(2)/2, 0]; d4 = [sqrt(2)/2, -sqrt(2)/2, 0];
v = [0,0,1];
% Initialize postions
positions{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);
list_u{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);
list_v{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);
% Set postions acording Z plane
for k = 1:length(zz)
    z = zz(k);
    % First lateral
    lateral_1 = zeros(numel(xx1), 3);
    lateral_1(:,1) = xx1';
    lateral_1(:,2) = s_lim(2, 1);
    lateral_1(:,3) = z;
    vec_u1 = repmat(u1, numel(xx1), 1);
    vec_v1 = repmat(v, numel(xx1), 1);
    % Second Lateral
    lateral_2 = zeros(numel(yy1), 3);
    lateral_2(:,1) = s_lim(1, 2);
    lateral_2(:,2) = yy1';
    lateral_2(:,3) = z;
    vec_u2 = repmat(u2, numel(yy1), 1);
    vec_v2 = repmat(v, numel(yy1), 1);
    % Third Lateral
    lateral_3 = zeros(numel(xx2), 3);
    lateral_3(:,1) = xx2';
    lateral_3(:,2) = s_lim(2, 2);
    lateral_3(:,3) = z;
    vec_u3 = repmat(u3, numel(xx2), 1);
    vec_v3 = repmat(v, numel(xx2), 1);
    % Fourth Lateral
    lateral_4 = zeros(numel(yy2), 3);
    lateral_4(:,1) = s_lim(1, 1);
    lateral_4(:,2) = yy2';
    lateral_4(:,3) = z;
    vec_u4 = repmat(u4, numel(yy2), 1);
    vec_v4 = repmat(v, numel(yy2), 1);
    % Diagonals
    diagonal_1 = [s_lim(1,1), s_lim(2,1) + diag_factor, z];
    diagonal_2 = [s_lim(1,2) - diag_factor, s_lim(2,1), z];
    diagonal_3 = [s_lim(1,2), s_lim(2,2) - diag_factor, z];
    diagonal_4 = [s_lim(1,1) + diag_factor, s_lim(2,2), z];
    % Set Postions
    positions{k} =  [diagonal_1; lateral_1;...
                     diagonal_2; lateral_2;...
                     diagonal_3; lateral_3;...
                     diagonal_4; lateral_4];
    list_u{k} = [d1; vec_u1; d2; vec_u2; d3; vec_u3; d4; vec_u4];
    list_v{k} = [v;  vec_v1; v;  vec_v2; v;  vec_v3; v;  vec_v4];
end
% Convert to matix
positions = cell2mat(positions(:));
list_u = cell2mat(list_u(:));
list_v = cell2mat(list_v(:));

end

