clear; close all; clc;


o_shape = [8, 58; 0, 15; 0, 50];
max_dim = [diff(o_shape(1,:)); diff(o_shape(2,:)); diff(o_shape(3,:))];

array_dim = [120, 120]/10; 
s_lim = o_shape;
s_lim(:,1) = o_shape(:,1)-10;
s_lim(:,2) = o_shape(:,2)+10;


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

diag_factor = array_dim(1)*sin(pi/4)/2;
u1 = [0,1,0];   u2 = [-1,0,0];   u3 = [0,-1,0];   u4 = [1,0,0];
d1 = [sqrt(2)/2, sqrt(2)/2, 0]; d2 = [-sqrt(2)/2, sqrt(2)/2, 0];
d3 = [-sqrt(2)/2, -sqrt(2)/2, 0]; d4 = [sqrt(2)/2, -sqrt(2)/2, 0];
v = [0,0,1];
positions{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);
list_u{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);
list_v{numel(zz)} = zeros(numel(xx1)+numel(yy1)+numel(xx2)+numel(yy2)+4,3);

for k = 1:length(zz)
    z = zz(k);
      
        
    lateral_1 = zeros(numel(xx1), 3);
    lateral_1(:,1) = xx1';
    lateral_1(:,2) = s_lim(2, 1);
    lateral_1(:,3) = z;
    vec_u1 = repmat(u1, numel(xx1), 1);
    vec_v1 = repmat(v, numel(xx1), 1);
    
    lateral_2 = zeros(numel(yy1), 3);
    lateral_2(:,1) = s_lim(1, 2);
    lateral_2(:,2) = yy1';
    lateral_2(:,3) = z;
    vec_u2 = repmat(u2, numel(yy1), 1);
    vec_v2 = repmat(v, numel(yy1), 1);
    
    lateral_3 = zeros(numel(xx2), 3);
    lateral_3(:,1) = xx2';
    lateral_3(:,2) = s_lim(2, 2);
    lateral_3(:,3) = z;
    vec_u3 = repmat(u3, numel(xx2), 1);
    vec_v3 = repmat(v, numel(xx2), 1);
    
    lateral_4 = zeros(numel(yy2), 3);
    lateral_4(:,1) = s_lim(1, 1);
    lateral_4(:,2) = yy2';
    lateral_4(:,3) = z;
    vec_u4 = repmat(u4, numel(yy2), 1);
    vec_v4 = repmat(v, numel(yy2), 1);
    
    diagonal_1 = [s_lim(1,1), s_lim(2,1) + diag_factor, z];
    diagonal_2 = [s_lim(1,2) - diag_factor, s_lim(2,1), z];
    diagonal_3 = [s_lim(1,2), s_lim(2,2) - diag_factor, z];
    diagonal_4 = [s_lim(1,1) + diag_factor, s_lim(2,2), z];
    
    positions{k} =  [diagonal_1; lateral_1;...
                     diagonal_2; lateral_2;...
                     diagonal_3; lateral_3;...
                     diagonal_4; lateral_4];
    list_u{k} = [d1; vec_u1; d2; vec_u2; d3; vec_u3; d4; vec_u4];
    list_v{k} = [v;  vec_v1; v;  vec_v2; v;  vec_v3; v;  vec_v4];
end

figure,
for c = 1:16
    plot(c,c,'*');
    hold on
end

positions = cell2mat(positions(:));
list_u = cell2mat(list_u(:));
list_v = cell2mat(list_v(:));


vec1 = positions+list_u;
vec2 = positions+list_v;

figure,
xlabel('x');ylabel('y');zlabel('z');
hold on
    
for c = 1:size(positions,1)
    
    plot3(positions(c,1), positions(c,2), positions(c,3), '*')
    
    xxx = [positions(c,1), vec1(c,1)];
    yyy = [positions(c,2), vec1(c,2)];
    zzz = [positions(c,3), vec1(c,3)];
    plot3(xxx, yyy, zzz)
    
    
    xxx = [positions(c,1), vec2(c,1)];
    yyy = [positions(c,2), vec2(c,2)];
    zzz = [positions(c,3), vec2(c,3)];
    plot3(xxx, yyy, zzz)
    
end

plot3([8, 8, 58, 58, 8], [0, 15, 15, 0, 0], [0, 0,0,0,0])
plot3([s_lim(1,1),s_lim(1,1),s_lim(1,2),s_lim(1,2),s_lim(1,1)],...
      [s_lim(2,1), s_lim(2,2), s_lim(2,2),s_lim(2,1),s_lim(2,1)],...
      [0, 0,0,0,0], '--')
axis equal

plot3(diagonal_1(1),diagonal_1(2),diagonal_1(3),'bo')
plot3(diagonal_2(1),diagonal_2(2),diagonal_2(3),'ro')
plot3(diagonal_3(1),diagonal_3(2),diagonal_3(3),'go')
plot3(diagonal_4(1),diagonal_4(2),diagonal_4(3),'yo')

