%function [ tri ] = wavefront( points )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all; clear; clc;

center = [0, 0.5, 0.75];
sensor_shape = [10, 4];
points = rectsensor(sensor_shape, sensor_shape/10,  center);
    
xs = [min(points(:,1)), max(points(:,1))];
ys = [min(points(:,2)), max(points(:,2))];
zs = [min(points(:,3)), max(points(:,3))];

[X,Y,Z] = meshgrid(xs,ys,zs);

x = reshape(X, numel(X),1);
y = reshape(Y, numel(Y),1);
z = reshape(Z, numel(Z),1);


if (xs(1) == xs(2))
    p = [y,z];
    p = unique(p,'rows');
    tri = delaunay(p);
    p = [xs(1)*ones(size(p,1),1), p];
elseif (ys(1) == ys(2))
    p = [x,z];
    p = unique(p,'rows');
    tri = delaunay(p);
    p = [xs(1)*ones(size(p,1),1), p];
elseif (zs(1) == zs(2))
    p = [x,y];
    p = unique(p,'rows');
    tri = delaunay(p);
    p = [xs(1)*ones(size(p,1),1), p];
else
    p = [x,y,z]; 
    tri = [];
end

trisurf(tri,p(:,1), p(:,2), p(:,3),'Facecolor','red','FaceAlpha',0.1)
axis equal


%end

