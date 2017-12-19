%function [ tri ] = wavefront( points )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all; clear; clc;


[X,Y,Z] = meshgrid(sensor(:,1),sensor(:,2),sensor(:,3));
    
xs = [min(points(:,1)), max(points(:,1))];
ys = [min(points(:,2)), max(points(:,2))];
zs = [min(points(:,3)), max(points(:,3))];

[X,Y,Z] = meshgrid(xs,ys,zs);

x = reshape(X, numel(X),1);
y = reshape(Y, numel(Y),1);
z = reshape(Z, numel(Z),1);

p = [x,y,z];
unique(p,'rows')
 
tri = [];
if (xs(1) ~= xs(2) && ys(1) ~= ys(2))
    xy = unique([x,y],'rows');
    tri = [tri; delaunay(xy)];
end



if (xs(1) ~= xs(2) && zs(1) ~= zs(2))
    xz = unique([x,z],'rows');
    tri = [tri; delaunay(xz)];
end

if (zs(1) ~= zs(2) && ys(1) ~= ys(2))
	yz = unique([y,z],'rows');
    tri = [tri; delaunay(yz)];
end

find(yz==p(:,2:3))

%end

