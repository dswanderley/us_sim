function [ rec ] = rectsensor(dim, max_dim,  center)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%close all; clear;clc;
c = center;%[0, 1, 1];
a = [0, 1, 0];
b = [0, 0, 1];

max_height = max_dim(1);% max_height = 0.7;
max_width = max_dim(2);% max_width = 0.3;
H = dim(1);% H = 10;
W = dim(2);% W = 7;

rows = linspace(-max_height, max_height, H);
cols = linspace(-max_width, max_width, W);

[X,Y] = meshgrid(cols, rows);

x = c(1) + (a(1) * X) + (b(1) * Y);
y = c(2) + (a(2) * X) + (b(2) * Y);
z = c(3) + (a(3) * X) + (b(3) * Y);

rec = [ reshape(x,numel(x),1), ...
        reshape(y,numel(y),1), ...
        reshape(z,numel(z),1)];

end
