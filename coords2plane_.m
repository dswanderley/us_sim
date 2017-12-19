close all; clear; clc;


% P = [ 1, 0, 2];
% Q = [-1, 1, 2];
% R = [ 5, 0, 3];
P = [ 0, 0, 1];
Q = [ 0, 0, 0];
R = [ 0, 1, 0];


% a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
% (x0, y0, z0) is a point on plane
% <a, b, c>: perpendicular to plane

PQ = Q-P;
PR = R-P;
M = [PQ; PR];


M_i = M(:,2:end);
M_j = [M(:,1), M(:,3)];
M_k = M(:,1:2);

a = det(M_i);
b = det(M_j);
c = det(M_k);


x = 0; y = 0.5; z = 0.5;

plan = a*(x - P(1)) + b*(y - P(2)) + c*(z - P(3));
disp(plan)
