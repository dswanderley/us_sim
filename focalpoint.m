function [ pf ] = focalpoint( pc, pm, ang, u, v )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% % Central Point
% pc = [0,0,0];
% % Distant point
% pm = [0,1,0];
% % Angle
% ang = pi*30/180;
% % Plan vectors
% u = [0, 1, 0];
% v = [0, 0, 1];

% CosTheta = dot(u,v)/(norm(u)*norm(v));
% ThetaInDegrees = acosd(CosTheta);
d = norm(pm - pc);      % Distance
cat = d / tan(ang);     % Magnitude of the adjacente cathetus

% Det
x = u(2) * v(3) - v(2) * u(3); % i
y = u(3) * v(1) - v(3) * u(1); % j
z = u(1) * v(2) - v(1) * u(2); % k

% Focal point
pf = pc - [x,y,z]*cat;



end

