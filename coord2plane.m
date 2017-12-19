function [const, bias] = coord2plane(matrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if size(matrix ~= [3,3])
    error('Error. \nmatrix must be a [3,3] element')
end

P = matrix(1,:);
Q = matrix(2,:);
R = matrix(3,:);

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

const = [a, b, c];
bias = [P(1), P(2), P(3)]';

end

