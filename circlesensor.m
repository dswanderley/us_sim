function p_vec = circlesensor(r_max, N_r, N_ang, center, u, v)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% close all; clear; clc;
% r_max=30; N_r=36; N_ang=360;

c = center;%[0, 0, 0];
% u = [0, 1, 0];
% v = [0, 0, 1];

r_stp = r_max / N_r;
rr = flip(0 : r_stp : r_max);
aa = flip(0:N_r);

stp = 2 * pi / (N_ang);
th_in = 0;
th_out = 2*pi - stp;

p_vec = center;
M{N_r} = 0;
for k = 1:N_r
    
    th = th_in : stp : th_out;
    
    r = rr(k);
    r_c = r * cos(th);
    r_s = r * sin(th);
    
    x = c(1) + (u(1) * r_c) + (v(1) * r_s);
    y = c(2) + (u(2) * r_c) + (v(2) * r_s);
    z = c(3) + (u(3) * r_c) + (v(3) * r_s);
    
    p = [x', y', z'];
    
    M{k} = p;
     
    stp = 2 * pi / (N_ang * aa(k+1)/N_r);
    th_in = th_in + stp;
    th_out = th_in + 2*pi - stp;

end

M = cell2mat(M(:));
p_vec = [M; p_vec];

%plot(p_vec(:,2),p_vec(:,3),'.')

end

