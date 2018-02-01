close all; clear; clc;
%% OBJ Settings %%
obj_str = 'torus-188';

if  exist([obj_str ,'.mat'],'file')
    obj = load([obj_str, '.mat']);
else
    obj = readObj([obj_str, '.obj']);
    save([obj_str, '.mat'],'-struct', 'obj');
end

scale = 0.5;
offsets = [8, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);
%%
%     z
%      |__ y
%   x /
%
%%
load('sim_all.mat');
total_points = all_sim;
%%
c_count = zeros(numel(total_points),1);
d_count = zeros(numel(total_points),1);
distances = zeros(numel(total_points),1);
%%
for c = 1:numel(total_points)
    
    sensor1 =total_points{1};
    sensor2 = total_points{c};
    
    colli=0;
    dvrgt=0;
    delta = 0;
    for s = 1:numel(sensor1)
        
        pts1 = sensor1{s};
        pts2 = sensor2{s};
        
        for l = 1:numel(pts1)
            p1 = pts1{l};
            p2 = pts2{l};
            if ~isempty(p1) && ~isempty(p2)
                colli = colli +1;
                dvrgt = dvrgt +1;       
                d = sqrt(sum((p2-p1).^2));
                delta = delta + exp(-d);
            elseif ~isempty(p1) && isempty(p2)
                colli = colli +1;
                delta = delta + 0;
            elseif isempty(p1) && ~isempty(p2)
                dvrgt = dvrgt +1;
                delta = delta + 0;
            end
        end
    end
    c_count(c) = colli;
    d_count(c) = dvrgt;
    distances(c) = delta;
end

metric = 2*distances./(c_count + d_count)