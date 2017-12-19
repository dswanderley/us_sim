close all; clear; clc;

obj_str = 'torus-188';

if  exist([obj_str ,'.mat'],'file')
    x_min = min(obj.v(:,1));
    x_offset = 5;
    obj.v(:,1) = obj.v(:,1) - x_min + x_offset;
    obj = load('torus.mat');
else
    obj = readObj([obj_str, '.obj']);
    save('obj.mat','-struct', 'obj');
end

scale = 0.01;
offsets = [1, 0, 0];
obj.v = setobjoffset(obj.v, offsets, scale);

trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','red','FaceAlpha',0.1)
xlabel('x'); ylabel('y'); zlabel('z');
% Feixe paralelo ao eixo x
%     z
%      |__ y
%   x /
p_ref = [ 0, 1, 1];          % Transducer center
p_foc = [offsets(1), 1, 1];   % Focal point
line = [p_ref; p_foc];

[ d_min, p_int,  face_int_idx ] = minvoxeldist( obj, line );

disp(['dist min: ', num2str(d_min)])
disp(['interception: ', num2str(p_int)])
