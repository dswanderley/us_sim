close all; clear; clc;
%% OBJ Settings %%
%obj_str = 'Estadio-do-dragao';
%obj_str = 'Igreja_Pampulha';
%obj_str = 'Apoteose';
%obj_str = 'torus-188';
obj_str = 'wheel_with_hat';

if  exist([obj_str ,'.mat'],'file')
    obj = load([obj_str, '.mat']);
else
    obj = readObj([obj_str, '.obj']);
    save([obj_str, '.mat'],'-struct', 'obj');   
end

%%
figure('units','normalized','outerposition',[0 0 1 1])
%trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','green','FaceAlpha',0.1)
trisurf(obj.f.v,obj.v(:,1),obj.v(:,2),obj.v(:,3),'Facecolor','blue','FaceAlpha',0.1)
axis equal
xlabel('x'), ylabel('y'), zlabel('z')
%view(3)
view([-37.5, -30])

% figure,
% plot3(obj.v(:,1),obj.v(:,2),obj.v(:,3),'*')
% axis equal
% view(2)
% xlabel('x'), ylabel('y'), zlabel('z')