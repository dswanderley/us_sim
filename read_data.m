
close all; clear; clc;

fid = fopen('torus.obj');

v = [];

% parse .obj file 
while 1    
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end  % exit at end of file 
     ln = sscanf(tline,'%s',1); % line type 
     %disp(ln)
    switch ln
        case 'v'   % mesh vertexs
            v = [v; sscanf(tline(2:end),'%f')'];
            
    end
end

fclose(fid);

v(:,1) = v(:,1) - min(v(:,1));
v(:,2) = v(:,2) - min(v(:,2)) + 2;
v(:,3) = v(:,3) - min(v(:,3));

x_max = max(v(:,1));
x_min = min(v(:,1));

y_max = max(v(:,2));
y_min = min(v(:,2));

z_max = max(v(:,3));
z_min = min(v(:,3));