function [v] = getobj(str)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(str);

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
        case 'vt'  % texture coordinate
            vt = [vt; sscanf(tline(3:end),'%f')'];
        case 'vn'  % normal coordinate
            vn = [vn; sscanf(tline(3:end),'%f')'];
        case 'f'   % face definition
            fv = []; fvt = []; fvn = [];
            
            str = textscan(tline(2:end),'%s'); str = str{1};
            disp(str);
    end
end

fclose(fid);

v(:,1) = v(:,1) - min(v(:,1));
v(:,2) = v(:,2) - min(v(:,2)) + 2;
v(:,3) = v(:,3) - min(v(:,3));
end

