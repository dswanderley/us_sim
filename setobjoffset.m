function [ M ] = setobjoffset( M, offsets, scale )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dim = length(offsets);

if size(M,2) == dim
    for k = 1:dim
        v_min = min(M(:,k));
        M(:,k) =  (M(:,k) - v_min) * scale + offsets(k);
    end 
end

end

