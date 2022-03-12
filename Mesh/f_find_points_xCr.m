function [ inds ] = f_find_points_xCr( points, rift, radius, radius2)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Oct 26 16:46:29 NZDT 2021
% This function finds all the points that are within the radius to the rift (or point)
  if isstruct(rift)
    coords = rift.coor;
  end
  if nargin == 3
    p_buff = polybuffer(coords,'lines',radius)
  elseif nargin ==4
    p_buff1 = polybuffer(coords,'lines',radius)
    p_buff2 = polybuffer([coords(1,:);coords(end,:)],'points',radius2)
    p_buff = union(p_buff1,p_buff2)
  end
  inds = find( isinterior(p_buff, points(:,1), points(:,2)) ) 
end
