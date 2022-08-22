function [ nu ] = f_norm( u, u_n, inds )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Dec  7 12:12:38 NZDT 2021
if nargin < 3
  inds = 1:length(u);
end
  
 nu = sum((u-u_n).^2)/sum((u.^2)+(u_n.^2));

end
