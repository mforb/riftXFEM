function [H] = f_extractHeight(x)
global FintH same_coords
% The applied stress field needs to be previously defined as three matlab interpolants
% QT is the rotation matrix between the field and model

if exist('same_coords') & same_coords
  xint = x
else
  xint = convertModeltoShelf(x);
end
H = FintH(xint); 
end
