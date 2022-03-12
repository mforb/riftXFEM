function [H] = f_extractHeight(x)
global FintH
% The applied stress field needs to be previously defined as three matlab interpolants
% QT is the rotation matrix between the field and model


xint = convertModeltoShelf(x);
H = FintH(xint); 
end
