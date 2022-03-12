function [sigma] = f_extractStress(x)
global FintXY FintX FintY QT
% The applied stress field needs to be previously defined as three matlab interpolants
% QT is the rotation matrix between the field and model


xint = convertModeltoShelf(x);
xy = FintXY(xint); 
St = [ FintX(xint) xy; xy FintY(xint) ];

% Model coordinates
Stp = QT'*St*QT;
sigma = [Stp(1,1), Stp(2,2), Stp(1,2)];
end
