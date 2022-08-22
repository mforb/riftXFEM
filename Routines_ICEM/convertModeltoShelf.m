function [xp] = convertModeltoShelf(x)
global Rtip QT xTip Tfact
% Model coordinates
x(:,1) = x(:,1)-xTip(1);
x(:,2) = x(:,2)-xTip(2);
xp = Rtip' + (QT*(x)').*Tfact;       
xp = (xp)';
end
