function [xp] = convertShelftoModel(x)
global Rtip QT xTip Tfact
% Model coordinates
x(:,1) = x(:,1)-Rtip(1);
x(:,2) = x(:,2)-Rtip(2);
xp = (QT'*x')./Tfact;      
xp(1,:) = xp(1,:) +xTip(1);
xp(2,:) = xp(2,:) +xTip(2);
xp = (xp)';
end
