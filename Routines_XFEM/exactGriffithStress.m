function [anaSigma] = exactGriffithStress(x,xTip,adv,sigmato,cracklength)

%adv - crack segment
%xTip - crack tip coordiante
%x - current point where the stress needs to be evaluated

% inclination of local coord
alfa = atan2(adv(2),adv(1));
% transpose of the rotation matrix
QT = [cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
% local coordinates
xp = QT*(x-xTip)';         
% local polar coordinates
[theta,r] = cart2pol(xp(1),xp(2));

sigma = sigmato;
KI    = sigmato*sqrt(pi*cracklength);   % exact KI

fac = KI/sqrt(r) ;

sigx = fac*cos(theta/2)*(1 - sin(theta/2)*sin(3*theta/2)) ;
sigy = fac*cos(theta/2)*(1 + sin(theta/2)*sin(3*theta/2)) ;
sigxy = fac*sin(theta/2)*cos(theta/2)*cos(3*theta/2) ;

anaSigma = [sigx; sigy; sigxy] ;