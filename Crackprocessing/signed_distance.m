function d = signed_distance(xCr,pt,truncation)

% Compute the signed distance from point x to crack xCr
% Inputs:
%     - xcR(2,2):    coordinates of points defining the crack
%     - pt(1,2):     coordinate of point 
%     - truncation:  1 to do the "truncation"   

%Author(s): Sundar, Stephane, Stefano, Phu, David
%Modified: Jan 2009
%--------------------------------------------------------------------


% select the right and left extremity of the crack to always have phi>0 on
% the same side of the crack
[val1,riga1] = min(xCr(:,1)); % min value on x
[val2,riga2] = max(xCr(:,1)); % max value on x

[val3,riga3] = min(xCr(:,2));   %min value on y
[val4,riga4] = max(xCr(:,2));   %max value on y

if(val3 == val4)
    x0  = val1; y0 = xCr(riga1,2);
    x1  = val2; y1 = xCr(riga2,2);
elseif(val1 == val2)
    x0 = xCr(riga3,1); y0 = val3;
    x1 = xCr(riga4,1); y1 = val4 ;
else
    x0  = val1; y0 = xCr(riga1,2);
    x1  = val2; y1 = xCr(riga2,2);
end

EPS = 1e-15;

x   = pt(1,1);
y   = pt(1,2);
%l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ;

if truncation == 1
    phi = Roundoffa(phi,6);
end

%d   = phi/l; 
d = phi;