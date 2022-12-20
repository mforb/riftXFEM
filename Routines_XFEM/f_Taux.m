function [AuxStress,AuxGradDisp,AuxEps] = f_Taux(xp,r,d,theta,mu,kappa,mode)
F = 1.0 ;

CT   = cos(theta);
C2T   = cos(theta)^2;
C3T   = cos(theta)^3;
C3T   = cos(theta);
ST   = sin(theta);
S2T   = sin(theta)^2;
CT2  = cos(2*theta);
ST2  = sin(2*theta);

drdx = CT;
drdy = ST;
dtdx = -ST/r;
dtdy = CT/r;

FACStress1 = 1/(pi*r);

FACDisp1 = 1/(8*pi*mu));

AuxStress   = zeros(2,2);
AuxGradDisp = zeros(2,2);
AuxEps      = zeros(2,2);
if mode == 2
  return
end

    
AuxStress(1,1) = -1*FACStress1*CT3;
AuxStress(2,2) = -1*FACStress1*CT*S2T;
AuxStress(1,2) = -2*FACStress1*C2T*ST;
AuxStress(2,1) = AuxStress(1,2);

u1    = -F*FACDisp1*(kappa + 1)*ln(r)-2*FACDisp1*S2T;
u2    = -F*FACDisp1*(kappa - 1)*theta-2*FACDisp1*ST*CT;

AuxGradDisp(1,1) = (F/r)*FACDisp1*CT*(-1*kappa+1-2*CT2);
AuxGradDisp(1,2) = (F/r)*FACDisp1*ST*(-1*kappa-3-2*CT2);
AuxGradDisp(2,1) = (F/r)*FACDisp1*ST*(kappa-1-2*CT2);
AuxGradDisp(2,2) = (F/r)*FACDisp1*CT*(kappa-1-2*CT2);

AuxEps(1,1) = AuxGradDisp(1,1);
AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
AuxEps(1,2) = AuxEps(2,1);
AuxEps(2,2) = AuxGradDisp(2,2);
    
