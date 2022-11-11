function [AuxStress,AuxGradDisp,AuxEps] = f_auxiliary(xp,r,theta,mu,kappa,mode)
K1 = 1.0 ;
K2 = K1  ;

SQR  = sqrt(r);
CT   = cos(theta);
ST   = sin(theta);
CT2  = cos(theta/2);
ST2  = sin(theta/2);
C3T2 = cos(3*theta/2);
S3T2 = sin(3*theta/2);

drdx = CT;
drdy = ST;
dtdx = -ST/r;
dtdy = CT/r;

FACStress1 = sqrt(1/(2*pi));
FACStress2 = FACStress1;

FACDisp1 = sqrt(1/(2*pi))/(2*mu);
FACDisp2 = FACDisp1;

AuxStress   = zeros(2,2);
AuxGradDisp = zeros(2,2);
AuxEps      = zeros(2,2);

if (mode == 1)
    
    AuxStress(1,1) = K1*FACStress1/SQR*CT2*(1-ST2*S3T2);
    AuxStress(2,2) = K1*FACStress1/SQR*CT2*(1+ST2*S3T2);
    AuxStress(1,2) = K1*FACStress1/SQR*ST2*CT2*C3T2;
    AuxStress(2,1) = AuxStress(1,2);
    
    u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
    du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
    du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);
    
    u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
    du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
    du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);
    
    AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
    AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
    AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
    AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
    
    AuxEps(1,1) = AuxGradDisp(1,1);
    AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
    AuxEps(1,2) = AuxEps(2,1);
    AuxEps(2,2) = AuxGradDisp(2,2);
    
elseif (mode == 2)
    AuxStress(1,1) = -K2*FACStress2/SQR*ST2*(2+CT2*C3T2);
    AuxStress(2,2) = K2*FACStress2/SQR*ST2*CT2*C3T2;
    AuxStress(1,2) = K2*FACStress2/SQR*CT2*(1-ST2*S3T2);
    AuxStress(2,1) = AuxStress(1,2);
    
    u1    = K2*FACDisp2*SQR*ST2*(kappa + 2 + CT);
    du1dr = K2*FACDisp2*0.5/SQR*ST2*(kappa + 2 + CT);
    du1dt = K2*FACDisp2*SQR*(0.5*CT2*(kappa + 2 + CT) - ST2*ST);
    
    u2    = -K2*FACDisp2*SQR*CT2*(kappa - 2 + CT);
    du2dr = -K2*FACDisp2*0.5*(1/SQR)*CT2*(kappa - 2 + CT);
    du2dt = -K2*FACDisp2*SQR*(-0.5*ST2*(kappa - 2 + CT) - CT2*ST);
    
    AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
    AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
    AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
    AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
    
    AuxEps(1,1) = AuxGradDisp(1,1);
    AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
    AuxEps(1,2) = AuxEps(2,1);
    AuxEps(2,2) = AuxGradDisp(2,2);
end
