function [ang] = f_SED_adj( KI,KII,T,k)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022
adt =[];
adtl =[];
adt2 = [];
adtl2 = [];
Keff = sqrt(KI*KI+KII*KII)
%T = T/sqrt(2*pi*d);
B = T/Keff;


for t = -pi/2:0.005:pi/2 
  ct = cos(t);
  st = sin(t);
  c2t = cos(2*t);
  ct2 = cos(t/2);
  st2 = sin(t/2);
  C1 = (1/16) * st*(2*ct-k+1);
  C2 = -1*(1/16) * st*(6*ct-k+1);
  C3 = (1/8) * (2*c2t-(k-1)*ct);
  C4 = -1*(1/16) *st2* (5*(c2t+ct)+(k+1));
  C5 = -1*(1/16) *ct2* (5*(c2t-ct)+(k+3));

  dt = C1*KI*KI+C2*KII+KII+C3*KI*KII + C4*T*KI +C5*T*KII;
  %dt = dt/sqrt(KI*KI+KII*KII);
  %dt = 1/2* cos(t/2) *(-3* KI* sin(t) - 9* KII* cos(t) + 3*KII - 4* T* sin(t/2) + 4* T* sin((3* t)/2));


  dtl2 = KI* cos(t/2)*(3*cos(t)-1) - KII* sin(t/2)*(9*cos(t)+5) - (16/3)*T*cos(2*t);
  adt = [adt,dt];
  adtl2 = [adtl2,dtl2];
end
th = -pi/2:0.005:pi/2;
inds = find(adtl2<0);
keyboard
th = th(inds);
adt = adt(inds);

[v,in] = min(abs(adt));
if abs(v)/sqrt(KI*KI+KII*KII)> 0.1
  ang = NaN;
else
  ang = th(in); 
end
