function [ang] = f_MTS_adj( KI,KII,T,Y,d)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022
ad =[];
adt =[];
adtl =[];
adt2 = [];
adtl2 = [];
%T = T*sqrt(2*pi*d);
%Y = Y*sqrt(2*pi*d);
KI = KI/(sqrt(2*pi*d));
KII = KII/(sqrt(2*pi*d));
for t = -pi/2:0.005:pi/2 
  d = (1/2)*cos(t/2)*(KI*cos(t/2)^2-(3/2)*KII*sin(t))+T*sin(t)^2+Y*cos(t)^2; 
  %dt = 1/4* cos(t/2)* (-3* KI* sin(t) - 9* KII* cos(t) + 3* KII - 8* T *sin(t/2) + 8* T* sin((3 *t)/2));
  dt = (-3/2)*KII *cos((1/2)* t)* cos(t) - 3/2* KI* cos((1/2)* t)^2 * sin((1/2)* t) + 2 *T* cos(t) *sin(t) - 2* Y *cos(t) *sin(t) + 3/4* KII* sin((1/2 )*t) *sin(t);

  %dt = dt/sqrt(KI*KI+KII*KII);
  %dt = 1/2* cos(t/2) *(-3* KI* sin(t) - 9* KII* cos(t) + 3*KII - 4* T* sin(t/2) + 4* T* sin((3* t)/2));
  dtl = KI*sin(t)+KII*(3*cos(t)-1) - (16/3)*T*cos(t)*sin(t/2)+ (16/3)*Y*cos(t)*sin(t/2);


  dt2 = (3/8)* cos(t/2) *(4* KI* sin(t/2)^2 + 5* KII* sin(t)) - 3/4* KI *cos(t/2)^3 + cos(t)* ((3/2)* KII* sin(t/2) + 4* T* sin(t)) - 2* t* T* sin(t)^2 + 2* t* T* cos(t)^2;
  dtl2 = KI* cos(t/2)*(3*cos(t)-1) - KII* sin(t/2)*(9*cos(t)+5) - (16/3)*T*cos(2*t) + (16/3)*Y*cos(2*t);
  ad = [ad,d];
  adt = [adt,dt];
  adtl = [adtl,dtl];
  adt2 = [adt2,dt2];
  adtl2 = [adtl2,dtl2];
end
%df = [adt(2:end),adt(end)]-adt
th = -pi/2:0.005:pi/2;
inds = find(adtl2>0);
th = th(inds);
adtl = adtl(inds);

[v,in] = min(abs(adtl));
if abs(v)/sqrt(KI*KI+KII*KII)> 0.1
  ang = NaN;
else
  ang = th(in); 
  if abs(ang) < 0.001
    ang = 0
  end
end
%figure(1)
%clf
%hold on
%plot(th,ad,'g-');
%plot(th,adt,'r-');
%plot(th,adtl,'b-');
%plot(th,adt2,'k--');
%plot(th,adtl2,'c--');
%keyboard
