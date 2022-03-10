function [ dist ] = f_dista2(iel,elem_crk,tip)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Nov 11 16:12:12 NZDT 2021

global node element epsilon

sctr = element(iel,:) ;

EPS = epsilon ;

if size(elem_crk,1)==1
  x0 = elem_crk(1) ; y0 = elem_crk(2) ;
  x1 = elem_crk(3) ; y1 = elem_crk(4) ;
else
  x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
  x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
end

if points_same_2d(tip,[x0,y0])
  x0 = x1;
  y0 = y1;
  x1 = tip(1);
  y1 = tip(2);
elseif ~points_same_2d(tip,[x1,y1])
  error('Tip and elem_crk dont match')
end

for i=1:size(sctr,2)
    x = node(sctr(i),1) ;
    y = node(sctr(i),2) ;
    psi = (x-x1)*(x1-x0)+(y-y1)*(y1-y0);
    if abs(psi) < EPS
     dist(i) = 0;
    else
     dist(i) = psi; 
    end
end

end
