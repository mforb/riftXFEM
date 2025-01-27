function [ phi ] = f_dista_point( pnt, iel, elem_crk  )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Jan  6 16:25:58 NZDT 2022
global epsilon
global node element

sctr = element(iel,:) ;

EPS = 1e-08 ;

x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;

x = pnt(1) ;
y = pnt(2) ;

phi = (y0-y1)*x + (x1-x0)*y + (x0*y1 - y0*x1) ;

if abs(phi) < epsilon 
  phi = 0 ;
end

end
