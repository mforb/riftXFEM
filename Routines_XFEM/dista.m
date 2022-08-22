function [dist] = dista(iel,elem_crk)

global node element  

sctr = element(iel,:) ;

EPS = 1e-8 ;
if size(elem_crk,1)==1
  x0 = elem_crk(1) ; y0 = elem_crk(2) ;
  x1 = elem_crk(3) ; y1 = elem_crk(4) ;
else
  x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
  x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
end

l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;

for i=1:size(sctr,2)
    x = node(sctr(i),1) ;
    y = node(sctr(i),2) ;
    %phi = (y0-y1)*x + (x1-x0)*y + (x0*y1 - y0*x1) ;
    phi = (y0-y1)*(x-x0) + (x1-x0)*(y-y0);
    phi = phi/l;
    if abs(phi) < EPS
     dist(i) = 0 ;
    else
     dist(i) = phi ; 
    end
end

