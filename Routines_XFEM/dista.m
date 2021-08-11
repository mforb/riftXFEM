function [dist] = dista(iel,elem_crk)

global node element

sctr = element(iel,:) ;

EPS = 1e-08 ;

x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;

for i=1:size(sctr,2)
    x = node(sctr(i),1) ;
    y = node(sctr(i),2) ;
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1 - y0*x1) ;
    if abs(phi) < EPS
        dist(i,1) = 0 ;
    else
        dist(i,1) = phi ;
    end
end