function [ cutEdge, node ] = f_edgedetect( node, corner, phi, psi )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Nov 25 15:49:14 NZDT 2021
cutEdge = [ ] ;
%loop on element edges
if nargin == 3
for iedge = 1:size(node,1)
    n1 = corner(iedge) ;
    n2 = corner(iedge+1) ;
    if( phi(n1)*phi(n2) < 0 ) 
        r = phi(n1)/(phi(n1)-phi(n2)) ;
        pnt = (1-r)*node(n1,:)+r*node(n2,:) ;
        node = [node; pnt] ;
        %         disp(['Edge is cut by crack     ',num2str(iedg
        cutEdge = [cutEdge iedge] ;
    end
end %end iedge loop

elseif nargin == 4
    for i = 1:size(node,1)
        n1 = corner(i);
        n2 = corner(i+1);
        if phi(n1)*phi(n2) < 0  & ( psi(n1) < 0 | psi(n2) < 0) 
          l1 = psi(n1);
          l2 = psi(n2);
          r1 = phi(n1);
          r2 = phi(n2);
          t1 = 0;
          if (l1 < 0 ) & (l2 < 0)
            t1 = 1;
          elseif (l1 < 0 )
            if abs(l1/r1) > abs(l2/r2)
              t1 = 1;
            end
          else
            if abs(l2/r2) > abs(l1/r1)
              t1 = 1;
            end
          end
            
          if t1
            r    = r1/(r1-r2);
            pnt  = (1-r)*node(n1,:)+r*node(n2,:);
            node = [node;pnt];
          end
        end
        
    end

end

end
