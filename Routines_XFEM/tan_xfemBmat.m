function [Bfem,Bxfem] = tan_xfemBmat(pt,e,tan_element,xCrl,crack_node)

%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon

sctr = tan_element(e,:);
nn   = length(sctr);
[N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
Gpt = N' * node(sctr,:);                  % GP in global coord, used

% standard B
Bfem = zeros(3,2*nn) ;
Bfem(1,1:2:2*nn) = dNdx(:,1)' ;
Bfem(2,2:2:2*nn) = dNdx(:,2)' ;
Bfem(3,1:2:2*nn) = dNdx(:,2)' ;
Bfem(3,2:2:2*nn) = dNdx(:,1)' ;

Bxfem = [];

for in = 1:nn
  ref_elem = e;
  xCre = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];                 %each element has its crack!
  % Enrichment function, H(x) at global Gauss point
  dist = signed_distance(xCre,Gpt,16);
  Hgp  = sign(dist);
  % Enrichment function, H(x) at node "in"
  dist = signed_distance(xCre,node(sctr(in),:),0);
  Hi  = sign(dist);
  if ismember(sctr(in),crack_node) 
    %Hi = sign(-1);
    %Hi = -1*Hgp;
    Hi = 0;
  end
  % Bxfem at node "in"
  BI_enr = [dNdx(in,1)*(Hgp - Hi) 0 ; 0 dNdx(in,2)*(Hgp - Hi) ;
      dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
  %else    %if the element is not cut by a crack, the enrichment is always 0 (NO LONGER TRUE)
      %BI_enr = [0 0 ; 0 0 ; 0 0];
  %end
  % Add to the total Bxfem
  Bxfem = [Bxfem BI_enr];
  clear BI_enr ;
end

%if e == 120
  %keyboard
%end
