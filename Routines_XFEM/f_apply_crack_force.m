function [ Fc ] = f_apply_crack_force( Fc, fd_xy, fu_xy, iel, xCr, xCrl, xTip, pos, type_elem, enr_node, n1, n2)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Jan  6 15:09:05 NZDT 2022
global elemType
global element node

sctr=element(iel,:);
nn = length(sctr);
if nargin == 11
  pn = 1;
elseif nargin==12
  pn = 1;
end

for p = 1:pn
  if p ==1 
    gpt = n1;
    Fux = fu_xy(1:2);
    Fdx = fd_xy(1:2);
  else
    gpt = n2;
    Fux = fu_xy(3:4);
    Fdx = fd_xy(3:4);
  end
  [N,dNdxi] = lagrange_basis(elemType,gpt);
  % if there is  a tip node then we have to calculated Branch function for the point
  if any(enr_node(sctr)==1)
    in = find(enr_node(sctr)==1,1);
    if type_elem(iel,1) == 1   %looking for the "tip" element
        ref_elem = iel;
    else  
        [sctrn,xx] = find(element == sctr(in));
        [ele,xx] = find(type_elem(sctrn,:)==1);
        ref_elem = sctrn(ele);
    end
    % compute branch functions at Gauss point
    if points_same_2d(xCrl(ref_elem,3:4),xTip(ref_elem,:),1e-6)   
      xCre  = [ xCrl(ref_elem,1:2); xCrl(ref_elem,3:4) ]; 
    else
      xCre  = [ xCrl(ref_elem,3:4); xCrl(ref_elem,1:2) ]; 
    end
    seg   = xCre(2,:) - xCre(1,:);
    alpha = atan2(seg(2),seg(1));
    Tip  = [xCre(2,1) xCre(2,2)];
    QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
    xp    = QT*(gpt-Tip)';           % local coordinates
    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
    theta = atan2(xp(2),xp(1));
    
    if ( theta > pi | theta < -pi)
        disp (['something wrong with angle ',num2str(thet)]);
    end
    if abs(abs(theta) - pi) < 0.001
     [Br_u,dBdx,dbdy] = branch_gp(r,pi,alpha);
     [Br_d,dBdx,dbdy] = branch_gp(r,-1*pi,alpha);
    else 
     [Br_u,dbdx,dbdy] = branch_gp(r,theta,alpha);
     Br_d = Br_u;
    end
  end

  for in = 1:nn
    nodeI = sctr(in) ;
    Ni = N(in);
    Fc(2*nodeI-1) = Fc(2*nodeI-1) + Ni*Fux(1) + Ni*Fdx(1);
    Fc(2*nodeI) = Fc(2*nodeI) + Ni*Fux(2) + Ni*Fdx(2);
    if (enr_node(nodeI) == 2) | (enr_node(nodeI) == 3)    % H(x) enriched node
      idx = pos(nodeI);
      Fc(2*idx-1) = Fc(2*idx-1) + Ni*Fux(1); 
      Fc(2*idx) = Fc(2*idx) + Ni*Fux(2);
    elseif(enr_node(nodeI) == 1)
      xp    = QT*(node(sctr(in),:)-Tip)';
      r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
      theta = atan2(xp(2),xp(1)); % presumably none of the nodes are on the crack near the tip
      % we should make this more general in the future
      
      
      if ( theta > pi | theta < -pi)
          disp (['something wrong with angle ',num2str(thet)]);
      end
      [BrI] = branch_node(r,theta);
      idx = pos(nodeI);
      for i = 1:4
        Fc(2*idx - 1) = Fc(2*idx-1) + Ni*(Br_u(i)-BrI(i))*Fux(1) + Ni*(Br_d(i)-BrI(i))*Fdx(1); 
        Fc(2*idx) = Fc(2*idx) + Ni*(Br_u(i)-BrI(i))*Fux(2) + Ni*(Br_d(i)-BrI(i))*Fdx(2); 
        idx = idx + 1;
      end
    end
  end
end

end
