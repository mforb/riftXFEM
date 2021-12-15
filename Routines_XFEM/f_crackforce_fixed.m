function [ Fc ] = f_crackforce_fixed( Fc,F_app, crack_lips, xCr, xCrl, xTip, pos,type_elem, enr_node, split_elem, vertex_elem, tip_elem )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon Nov 29 11:58:33 NZDT 2021
global node element elemType

elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);
% no point in including tip_elem. the forces are going to be the same and cancel out

for kk = 1:size(xCr,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii) ;
    sctr=element(iel,:);
    nn = length(sctr);
    p1 = crack_lips(ii,1:2,4,kk);
    p2 = crack_lips(ii,3:4,4,kk);
    p3 = crack_lips(ii,5:6,4,kk);
    n1 = crack_lips(ii,1:2,1,kk);
    n2 = crack_lips(ii,3:4,1,kk);
    n3 = crack_lips(ii,5:6,1,kk);
    dup1 = crack_lips(ii,1:2,2,kk);
    dup2 = crack_lips(ii,3:4,2,kk);
    dup3 = crack_lips(ii,5:6,2,kk);
    ddown1 = crack_lips(ii,1:2,3,kk);
    ddown2 = crack_lips(ii,3:4,3,kk);
    ddown3 = crack_lips(ii,5:6,3,kk);

    [fd_xy,fu_xy] = f_calc_crack_force([p1,p2],[ddown1,ddown2], [dup1,dup2],F_app);
    for p = 1:2
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
        else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
            [sctrn,xx] = find(element == sctr(in));
            [ele,xx] = find(type_elem(sctrn,:)==1);
            ref_elem = sctrn(ele);
        end
        % compute branch functions at Gauss point
        xCre  = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];
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

    if ismember(iel,vertex_elem)
      [fdx,fux] = f_calc_crack_force([p2,p3],[c_down2,cdown3], [c_up2, c_up3]);

      for p = 1:2
        if p ==1 
          [N,dNdxi] = lagrange_basis(elemType,n2);
          Fux = fu_xy(1:2);
          Fdx = fd_xy(1:2);
        else
          [N,dNdxi] = lagrange_basis(elemType,n3);
          Fux = fu_xy(3:4);
          Fdx = fd_xy(3:4);
        end
        for in = 1:nn
          nodeI = sctr(in) ;
          Ni = N(in);
          Fc(2*nodeI-1) = Fc(2*nodeI-1) + Ni*Fux(1) + Ni*Fdx(1)
          Fc(2*nodeI) = Fc(2*nodeI) + Ni*Fux(2) + Ni*Fdx(2)
          if (enr_node(nodeI) == 2) | (enr_node(nodeI) == 3)    % H(x) enriched node
            idx = pos(nodeI);
            Fc(2*idx-1) = Fc(2*idx-1) + Ni*Fux(1) 
            Fc(2*idx) = Fc(2*idx) + Ni*Fux(2) 
          elseif(enr_node(nodeI) == 1)
            idx = pos(nodeI);
            idxs = idx:1:idx+3;
            for id = 1:4
              idx = idxs(id);
              Fc(2*idx - 1) = Fc(2*idx-1) + Ni*Fux(1) + Ni*Fdx(1) 
              Fc(2*idx) = Fc(2*idx) + Ni*Fux(2) + Ni*Fdx(1)
            end
          end
        end
      end
    end
end

end
