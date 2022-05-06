function  [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enrich_node)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: April 3 2022

global node element skip_branch

sctr = element(iel,:) ;
nn = length(sctr) ;
n1 = zeros(1,nn);
A = [];
BrI = [];
Tip = [];
alpha = 0;
QT = eye(2);
for nI = 1:nn
  nodeI = sctr(nI);
  if (enrich_node(nodeI) == 2)     % H(x) enriched node
      AA = [2*pos(nodeI)-1;2*pos(nodeI)];
      A  = [A;AA];
  elseif(enrich_node(nodeI) == 3) % H(x) due to material
      AA = [2*pos(nodeI)-1; 2*pos(nodeI)];
      A = [A;AA];
  elseif (enrich_node(nodeI) == 1) % B(x) enriched node only the first branch function is relevant
     if ~skip_branch
        AA = [2*pos(nodeI)-1;
            2*pos(nodeI);
            %2*(pos(nodeI)+1)-1;
            %2*(pos(nodeI)+1);
            %2*(pos(nodeI)+2)-1;
            %2*(pos(nodeI)+2);
            %2*(pos(nodeI)+3)-1;
            %2*(pos(nodeI)+3);
            ];
        A  = [A;AA];
        if type_elem(iel,1) == 1   %looking for the "tip" element
            ref_elem = iel;
        else  
            [sctrn,xx] = find(element == nodeI);
            [ele,xx] = find(type_elem(sctrn,:)==1);
            ref_elem = sctrn(ele);
        end
        % compute branch functions
        if ~any(n1) % do this only once
          xCre  = [elem_crk(ref_elem,1:2); elem_crk(ref_elem,3:4)];
          seg   = xCre(2,:) - xCre(1,:);
          alpha = atan2(seg(2),seg(1));
          Tip  = [xCre(2,1) xCre(2,2)];
          QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
        end
        n1(nI) = 1;
        xp    = QT*(node(nodeI,:)-Tip)';
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1)); % presumably none of the nodes are on the crack near the tip
        BrI = [BrI; branch_node(r,theta)];
    end
  end
end
