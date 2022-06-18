function [ flag1, flag2, flag3, crack_node, d_r ] = crack_interact_element( q, e, crack_node  )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue March 15 NZDT 2022

global node element  epsilon

sctr = element(e,:) ;
q1 = q(1:2);
q2 = q(3:4);
[ psi1, mps1 ] = f_dista2(e,[q1,q2],q1);
[ psi2, mps2 ] = f_dista2(e,[q2,q1],q2);
flag1 = 0;
flag2 = 0;
flag3 = 0;
ni = [];

if any(psi1 < 0 ) & any(psi2 < 0 ) 
  [ phi  ] = dista(e,[q1,q2]);
  for i = 1:length(sctr)
    if abs(phi(i)) < epsilon
      if psi1(i) < 0 & psi2(i) < 0
        ni = [ni,i];
        flag2 = 1;
        crack_node = [ crack_node, sctr(i) ]; % the element has a cracknode within the segment
      end
    end
  end
  phi(ni) = []; % this prevents crack_nodes from being used in determining if an element spans both sides of the crack
  if any(diff(sign(phi)))
    flag1 = 1; %this element does interat with this segment
    if any(psi1 > 0)
      vv = node(sctr,:);
      flag3 = inhull(q1,vv,[]) ; % one of the crack segment ends is within the element 
    elseif any(psi2>0) 
      vv = node(sctr,:);
      flag3 = inhull(q2,vv,[]) ; % one of the crack segments ends is within the element 
    end
    if sum(psi1 > 0) > 1 | sum(psi2 > 0) > 1
      if ~flag3 
        flag1 = 0;
      end
    end
  end
end

if nargout > 4
  d_r = [ mps1/(mps1+mps2), mps2/(mps1 + mps2)];
end
