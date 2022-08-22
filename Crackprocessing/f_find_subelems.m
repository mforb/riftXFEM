function [subelems, tangent_elem, corner_elem, type_elem, elem_crk, crack_node, kj_track ] = find_sub_elems(e,xCr,kk, subelems, tangent_elem, corner_elem, type_elem, elem_crk, crack_node ) ;
kj_track = [];
t = size(xCr.coor,1);
for kj = 1:size(xCr.coor,1)-1       %loop over the elements of the fracture
  t = t - 1;
  if t == 0 
    break;
  end
  q1 = xCr.coor(kj,:); 
  q2 = xCr.coor(kj+1,:);
  [flag1,flag2,flag3,crack_node] = crack_interact_element([q1,q2],e,crack_node);
  if flag1
    subelems = [subelems,e];
    if ~flag3
      kj_track = [kj_track, kj];
    end
    t = 2;
    if flag2
      corner_elem = [corner_elem,e];
    end
  elseif flag2 % but not flag 1 !
    kj_track = [kj_track, kj];
    t = 2;
    tangent_elem = [tangent_elem,e]; % this already takes care of all the tangent elements and (most likely?) all the crack nodes
    type_elem(e,kk) = 4; 
    elem_crk(e,:) = [q1,q2];
  end
end
%if ( e == 1715 | e == 2183 )  
  %keyboard
%end

