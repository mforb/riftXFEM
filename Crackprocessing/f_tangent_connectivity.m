function [tan_elem, tan_elem_crk] = f_tangent_connectivity( tangent_elem,crack_node,enrich_node,elem_crk)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Wed Nov 10 11:43:15 NZDT 2021

global node element numnode numelem elemType
global orig_nn
% crack_nodes has columns: detected node, elem_crack coord 1 x, elem_crack coord 1 y, elem_crack coord 2 x, elem_crack coord 2 y,
orig_nn = numnode;
n_red = 0;
already_done = [];
tan_elem = [];
tan_elem_crk = [];
for i = 1:length(tangent_elem)
  % test if adjacent to other tangent_elem
  if ~ismember(i,already_done)
    sctr = element(tangent_elem(i),:);
    nsctr = setdiff(sctr,crack_node);
    ri = [];
    for n = 1:length(nsctr)
      [r,c] = find(element(tangent_elem,:)==nsctr(n));
      ri = [ri;r];
    end
    ri = unique(ri);
    already_done = [already_done, ri'];
    % if so then merge
    if length(ri)>1
      all_nodes = unique(element(tangent_elem(ri),:));
      osctr = [];
      for ns = 1:length(all_nodes)
        if length(find(element(tangent_elem(ri),:) == all_nodes(ns))) == 1
          osctr =  [ osctr, all_nodes(ns) ]; 
        elseif ismember(all_nodes(ns),crack_node)
          osctr =  [ osctr, all_nodes(ns) ]; 
        end
      end
      if length(osctr) ~= 3
        keyboard
        error('tangent element merge resulted in non-triangular element')
      end
      tri = [ 1 2 3];  
      tri=tricheck(node(osctr,:),tri)
      osctr = osctr(tri)
    else
      osctr = sctr;
    end
    tan_elem = [tan_elem; osctr]
    tan_elem_crk = [tan_elem_crk; elem_crk(tangent_elem(i),:)]
  end
end
