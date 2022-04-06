function [node_iso] = f_tangent_iso_node( tangent_elem,crack_node )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Wed Nov 10 11:43:15 NZDT 2021

global node element 
node_iso = [];
for i = 1:length(tangent_elem)
  % test if adjacent to other tangent_elem
  sctr = element(tangent_elem(i),:);
  nsctr = setdiff(sctr,crack_node);
  for n = 1:length(nsctr)
    [r,c] = find(element(tangent_elem,:)==nsctr(n));
    if length(r) > 1  
      node_iso = [node_iso, nsctr(n) ] ;
    end
  end
end
node_iso = unique(node_iso);
