function [enrich_node, n_red] = f_phantomNode( crack_nodes, elem_crk, split_elem, tip_elem, vertex_elem, enrich_node )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Wed Nov 10 11:43:15 NZDT 2021

global node element numnode numelem elemType
global orig_nn
% crack_nodes has columns: detected node, elem_crack coord 1 x, elem_crack coord 1 y, elem_crack coord 2 x, elem_crack coord 2 y,
orig_nn = numnode;
n_red = 0;

for i = 1:length(crack_nodes)
  n1 = crack_nodes(i,1);
  n2 = numnode+1;
  node(n2,:) = node(n1,:);
  enrich = enrich_node(n1);
  if enrich == 2
    enrich_node(n2) = 0; % this could change depending on implementation 
    n_red = n_red + 1;
  else
    enrich_node(n2) = 0;
  end

  % now we deal with elements that have corner nodes. 
  % find list of elements that support this node
  [els,n_ind] = find(element==n1);
  for j = 1:length(els)
    e = els(j);
    test = [ ismember(e,split_elem),ismember(e,tip_elem), ismember(e,vertex_elem) ] ;
    if ~any( test )  
      xCr = [elem_crk(e,1:2); elem_crk(e,3:4)];
      me = mean(node(element(e,:),:));
      d = signed_distance(xCr,me,0);  
      if d > 0  
         np = n_ind(j);
         element(e,np) = n2;
      end
    end
  end
end

numnode = size(node,1);


end
