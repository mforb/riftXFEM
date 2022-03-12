function [p] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,crack_node)
global node, epsilon

sctr = element(iel,:) ;
vv = node(sctr,:);
[phi] = dista(iel,elem_crk) ;
if ismember(iel, tip_elem) % for now we wont deal with this element
  tip = xTip(iel,:);
  ntip = f_naturalpoint(tip,vv,20,1e-6);
  psi = f_dista2(iel,elem_crk,tip);
  [cutEdge,nnodes] = f_edgedetect(nnode, corner,  phi, psi) ;
  p = [nnodes(end,:); ntip ];
else 
  [cutEdge, nnodes] = f_edgedetect(nnode, corner,  phi) ;
  nEdge = length(cutEdge);
  if nEdge==1 % then one of the nodes must be a crack_node
    crack_n = intersect(crack_node,sctr);
    if isempty(crack_n)
      error("only one edge but no crack node in elem ,",num2str(iel),", when evaluating crack lips")
    end
    crack_c = find(sctr==crack_n);
    nnodes = [ nnodes; nnodes(crack_c,:) ] ;
  end

  p = [nnodes(end-1,:);nnodes(end,:)];

  if ismember(iel,vertex_elem) % we are jsut going to draw a ligne across the two intersections
    tip = xVertex(iel,:);
    ntip = f_naturalpoint(tip,vv,20,1e-6);
    %p = [p(1,:) ; ntip ; p(2,:) ];
  end
end
