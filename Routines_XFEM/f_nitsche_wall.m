function [p,pg,pns,nns] = f_nitsche_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node)
global node element epsilon

sctr = element(iel,:) ;
vv = node(sctr,:);
[phi] = dista(iel,elem_crk) ;
[cutEdge, nnodes,pns,nns] = f_edgedetect(nnode, corner,  phi) ;
nEdge = length(cutEdge);
if nEdge==1 % then one of the nodes must be a crack_node
  crack_n = intersect(crack_node,sctr);
  if isempty(crack_n)
    keyboard
    error("only one edge but no crack node in elem ,",num2str(iel),", when evaluating crack lips")
  end
  crack_c = find(sctr==crack_n);
  nnodes = [ nnodes; nnodes(crack_c,:) ] ;
end

p = [nnodes(end-1,:);nnodes(end,:)];
pg = [elem_crk(iel,1:2);elem_crk(iel,3:4)];

%if ismember(iel,vertex_elem)
  %tip = xVertex(iel,:);
  %ntip = f_naturalpoint(tip,vv,20,1e-6);
  %p = [p(1,:) ; ntip ; p(2,:) ];
  %pg = [pg(1,:);tip;pg(2,:) ];
%end

end
