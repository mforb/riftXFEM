function [sigma] = f_calc_melange_stress(pt,e,u,C,type_elem,enrich_node,elem_crk,xVertex,xyTip,crack_nodes,pos,k,varargin)

global node element elemType E nu Cm1

if nargin == 13
  QT = varargin{1};
else
  QT = eye(2);
end

[B, dJ] = xfemBmel(pt,e,type_elem,enrich_node(:,k),elem_crk,xVertex,xyTip,crack_nodes,kn,QT,mT,mE);


[A,~,~,~,~] = f_enrich_assembly(e,pos,type_elem,elem_crk,enrich_node);
U = u(A);

epsilon = B*U ;
sigma   = C*epsilon;

voit2ind    = [1 3;3 2];
%sigma   = QT*sigma(voit2ind)*QT';
sigma   = sigma(voit2ind);
