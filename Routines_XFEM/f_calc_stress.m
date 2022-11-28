function [sigma] = f_calc_stress(pt,e,u,C,type_elem,enrich_node,elem_crk,xVertex,xyTip,crack_nodes,k,varargin)

global node element elemType E nu Cm1

if nargin == 12
  QT = varargin{1}
else
  QT = eye(2);
end

B = [];
for k = 1:size(xCr,2)
    B = [B xfemBmat(pt,e,type_elem,enrich_node(:,k),elem_crk,xVertex,xyTip,crack_nodes,k)];
end

leB = size(B,2);

U = [];
for k = 1:size(xCr,2)
    U     = [U; element_disp(e,pos(:,k),enrich_node(:,k),u,k)];
end

epsilon = B*U ;
sigma   = C*epsilon;

voit2ind    = [1 3;3 2];
sigma   = QT*sigma(voit2ind)*QT';
